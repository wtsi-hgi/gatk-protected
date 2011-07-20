import os.path
import sys
from optparse import OptionParser
import subprocess
from itertools import *
import multiprocessing
import time
import Queue

# 
#
# A script for downloading GATK S3 logs
#
# Structure: first get an ls of the bucket
#   then either get the results of this ls or move (get + delete on server) these files
#   supports upload feature for testing
#
#

MODES = dict()
DEVNULL = open("/dev/null", 'w')

def s3bucket():
    return "s3://" + OPTIONS.bucket

def execS3Command(args, stdout = None):
    """Executes the S3cmd command, putting results into stdout, if provided"""
    executionString = " ".join([OPTIONS.S3CMD] + args)
    if OPTIONS.dryRun:
        if OPTIONS.verbose: print 'DRY-RUN:', executionString
        return
    try:
        retcode = subprocess.call(executionString, shell=True, stdout=stdout)
        if retcode < 0:
            print >>sys.stderr, "Child was terminated by signal", -retcode
    except OSError, e:
        print >>sys.stderr, "Execution failed:", e

def putFilesToBucket(args):
    execS3Command(["put"] + args + [s3bucket()]) 

def lsBucket(args):
    print 'Logging ls to', args[0]
    execS3Command(["ls", s3bucket()], stdout=open(args[0], 'w'))
    if OPTIONS.verbose:
        print 'ls:', args[0]
        for line in open(args[0]): print line,

def getFilesFromS3LSByGroup(file):
    def fileStream():
        for line in open(file):
            yield line.split()[3]
    return grouper(OPTIONS.GROUP_SIZE, fileStream())
            
def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)    

class GetWorker(multiprocessing.Process):
    MAX_BLOCKING_TIME = 60 # seconds
    
    def __init__(self, work_queue, result_queue, delete, alreadyGot, alreadyDel):
        # base class initialization
        multiprocessing.Process.__init__(self)
 
        # job management stuff
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.kill_received = False
        self.delete = delete
        self.alreadyGot = alreadyGot
        self.alreadyDel = alreadyDel

    def filterExistingFiles(self, files):
        def alreadyExists(file):
            if file in self.alreadyGot or OPTIONS.FromScratch: 
                return True
            elif OPTIONS.checkExistsOnDisk:
                destFile = os.path.join(OPTIONS.DIR, file.replace(s3bucket() + "/", ""))
                return os.path.exists(destFile)
            else:
                return False
        return filter(lambda x: not alreadyExists(x), files)

    def getOutStream(self):
        if OPTIONS.verbose:
            return sys.stdout
        else:
            return DEVNULL
        
    def processGroup(self, filesInGroupRaw):
        outstream = self.getOutStream()
        if OPTIONS.debug: print 'process id:', os.getpid()
        filesInGroup = filter(lambda x: x != None, list(filesInGroupRaw))
        if OPTIONS.debug: print '\ngroup:', len(filesInGroup), 'files'
        filesToGet = self.filterExistingFiles(filesInGroup)
        filesToDel = [] # by default we aren't deleting anything
        if OPTIONS.debug: print 'to get:', len(filesToGet), 'files'
        if filesToGet != []:
            destDir = OPTIONS.DIR
            if OPTIONS.debug: print 'Getting files', filesToGet, 'to', destDir
            execS3Command(["get", "--force"] + filesToGet + [destDir], stdout=outstream)
        if self.delete:
            filesToDel = [file for file in filesInGroup if file not in self.alreadyDel]
            if OPTIONS.debug: print 'Deleting remotes', filesToDel
            execS3Command(["del"] + filesToDel, stdout=outstream ) 
        return os.getpid(), filesToGet, filesToDel
 
    def run(self):
        while not self.kill_received and not self.work_queue.empty():
            # get a task
            try:
                filesInGroup = self.work_queue.get(True, self.MAX_BLOCKING_TIME)
                result = self.processGroup(filesInGroup)
            except Queue.Empty:
                print 'Stopping run()'
                break
 
            # the actual processing
            self.result_queue.put(result)

def getFilesInBucket(args, delete=False):
    lsFile, logFile = args
    
    # load up the sets of already downloaded and deleted
    if os.path.exists(logFile):
        #
        logLines = [line.split() for line in open(logFile)]
    else:
        logLines = []
    alreadyGot = set([parts[1] for parts in logLines if parts[0] == "get"])
    alreadyDel = set([parts[1] for parts in logLines if parts[0] == "del"])
    print 'alreadyGot', len(alreadyGot)
    print 'alreadyDel', len(alreadyDel)

    # Logging progress to file        
    log = open(logFile, 'a')
    def writeLog(action, skipSet, files):
        if not OPTIONS.dryRun:
            for file in files: 
                if file not in skipSet: 
                    print >> log, action, file
          
    # parallel processing
    # load up work queue
    work_queue = multiprocessing.Queue()
    nFilesToProcess = 0
    jobs = list(getFilesFromS3LSByGroup(lsFile))
    for filesGroup in jobs:
        work_queue.put(filesGroup)
        nFilesToProcess += len(filter(lambda x: x != None, filesGroup))
    print 'Number of work units', len(jobs)
 
    # create a queue to pass to workers to store the results
    result_queue = multiprocessing.Queue()
 
    # spawn workers
    for i in range(OPTIONS.N_PARALLEL_PROCESSES):
        worker = GetWorker(work_queue, result_queue, delete, alreadyGot, alreadyDel)
        if OPTIONS.debug: print 'Starting worker', worker
        worker.start()
 
    # collect the results off the queue
    results = []
    nGot, nDel = 0, 0
    while len(results) < len(jobs):
        pid, filesGot, filesDel = result_queue.get()
        results.append([filesGot, filesDel])
        if OPTIONS.debug:
            print '\nPID          ', pid
        writeLog('get', alreadyGot, filesGot)
        writeLog('del', alreadyDel, filesDel)
        nGot = nGot + len(filesGot)
        nDel = nDel + len(filesDel)
        print 'Total %d got, %d deleted of %d overall. This work unit got %d files, deleted %d files in group' % (nGot, nDel, nFilesToProcess, len(filesGot), len(filesDel))

    print '\nDownloading complete'
    print 'No. files downloaded   :', nGot
    print 'No. files deleted at S3:', nDel
        
            
# Create the mode map
MODES["upload"] = putFilesToBucket
MODES["ls"] = lsBucket
MODES["get"] = getFilesInBucket
MODES["move"] = lambda x: getFilesInBucket(x, delete=True)

if __name__ == "__main__":
    usage = "usage: %prog [options] mode args"
    parser = OptionParser(usage=usage)
    parser.add_option("-b", "--bucket", dest="bucket",
                        type='string', default="GATK_Run_Reports",
                        help="If true, will use")
    parser.add_option("-d", "--dir", dest="DIR",
                        type='string', default='.',
                        help="Path to write local logs to")
    parser.add_option("-s", "--s3cmd", dest="S3CMD",
                        type='string', default="/Users/depristo/Desktop/broadLocal/s3/s3cmd-1.0.1/s3cmd",
                        help="Path to s3cmd executable")
    parser.add_option("-g", "--groupSize", dest="GROUP_SIZE",
                        type='int', default=100,
                        help="Number of elements to get at the same time")
    parser.add_option("-r", "--fromScratch", dest="FromScratch",
                        action='store_true', default=False,
                        help="If provided, we will redownload files already present locally")
    parser.add_option("", "--checkExistsOnDisk", dest="checkExistsOnDisk",
                        action='store_true', default=False,
                        help="If provided, we will check the file system for a file before we download it")
    parser.add_option("", "--dryRun", dest="dryRun",
                        action='store_true', default=False,
                        help="If provided, we will not actually execute any s3 commands")
    parser.add_option("-v", "--verbose", dest="verbose",
                        action='store_true', default=False,
                        help="If provided, will print verbose info about activities")
    parser.add_option("", "--debug", dest="debug",
                        action='store_true', default=False,
                        help="If provided, will print debugging info")
    parser.add_option("-p", "--parallel", dest="N_PARALLEL_PROCESSES",
                        type='int', default=2,
                        help="Number of parallel gets to execute at the same time")
                        
    (OPTIONS, args) = parser.parse_args()
    if len(args) < 1:
        parser.error("incorrect number of arguments")

    if not os.path.exists(OPTIONS.DIR):
        os.makedirs(OPTIONS.DIR)

    MODES[args[0]](args[1:])
    

    
