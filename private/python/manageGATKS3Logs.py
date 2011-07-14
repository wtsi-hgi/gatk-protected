import os.path
import sys
from optparse import OptionParser
import subprocess
from itertools import *
import multiprocessing
import time

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

def s3bucket():
    return "s3://" + OPTIONS.bucket

def execS3Command(args, stdout = None):
    """Executes the S3cmd command, putting results into stdout, if provided"""
    executionString = " ".join([OPTIONS.S3CMD] + args)
    if OPTIONS.dryRun:
        print 'DRY-RUN:', executionString
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

# TODO -- get by groups
def getFilesInBucket(args, delete=False):
    lsFile, logFile = args
    logLines = [line.split() for line in open(logFile)]
    alreadyGot = set([parts[1] for parts in logLines if parts[0] == "get"])
    alreadyDel = set([parts[1] for parts in logLines if parts[0] == "del"])
    print 'alreadyGot', len(alreadyGot)
    print 'alreadyDel', len(alreadyDel)
    log = open(logFile, 'a')

    def writeLog(action, skipSet, files):
        for file in files: 
            if file not in skipSet: 
                print >> log, action, file

    def filterExistingFiles(files):
        def alreadyExists(file):
            if file in alreadyGot or OPTIONS.FromScratch: 
                return True
            elif OPTIONS.checkExistsOnDisk:
                destFile = os.path.join(OPTIONS.DIR, file.replace(s3bucket() + "/", ""))
                return os.path.exists(destFile)
            else:
                return False
        return filter(lambda x: not alreadyExists(x), files)

    def processGroup(filesInGroupRaw):
        filesInGroup = filter(lambda x: x != None, list(filesInGroupRaw))
        print '\ngroup:', len(filesInGroup), 'files'
        filesToGet = filterExistingFiles(filesInGroup)
        print 'to get:', len(filesToGet), 'files'
        if filesToGet != []:
            destDir = OPTIONS.DIR
            if OPTIONS.verbose: print 'Getting files', filesToGet, 'to', destDir
            execS3Command(["get", "--force"] + filesToGet + [destDir])
            writeLog('get', alreadyGot, filesToGet)
        if delete:
            filesToDel = [file for file in filesInGroup if file not in alreadyDel]
            if OPTIONS.verbose: print 'Deleting remotes', filesToDel
            execS3Command(["del"] + filesToDel) 
            writeLog('del', alreadyDel, filesToDel)
        return 'Complete'

    print map( processGroup, getFilesFromS3LSByGroup(lsFile) )
            
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
    

    
