import os.path
import sys
from optparse import OptionParser
import subprocess
from itertools import *

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
    def filterExistingFiles(files):
        def alreadyExists(file):
            destFile = os.path.join(OPTIONS.DIR, file.replace(s3bucket() + "/", ""))
            return os.path.exists(destFile) or OPTIONS.FromScratch
        return filter(lambda x: not alreadyExists(x), files)
            
    for filesInGroupRaw in getFilesFromS3LSByGroup(args[0]):
        filesInGroup = filter(lambda x: x != None, list(filesInGroupRaw))
        print '\ngroup:', len(filesInGroup), files
        filesToGet = filterExistingFiles(filesInGroup)
        print '\nto get:', len(filesToGet), files
        if filesToGet != []:
            destDir = OPTIONS.DIR
            print 'Getting files', filesToGet, 'to', destDir
            execS3Command(["get", "--force"] + filesToGet + [destDir]) 
        if delete:
            print 'Deleting remotes', filesInGroup
            execS3Command(["del"] + filesInGroup) 
            
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
    parser.add_option("", "--dryRun", dest="dryRun",
                        action='store_true', default=False,
                        help="If provided, we will not actually execute any s3 commands")
                        
    (OPTIONS, args) = parser.parse_args()
    if len(args) < 1:
        parser.error("incorrect number of arguments")

    if not os.path.exists(OPTIONS.DIR):
        os.makedirs(OPTIONS.DIR)

    MODES[args[0]](args[1:])
    

    
