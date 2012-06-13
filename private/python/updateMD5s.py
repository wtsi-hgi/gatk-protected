from optparse import OptionParser
import itertools
import os

def main():
    global OPTIONS
    usage = """usage: %prog [options] mode md5mismatches ... sourceDir
    Update the MD5s across multiple files.
    
    The first argument is the database of MD5 updates produced automatically by the
    integrationtests, which lives in integrationtest/md5mismatches.txt.
    
    The second (and more) argument(s) is the path to the source directory containing the java / scala
    code with md5s to update.  This program crawls that directory looking for .java
    and .scala files to update, id's those with MD5s in the md5mismatches.txt files and
    updates these files.
    
    It is highly recommended that you commit your code before updating the MD5s to ensure
    that if things go wrong with this script that you don't lose any work."""
    
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v", "--verbose", dest="verbose",
                        action='store_true', default=False,
                        help="If provided, verbose progress will be enabled")         

    parser.add_option("-m", "--allowMissingMismatches", dest="allowMissingMismatches",
                        action='store_true', default=False,
                        help="If provided, we will tolerate not having a source file to update for every mismatch")         

    parser.add_option("-n", "--dry", dest="dry",
                        action='store_true', default=False,
                        help="If provided, don't actually update anything")         

    (OPTIONS, args) = parser.parse_args()
    if len(args) < 2:
        parser.error("Requires a md5mismatch db and at least one path to find files to update")
        
    print 'Reading md5 mismatches from', args[0], '...' 
    md5DB = readMD5Mismatches(args[0])
    print '    Read', len(md5DB), 'mismatches' 

    print '\nLooking for files to update...'
    filesToUpdate = []
    for file in findFilesInDirectories(args[1:]):
        found = findMismatchSources(md5DB, file)
        if found:
            filesToUpdate.append(found)
            
    print '\nEnumerating files we will update...'
    for mismatch in md5DB:
        if mismatch.getSources() == [] and not OPTIONS.allowMissingMismatches:
            raise Exception("Didn't find a source file to update for mismatch %s" % str(mismatch))
        print '    Mismatch %s will update file(s) %s' % (mismatch.oldMD5, mismatch.getSources())

    print '\nUpdating files...'
    for mismatch in md5DB:
        print '    Updating file(s) for', mismatch
        updateMismatch(mismatch)

class Mismatch:
    def __init__(self, oldMD5, newMD5, name):
        self.oldMD5 = oldMD5
        self.newMD5 = newMD5
        self.name = name
        self.sources = list()
        
    def addSource(self, file):
        self.sources.append(file)

    def getSources(self):
        return self.sources
        
    def __str__(self):
        return "[Mismatch md5=%s name=%s]" % (self.oldMD5, self.name) 
        
def readMD5Mismatches(file):
    def parseOne(line):
        parts = line.strip().split("\t")
        if len(parts) != 3:
            raise Exception("Unexpected number of elements in md5 file at line: " + line)
        return Mismatch(parts[0], parts[1], parts[2]) 
        
    return [parseOne(line) for line in open(file) if not line.startswith("expected")]
    
def findMismatchSources(db, file):
    txt = open(file).read()
    found = False
    for mismatch in db:
        if txt.find(mismatch.oldMD5) != -1:
            mismatch.addSource(file)
            print '    MD5 update found: expected=%s => new=%s in file %s' % ( mismatch.oldMD5, mismatch.newMD5, file )
            found = True
    return False
    
def findFilesInDirectories(directoriesOrFiles):
    if len(directoriesOrFiles) == 0:
        raise Exception("BUG: no directories provided")

    filesToConsiderUpdating = []
    for top in directoriesOrFiles:
        for root, dirs, files in os.walk(top):
            for name in files:
                if name.endswith(".java") or name.endswith(".scala"):
                    filesToConsiderUpdating.append(os.path.join(root, name))
    return filesToConsiderUpdating

def updateMismatch(mismatch):
    for file in mismatch.getSources():
        print '        => Updating file', file
        fd = open(file)
        txt = fd.read()
        txt = txt.replace(mismatch.oldMD5, mismatch.newMD5)
        fd.close()
        print '        => Updated text'
        if not OPTIONS.dry:
            fd = open(file, "w")
            fd.write(txt)
            fd.close()
            print '      => Updated file on disk'

if __name__ == "__main__":
    main()
