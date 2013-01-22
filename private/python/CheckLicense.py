import sys
import re
import logging

def blankLine(line):
    return re.match("^$", line) is not None

def lineIsCommentedOut(line):
    return line.startswith("//") or line.startswith("/*")

def extractLicenseFromSource(file):
    extractedLicense = ""
    inCommentBlock = 0
    for line in file.readlines():
        strippedLine = line.strip()
        if strippedLine.startswith("package"):
            return extractedLicense # found package line, return the license.

        if blankLine(strippedLine):
            continue                # just skip blank lines

        if strippedLine.startswith("/*"):
            inCommentBlock += 1     # mark start of a comment block

        if inCommentBlock or lineIsCommentedOut(strippedLine):
            extractedLicense += strippedLine.strip(" * ").strip(" // ")
        else:
            break

        if strippedLine.endswith("*/"):
            inCommentBlock -= 1     # mark end of a comment block

    return False

def extractLicense(file):
    license = ""
    for line in file.readlines():
        license += line.strip()
    return license

licenseFile = open(sys.argv[1])
license = extractLicense(licenseFile)
exitStatus = 0
for i in range(2,len(sys.argv)):
    sourceFile = open(sys.argv[i])
    source = extractLicenseFromSource(sourceFile)
    if license != source:
        logging.error("Wrong license for file " + sys.argv[i])
        exitStatus = 1

exit(exitStatus)