#!/usr/bin/python
import sys

def lineIsNotCommentedOut(line):
    return not line == "" and not line.startswith(" *") and not line.startswith("*") and not line.startswith("/*") and not line.startswith("//") and not line == "\n" and not line == "\r\n"

licenseFile = open(sys.argv[2])
sys.stdout.write("/*\n")
for line in licenseFile.readlines():
    sys.stdout.write("* " + line)

skipLicense = True
sourceFile = open(sys.argv[1])
for line in sourceFile.readlines():
    strippedLine = line.strip()
    if skipLicense and strippedLine.startswith("package"):
        sys.stdout.write("*/\n\n")
        skipLicense = False
    elif skipLicense and lineIsNotCommentedOut(strippedLine):
        sys.stderr.write("***ERROR***: Couldn't find package information for this file\n")
        exit(1)

    if not skipLicense:
        sys.stdout.write(line)

