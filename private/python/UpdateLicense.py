#!/usr/bin/python
import sys
import logging

def lineIsNotCommentedOut(line):
    return not line == "" and not line.startswith(" *") and not line.startswith("*") and not line.startswith("/*") and not line.startswith("//") and not line == "\n" and not line == "\r\n"

def isSourceFile(filename):
    return str.endswith(filename, ".java") or str.endswith(filename, ".scala")

def getAppropriateLicense(filename):
    licenseFileName = ""
    if filename.startswith("private"):
        licenseFileName = "licensing/private_license.txt"
    elif filename.startswith("protected"):
        licenseFileName = "licensing/protected_license.txt"
    elif filename.startswith("public"):
        licenseFileName = "licensing/public_license.txt"
    else:
        logging.error(filename + " is not in public, private or protected\n")
        exit(1)
    return licenseFileName

logging.basicConfig(format="%(levelname)s: %(message)s", level = logging.INFO)
for filename in sys.stdin.readlines():
    filename = filename.strip()

    if isSourceFile(filename):
        licenseFile = open(getAppropriateLicense(filename))
        sourceFile = open(filename)
        updatedSource = "/*\n"
        for line in licenseFile.readlines():
            updatedSource += "* " + line

        skipLicense = True
        for line in sourceFile.readlines():
            strippedLine = line.strip()
            if skipLicense and strippedLine.startswith("package"):
                updatedSource += "*/\n\n"
                skipLicense = False
            elif skipLicense and lineIsNotCommentedOut(strippedLine):
                logging.error(filename + " is missing package information")
                exit(2)

            if not skipLicense:
                updatedSource += line

        sourceFile.close()
        sourceFile = open(filename, "w")
        sourceFile.write(updatedSource)
        sourceFile.close()

        logging.info(filename + " [license successfully added]")
    else:
        logging.info(filename + " doesn't require a license")
