import sys
import logging

import LicenseUtils

logging.basicConfig(format="[CheckLicense]: %(message)s", level = logging.INFO)

exitStatus = 0
for filename in sys.stdin.readlines():
    filename = filename.strip()
    if LicenseUtils.isSourceFile(filename):
        license = LicenseUtils.extractLicense(LicenseUtils.getAppropriateLicense(filename))
        sourceFile = open(filename)
        source = LicenseUtils.extractLicenseFromSource(sourceFile)
        if license != source:
            logging.info("Wrong license for file " + filename)
            exitStatus = 1
exit(exitStatus)