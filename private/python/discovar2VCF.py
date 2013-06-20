from optparse import OptionParser
import itertools
import os
import sys
import re

def main():
    global OPTIONS
    usage = """usage: %prog [options] discover.variants ...
    Convert discovar variant files to VCF
    """
    
    parser = OptionParser(usage=usage)

    (OPTIONS, args) = parser.parse_args()
    if len(args) < 1:
        parser.error("Requires at least 1 discovar.variants file")

    print header()
    for file in args:
        print >> sys.stderr, 'Converting ', file, '...'
        for line in open(file):
            discovarLineToVCF(line)

def header():
    return """##fileformat=VCFv4.1
##source=discovar
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"""

# VAR    225 20:10096905 TA         T          -|+       PROB= 1
# VAR    226 20:10096933 G          C          -|+       PROB= 1
# VAR    230 20:10097436 CTTTTCTTT+ C          -|+       PROB= 1 REF=CTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTT;ALT=C
# VAR    231 20:10097437 TTTTC      T          +|-       PROB= 1
def discovarLineToVCF(line):
    parts = line.split()
    if line.strip() != "" and parts[0] == "VAR":
        chr, pos = parts[2].split(":")
        ref, alt = parseRefAlt(parts)
        print "\t".join([chr, pos, ".", ref, alt, "99", "PASS", "."])

pattern = re.compile("REF=(.*);ALT=(.*)")
def parseRefAlt(parts):
    ref, alt = parts[3], parts[4]
    if "+" in ref or "+" in alt:
        fullAlleles = pattern.match(parts[8])
        if fullAlleles == None:
            raise Exception("Couldn't match", parts[8])
        return fullAlleles.group(1), fullAlleles.group(2)
    else:
        return ref, alt

if __name__ == "__main__":
    main()
