import StingConversionUtils as utils
import PlinkBedReader as pb
import argparse 

class VCFHeaderConstants:
 fileFormat = "##fileformat=VCFv4.1"
 sourceFormat = "##bedToVCF.py -vcf %s -bedBase %s\n"
 alleleCountFormat = '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">'
 alleleFrequencyFormat = '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">'
 alleleNumberFormat = '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">'

 standardFields = ["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]

class VCFBodyConstants:
 gtFieldFormat = "GT"
 filterField = "."
 qualField = "."
 acKey = "AC"
 anKey = "AN"
 afKey = "AF"

 def getStandardFields(chr,pos,id,ref,alt,info):
  return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chr,pos,id,ref,alt,VCFBodyConstants.qualField,VCFBodyConstants.filterField,info,VCFBodyConstants.gtFieldFormat)

def parseArgs():
 parser = argparse.ArgumentParser(description='Parse out arguments for the generalized burden test')
 parser.add_argument("-vcf",action='store',default=None,dest="vcf",
               required=True,help="The variants and genotypes on which to run the test, in VCF format")
 parser.add_argument("-bedBase",action='store',default=None,dest="bedBase",
               required=True,help="The base of the bed (so 'xxxx' if file named 'xxxx.bed'. Requires a bim and fam with the same base name.")
 return parser.parse_args()

def writeHeader(outVCF,inBed,args):
 outVCF.write("%s\n" % (VCFHeaderConstants.fileFormat))
 outVCF.write(VCFHeaderConstants.sourceFormat % (args.vcf,args.bedBase))
 outVCF.write("%s\n%s\n%s\n" % (VCFHeaderConstants.alleleCountFormat,VCFHeaderConstants.alleleFrequencyFormat,VCFHeaderConstants.alleleNumberFormat))
 fields = VCFHeaderConstants.standardFields
 fields.extend(map(lambda x: x.individual_id,pb.getSamples(inBed)))
 outVCF.write("#%s\n" % ("\t".join(fields)))

def writeBody(outVCF,inBed):
 curVariant = None
 formattedGenotypes = list()
 alleleCount = 0
 alleleNumber = 0
 for genotype in inBed:
  variant = genotype.variant
  if ( curVariant != None and variant.chr == curVariant.chr and variant.pos == curVariant.pos ):
   # sanity check
   assert variant.id == curVariant.id
   alleleCount += genotype.getDosage() 
   alleleNumber += 0 if genotype.isNoCall() else 2
   formattedGenotypes.append(utils.plinkGenotypeToVCFString(genotype))
  else:
   # dump the logged variant
   if ( curVariant != None ):
    info = "%s;%s;%s" % ("%s=%d" %   (VCFBodyConstants.acKey,alleleCount),
                         "%s=%.3f" % (VCFBodyConstants.afKey,float(alleleCount)/alleleNumber),
                         "%s=%d" %   (VCFBodyConstants.anKey,alleleNumber))
    outVCF.write("%s\t%s\n" % 
       (VCFBodyConstants.getStandardFields(variant.chr,variant.pos,variant.id,variant.ref,variant.alt,info),
        "\t".join(formattedGenotypes)))
   # reset data
   formattedGenotypes = list()
   alleleCount = 0
   alleleNumber = 0
   # process the new variant
   curVariant = variant
   alleleCount += genotype.getDosage()
   alleleNumber += 0 if genotype.isNoCall() else 2
   formattedGenotypes.append(utils.plinkGenotypeToVCFString(genotype))
 if ( curVariant != None ):
  info = "%s;%s;%s" % ("%s=%d" %   (VCFBodyConstants.acKey,alleleCount),
                       "%s=%.3f" % (VCFBodyConstants.afKey,float(alleleCount)/alleleNumber),
                       "%s=%d" %   (VCFBodyConstants.anKey,alleleNumber))
  outVCF.write("%s\t%s\n" % (VCFBodyConstants.getStandardFields(variant.chr,variant.pos,variant.id,variant.ref,variant.alt,info),"\t".join(formattedGenotypes)))

def runMain():
 args = parseArgs()
 vcfOutput = open(args.vcf,'w')
 bedReader = pb.getReader(args.bedBase)
 writeHeader(vcfOutput,bedReader,args)
 writeBody(vcfOutput,bedReader)

if __name__ == '__main__':
 runMain()
