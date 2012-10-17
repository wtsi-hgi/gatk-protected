import StingConversionUtils as util
import argparse
import random
import math
import numpy

def getCorrelatedTestGenotypes(seed,nSamples,nVariants,distBetweenVar,pMissing,nPreviousPredictors = 3):
 """ this is not an attempt to simulate any kind of natural LD structure, just a means
     of generating genotypes that are correlated so that there can be positive controls
     for testing. One such divergence is that we don't care about the distance for the genotypes.
 """
 random.seed(seed)
 samID = list(map(lambda u: "sample_%d" % (u),range(nSamples)))
 plinkSamples = list(map(lambda x: util.Plink.Sample("%s\t%s\t%s\t%s\t%s\t%s" % (x,x,"0","0",-9,-9)),samID))
 pos = 10000
 chr = "1"
 genotypesByVariant = list()
 gbvIdx = 0
 for varNum in range(nVariants):
  varID = "var_%s_%d" % (chr,pos)
  variant = util.Plink.Variant("%s %s %s %d %s %s" % (chr,varID,"-1",pos,"0","1"))
  # generate the frequency OFFSET for logistic
  alpha = random.gauss(0,0.5)
  # going back nPreviousPredictors, generate beta values
  beta = []
  for t in range(min(nPreviousPredictors,len(genotypesByVariant))):
   # generate a beta value
   beta.append(random.gauss(0,2))
  # generate per-sample genotypes
  genotypes = list()
  for idx in range(len(plinkSamples)):
   sample = plinkSamples[idx]
   logitsum = alpha
   for t in range(min(nPreviousPredictors,len(genotypesByVariant))):
    logitsum += genotypesByVariant[gbvIdx-1-t][idx].getDosage()*beta[t]
   freq = 1/(1+math.exp(logitsum))
   r = random.random()
   if ( r < freq**2 ):
    genotype = 3
   elif ( r < freq**2 + 2*freq*(1-freq) ):
    genotype = 2
   else:
    genotype = 0
   genotypes.append(util.Plink.Genotype(variant,sample,genotype))
  genotypesByVariant.append(genotypes)
  gbvIdx += 1
 # now at the end we generate some missingness
 if ( pMissing > 0 ):
  for genoList in genotypesByVariant:
   for geno in genoList:
    r = random.random()
    if ( r < pMissing ):
     geno.setType(util.Plink.Genotype.Type.NO_CALL)
     assert geno.getDosage() == -1
 return genotypesByVariant

def getCorrelatedDosages(seed,nSamples,nVariants,distanceBetweenVar=100,pMissing=-1,nPreviousPredictors=3):
 genotypes = getCorrelatedTestGenotypes(seed,nSamples,nVariants,distanceBetweenVar,pMissing,nPreviousPredictors)
 dosages = list()
 for geno in genotypes:
  innerDosages = list()
  missing = list()
  sumDosage = 0.0
  idx = 0
  for g in geno:
   d = float(g.getDosage())
   innerDosages.append(d)
   if ( d < 0 ):
    missing.append(idx)
   else:
    sumDosage += d
   idx += 1
  freq = sumDosage/(2*(len(innerDosages)-len(missing)))
  var = geno[0].variant
  var.frequency = freq
  dosageTriple = (var,numpy.array(innerDosages),missing)
  dosages.append(dosageTriple)
 return dosages


def getIndependentTestGenotypes(seed,nSamples,nVariants,distBetweenVar,pMissing):
 random.seed(seed)
 samID = list(map(lambda u: "sample_%d" % (u),range(nSamples)))
 plinkSamples = list(map(lambda x: util.Plink.Sample("%s\t%s\t%s\t%s\t%s\t%s" % (x,x,"0","0",-9,-9)),samID))
 pos = 10000
 chr = "1"
 genotypesByVariant = list()
 for varNum in range(nVariants):
  varID = "var_%s_%d" % (chr,pos)
  variant = util.Plink.Variant("%s %s %s %d %s %s" % (chr,varID,"-1",pos,"0","1"))
  # generate the frequency - a beta distribution
  freq_num = -math.log(random.random())
  freq_denom = -math.log(random.random())
  for t in range(4):
   freq_denom = -math.log(random.random())
  freq = freq_num/(freq_num+freq_denom)
  # map a genotype over all the samples
  genotypes = list()
  for sample in plinkSamples:
   r = random.random()
   if ( r < freq**2 ):
    genotype = 3
   elif ( r < freq**2 + 2*freq*(1.0-freq) ):
    genotype = 2
   else:
    genotype = 0
   if ( pMissing > 0 ):
    r = random.random()
    if ( r < pMissing ):
     genotype = 1
   genotypes.append(util.Plink.Genotype(variant,sample,genotype))
  genotypesByVariant.append(genotypes)
  pos += distBetweenVar
 return genotypesByVariant

def getIndependentDosages(seed,nSamples,nVariants,distance = 100, pMissing = -1):
 genotypes = getIndependentTestGenotypes(seed,nSamples,nVariants,distance,pMissing)
 dosages = list()
 for geno in genotypes:
  innerDosages = list()
  missing = list()
  sumDosage = 0.0
  idx = 0
  for g in geno:
   d = float(g.getDosage())
   innerDosages.append(d)
   if ( d < 0 ):
    missing.append(idx)
   else:
    sumDosage += d
   idx += 1
  freq = sumDosage/(2*(len(innerDosages)-len(missing)))
  var = geno[0].variant
  var.frequency = freq
  dosageTriple = (var,numpy.array(innerDosages),missing)
  dosages.append(dosageTriple)
 return dosages

def getStandardIntervals():
 return [("1",11000,12000),("1",13000,14000)]

def getArgs(localBP = 0,bedBase = "./test/1000G_subset.chr20.79234"):
 parser = argparse.ArgumentParser(description='Parse out arguments for the generalized burden test')
 parser.add_argument("-bedBase",action='store',default=bedBase,dest="bedBase",
                required=False,help="The base of the bed (so 'xxxx' if file named 'xxxx.bed'). Requires a bim and fam with the same base name.")
 parser.add_argument("-intervals",action='store',default=None,dest="intervals",required=False,help="A file listing intervals over which the GRM should be calculated")
 parser.add_argument("-localCorrectionBP",action='store',default=localBP,type=int,dest="localCorrectionBP",required=False,
          help="Use variants outside of the given intervals but within this many BP to correct for LD structure.")
 parser.add_argument("-globalCorrection",action='store',default=None,dest="globalCorrection",required=False,
          help="Include this list of variant IDs in LD correction")
 parser.add_argument("-minFrequency",action='store',default=1e-7,type=float,dest="minFrequency",required=False,
          help="The minimum frequency of a variant to use in the GRM calculation.")
 parser.add_argument("-maxFrequency",action='store',default=1.0+1e-7,type=float,dest="maxFrequency",required=False,
          help="The maximum frequency of a variant to use in the GRM calculation.")
 parser.add_argument("-minFrequencyCorrect",action='store',default=1e-7,type=float,dest="minFrequencyCorrect",required=False,
          help="The minimum frequency of a variant to use for local/global LD correction")
 parser.add_argument("-maxFrequencyCorrect",action='store',default=1.0-1e-7,type=float,dest="maxFrequencyCorrect",required=False,
          help="The maximum frequency of a variant to use for local/global LD correction")
 parser.add_argument("-naCorrect",action='store_true',dest="naCorrect",required=False,
          help="Replace missing genotype values with the mean dosage (2 x p) rather than ignoring them")
 return parser.parse_args()
