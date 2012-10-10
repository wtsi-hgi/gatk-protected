import calculateGRM
import PlinkBedReader as PlinkReader
import argparse
import numpy
import math
import linear_old as linear

def testArgs(testNo):
 """ Returns the arguments for the test args as a parser.
     @testNo - the test number.
 """
 if ( testNo == 1 ):
  parser = argparse.ArgumentParser(description='Parse out arguments for the generalized burden test')
  parser.add_argument("-bedBase",action='store',default="./test/1000G_subset.chr20.79234",dest="bedBase",
                required=False,help="The base of the bed (so 'xxxx' if file named 'xxxx.bed'). Requires a bim and fam with the same base name.")
  parser.add_argument("-intervals",action='store',default=None,dest="intervals",required=False,help="A file listing intervals over which the GRM should be calculated")
  parser.add_argument("-localCorrectionBP",action='store',default=0,type=int,dest="localCorrectionBP",required=False,
          help="Use variants outside of the given intervals but within this many BP to correct for LD structure.")
  parser.add_argument("-globalCorrection",action='store',default=None,dest="globalCorrection",required=False,
          help="Include this list of variant IDs in LD correction")
  parser.add_argument("-minFrequency",action='store',default=1e-7,type=float,dest="minFrequency",required=False,
          help="The minimum frequency of a variant to use in the GRM calculation.")
  parser.add_argument("-maxFrequency",action='store',default=1.0-1e-7,type=float,dest="maxFrequency",required=False,
          help="The maximum frequency of a variant to use in the GRM calculation.")
  parser.add_argument("-minFrequencyCorrect",action='store',default=1e-7,type=float,dest="minFrequencyCorrect",required=False,
          help="The minimum frequency of a variant to use for local/global LD correction")
  parser.add_argument("-maxFrequencyCorrect",action='store',default=1.0-1e-7,type=float,dest="maxFrequencyCorrect",required=False,
          help="The maximum frequency of a variant to use for local/global LD correction")
  parser.add_argument("-naCorrect",action='store_true',dest="naCorrect",required=False,
          help="Replace missing genotype values with the mean dosage (2 x p) rather than ignoring them")
 if ( testNo == 2 ):
  parser = argparse.ArgumentParser(description='Parse out arguments for the generalized burden test')
  parser.add_argument("-bedBase",action='store',default="./test/1000G_subset.chr20",dest="bedBase",
                required=False,help="The base of the bed (so 'xxxx' if file named 'xxxx.bed'). Requires a bim and fam with the same base name.")
  parser.add_argument("-intervals",action='store',default="./test/1000G_subset.chr20.test2.txt",dest="intervals",required=False,help="A file listing intervals over which the GRM should be calculated")
  parser.add_argument("-localCorrectionBP",action='store',default=1000,type=int,dest="localCorrectionBP",required=False,
          help="Use variants outside of the given intervals but within this many BP to correct for LD structure.")
  parser.add_argument("-globalCorrection",action='store',default=None,dest="globalCorrection",required=False,
          help="Include this list of variant IDs in LD correction")
  parser.add_argument("-minFrequency",action='store',default=1e-7,type=float,dest="minFrequency",required=False,
          help="The minimum frequency of a variant to use in the GRM calculation.")
  parser.add_argument("-maxFrequency",action='store',default=1.0-1e-7,type=float,dest="maxFrequency",required=False,
          help="The maximum frequency of a variant to use in the GRM calculation.")
  parser.add_argument("-minFrequencyCorrect",action='store',default=1e-7,type=float,dest="minFrequencyCorrect",required=False,
          help="The minimum frequency of a variant to use for local/global LD correction")
  parser.add_argument("-maxFrequencyCorrect",action='store',default=1.0-1e-7,type=float,dest="maxFrequencyCorrect",required=False,
          help="The maximum frequency of a variant to use for local/global LD correction")
  parser.add_argument("-naCorrect",action='store_true',dest="naCorrect",required=False,
          help="Replace missing genotype values with the mean dosage (2 x p) rather than ignoring them")
 return parser.parse_args()



def runTest():
 ## test 1 -- push a single variant through and test the returned values
 test1arg = testArgs(1)
 reader = PlinkReader.getReader(test1arg.bedBase)
 assert reader.snpMajor()
 assert reader.numGenotypesPerMajor == 12
 genotypeDosages = calculateGRM.getNextVariant(reader,test1arg)
 assert calculateGRM.getNextVariant(reader,test1arg) == None 
 assert str(genotypeDosages) == "(Variant: (20,79234), array([ 0.,  0.,  2.,  0.,  2.,  2.,  1.,  0.,  1.,  1.,  1.,  1.]), [])"
 dosages = dict()
 regressor = list()
 variants = list()
 missing = list()
 calculateGRM.accumulateVariant(genotypeDosages,regressor,variants,[("20",1,66000000)],test1arg.localCorrectionBP,test1arg,missing,list())
 assert str(variants) == "[(Variant: (20,79234), array([ 0.,  0.,  2.,  0.,  2.,  2.,  1.,  0.,  1.,  1.,  1.,  1.]), [])]" 
 assert regressor == [] 
 completed = calculateGRM.findCompletedVariants([("20",1,66000000)],regressor,variants,test1arg.localCorrectionBP)
 assert str(completed) == "[((Variant: (20,79234), array([ 0.,  0.,  2.,  0.,  2.,  2.,  1.,  0.,  1.,  1.,  1.,  1.]), []), [])]" 
 ## so that's the variant, the dosages, the index of missing values, and no regressor variants
 dosCorrect = calculateGRM.getCorrectedDosagesNoRegression(completed[0][0],completed[0][1])
 assert str(numpy.sum(dosCorrect[1])) == "1.11022302463e-15"
 assert str(numpy.prod(dosCorrect[1])) == "0.000240749030421"
 assert str(numpy.prod(1+dosCorrect[1])) == "0.234156983841"
 dosages[dosCorrect[0]]=dosCorrect[1]
 dosMatrix = numpy.matrix([dosages[t] for t in dosages])
 d1 = calculateGRM.calcDistance(2,3,dosMatrix,dict(),test1arg)
 assert str(d1) == "(1, matrix([[-2.]]))"
 d2 = calculateGRM.calcDistance(1,5,dosMatrix,dict(),test1arg)
 assert str(d2) == "(1, matrix([[-2.]]))"
 d3 = calculateGRM.calcDistance(4,4,dosMatrix,dict(),test1arg)
 assert str(d3) == "(1, matrix([[ 2.36363636]]))"

 ## test 2 -- get a chunk of variants with correction variants
 test2arg = testArgs(2)
 reader = PlinkReader.getReader(test2arg.bedBase)
 assert reader.snpMajor()
 assert reader.numGenotypesPerMajor == 12
 # pull in the variants
 dependentGenotypeIntervals = calculateGRM.getIntervals(test2arg.intervals)
 assert dependentGenotypeIntervals == [('20', 11000000, 11100000)] 
 regressorVariants = list()
 variantsToCorrect = list()
 dosages = dict()
 individualMissingCalculate = dict()
 individualMissingCorrect = dict()
 genotypes = calculateGRM.getNextVariant(reader,test2arg)
 # this is a 3-ple: (PlinkVariant,List[PlinkGenotype],List[Int])
 siteNo = 0
 completedVariants = list()
 while ( genotypes != None and len(dependentGenotypeIntervals) > 0):
  if ( siteNo % 10000 == 1 ):
   print("read : "+str(siteNo) + " at : " + str(genotypes[0].pos) + " accum : "+str(len(regressorVariants))+ " compute : "+str(len(variantsToCorrect)))
  calculateGRM.accumulateVariant(genotypes,regressorVariants,variantsToCorrect,dependentGenotypeIntervals,test2arg.localCorrectionBP,test2arg,individualMissingCalculate,individualMissingCorrect)
  doneIntervals = calculateGRM.removeStaleIntervals(genotypes[0],regressorVariants,dependentGenotypeIntervals,test2arg.localCorrectionBP)
  completedVariants = calculateGRM.findCompletedVariants(doneIntervals,regressorVariants,variantsToCorrect,test2arg.localCorrectionBP)
  calculateGRM.removeStaleVariants(doneIntervals,dependentGenotypeIntervals,regressorVariants,variantsToCorrect,test2arg.localCorrectionBP)
  if ( len(completedVariants) > 0 ):
   break
  genotypes = calculateGRM.getNextVariant(reader,test2arg)
  siteNo += 1
 assert str(completedVariants[0][0]) == '(Variant: (20,11000515), array([ 0.,  0., -1.,  0.,  2.,  1.,  2.,  2., -1., -1.,  2.,  2.]), [2, 8, 9])' 
 completed = completedVariants[1]
 toCorrect = completed[0]
 predictors = completed[1]
 corDos = calculateGRM.getCorrectedDosages(toCorrect,predictors,test2arg)
 assert str(corDos) == '(Variant: (20,11000608), [-0.28453360070067379, -1.0125454074933884, -0.11696936309431483, -1.0, -0.054021947254719954, 0.79203273634811333, -2.0, -0.0068409159513448305, 2.2740522457002408, 1.1101062165440048, -1.0, -1.3590869384833639], [3, 10])' 

if ( __name__ == '__main__' ):
 runTest()
