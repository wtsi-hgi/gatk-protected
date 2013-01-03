import calculateGRM
import PlinkBedReader as PlinkReader
import argparse
import numpy
import math
import linear_old as linear
import grmTestDataProvider as provider

class TestException(BaseException):
 def __init__(self,value):
  self.value = value
 def __str__(self):
  return repr(self.value)

def unitTest():
 print("Testing bed file IO...")
 testBedReading()
 print("Bed file IO passed. Testing variant partitioning...")
 accumulateVariantsUnitTest()
 print("Variant partitioning passed. Testing interval processing...")
 intervalsSubsetUnitTest()
 print("Interval processing passed. Testing variant correction...")
 localCorrectionUnitTest()
 print("Variant correction passed. Testing GRM calculation...")
 grmCalculationUnitTest()

def grmCalculationUnitTest():
 args = provider.getArgs(100)
 variants = provider.getIndependentDosages("grmCalculationUnitTest",4,1000,-1)
 dosages = list()
 for v in variants:
  dosages.append(calculateGRM.getCorrectedDosagesNoRegression(v,None))
 dosageMatrix = numpy.matrix([t[1] for t in dosages])
 sdos0 = [t[1][0] for t in dosages]
 sdos1 = [t[1][1] for t in dosages]
 sdos2 = [t[1][2] for t in dosages]
 sdos3 = [t[1][3] for t in dosages]
 dist01 = sum(map(lambda i: sdos0[i]*sdos1[i],range(len(sdos0))))/1000
 dist01prime = calculateGRM.calcDistance(0,1,dosageMatrix,dict(),args)[1]
 if ( dist01 != dist01prime ):
  raise TestException("Error in GRM calculation. Expected: %e, Observed: %e" % (dist01,dist01prime))
 dist22 = sum(map(lambda i: sdos2[i]**2,range(len(sdos2))))/1000
 dist22prime = calculateGRM.calcDistance(2,2,dosageMatrix,dict(),args)[1]
 if ( dist22 != dist22prime ):
  raise TestException("Error in GRM calculation. Expected: %e, Observed: %e" % (dist22,dist22prime))
 dist13 = sum(map(lambda i: sdos1[i]*sdos3[i] if i>99 else 0, range(len(sdos2))))/900
 dist13prime = calculateGRM.calcDistance(1,3,dosageMatrix,dict([[1,range(100)[0:49]],[3,range(100)[49:100]]]),args)
 if ( dist13prime[0] != 900 ):
  raise TestException("Unexpected number of non-missing variants. Expected: 900, observed: %d" % (dist13prime[0]))
 if ( dist13prime[1] != dist13 ):
  raise TestException("Error in GRM calculation. Expected: %e, Observed: %e" % (dist13,dist13prime[1]))

def localCorrectionUnitTest():
 testCorrectionSansRegression()
 testRegressionCorrection()
 testEdgeCaseRegression()
 largeScaleIntegrityTest()

def largeScaleIntegrityTest():
 """ Moving from python to c-bindings creates potential problems with memory issues. Doing lots of regressions
     to test for memory integrity will be useful.
 """
 print("Running large-scale memory integrity test.")
 tolerance = 1e-3
 nTestFiles = 1
 nRunsPerFile = 8
 expectedCoeff = [
   [2.77075768,0.53143971,-0.72988373,-0.26440337,-0.89550549,-1.87032901,0.63513517,-0.78523482,-0.69541818,0.32939047,0.13667907,0.59443712,
    0.63398860,-0.68335889,0.54900151,-0.55853686,0.49036169,-0.35506985,0.56495774,-0.39632489,-0.53151681,-0.48094469,0.03851889,-0.62350338,7.14666839]
   ]
 for fileNo in range(nTestFiles):
  predict = numpy.matrix(list(map(lambda x: list(map(lambda y: float(y),x.strip().split("\t"))),open("test/linLargeScaleP%d.txt" % (fileNo+1)).readlines())))
  response = numpy.array(list(map(lambda x: float(x.strip()),open("test/linLargeScaleR%d.txt" % (fileNo+1)).readlines())))
  for runNo in range(nRunsPerFile):
   result = linear.GLM.Logistic.Fit.newton(response,predict,2)
   coef = result.coefficients
   for i in range(len(coef)):
    if ( abs(coef[i]-expectedCoeff[fileNo][i]) > tolerance ):
     print("Error in run number %d of test %d. Expected: %e   Observed:  %e" % (1+runNo,1+fileNo,coef[i],expectedCoeff[fileNo][i]))

def testEdgeCaseRegression():
 """ Tests some edge cases of regressions that have been found during normal operation of the code.
     These cases tend to yield problems even when the data is exported and the model run in R.
 """
 resp1 = numpy.array(list(map(lambda x: float(x),open("test/resp_break.tsv").readline().strip().split("\t"))))
 pred1 = numpy.matrix(list(map(lambda x: list(map(lambda y: float(y),x.strip().split("\t"))),open("test/pred_break.tsv").readlines())))
 result1 = linear.GLM.Logistic.Fit.newton(resp1,pred1,2)
 if ( abs(result1.residuals[0] - 1.000507 ) < 1e-5 ):
  raise TestException("Bad residual. Expected: %e, observed: %e" % (1.000507,result1.residuals[0]))
 if ( abs(result1.residuals[211] - (-0.427816)) < 1e-5 ):
  raise TestException("Bad residual. Expected: %e, observed: %e" % (-0.427816,result1.residuals[211]))
 if ( abs(result1.residuals[2099] - (-0.7618973)) < 1e-5 ):
  raise TestException("Bad residual. Expected: %e, observed: %e" % (-0.7618973,result1.residuals[2099]))

def testRegressionCorrection():
 """ Tests the accuracy of getCorrectedDosages method. Expected values established by performing appropriate calculation in R. Note that
     the corrected residuals are normalized by the predicted variance, so from R (where you must divide the response by 2 prior to glm):
     (l is the logistic fit, gg is the data frame)
     > predict = 2*predict.glm(l,gg,type="response")[1:10]
     > resid = resp[1:10]-predict
     > predF = predict/2
     > pred_var = sqrt(2*predF*(1-predF))
     > norm_resid = resid/pred_var
 """
 variants = provider.getCorrelatedDosages(0,300,12,100,-1)
 pred = variants[0:3]
 resp = variants[4]
 args = provider.getArgs(100)
 corrected = calculateGRM.getCorrectedDosages(resp,pred,args)
 expectedFirstTenResid = [0.06730088,-0.40206308,-0.40206308,-0.40206308,-0.78505961,2.28614035,-0.40206308,-0.40206308,2.54757725,-0.40206308]
 corResid = corrected[1]
 for idx in range(len(expectedFirstTenResid)):
  if ( abs(corResid[idx] - expectedFirstTenResid[idx]) > 1e-7 ):
   raise TestException("Residual mismatch error. Expected: %e,  Observed: %e" % (expectedFirstTenResid[idx],corResid[idx]))
 variantsWithMissing = provider.getCorrelatedDosages(0,300,5,100,0.01)
 pred = variantsWithMissing[0:3]
 resp = variantsWithMissing[4]
 expectedFirstTenResid = [0.08842582,-0.40163196,-0.40163196,-0.40163196,-0.79978657,2.28902574,-0.40163196,-0.40163196,2.50066716,-0.40163196]
 corrected = calculateGRM.getCorrectedDosages(resp,pred,args)
 corResid = corrected[1]
 for idx in range(len(expectedFirstTenResid)):
  if ( abs(corResid[idx]-expectedFirstTenResid[idx]) > 1e-4 ):
   raise TestException("Residual mismatch error (with missingness). Expected: %e, Observed: %e" % (expectedFirstTenResid[idx],corResid[idx]))
 bigVarWithMissing = provider.getCorrelatedDosages(0,1000,5,100,0.01)
 pred = bigVarWithMissing[0:3]
 resp = bigVarWithMissing[3] 
 #print(pred[0][1])
 #print("-----------")
 #print(pred[0][2])
 #print("===========")
 #print(pred[1][1])
 #print("-----------")
 #print(pred[1][2])
 #print("===========")
 #print(pred[2][1])
 #print("-----------")
 #print(pred[2][2])
 #print("=-=-=-=-=-=")
 #print(resp[1])
 #print("@@@@@@@@@@@")
 #print(resp[2])
 corrected = calculateGRM.getCorrectedDosages(resp,pred,args)
 corResid = corrected[1]
 expected = [-0.06682232,-0.17556269,-0.17556269,-0.46125693,-0.29639402,-0.06682232,-0.09449757,-0.06682232,-0.15953570,-0.17556269,
             -0.46125693,2.17621124,-0.17556269,-0.06682232,-0.29639402,-0.29639402,0.35745650,2.56832439,-0.46125693,-0.22560929] 
 for idx in range(len(expected)):
  if ( abs(corResid[idx]-expected[idx]) > 4e-4 ):
   raise TestException("Residual mismatch error. Expected: %e, Observed: %e" %(expected[idx],corResid[idx]))

def testCorrectionSansRegression():
 """ Tests the accuracy of the getCorrectedDosagesNoRegression method. Expected values established by performing appropriate calculation in R.
 """
 response = provider.getIndependentDosages("testCorrectionSansRegression",50,1,1,-1)[0]
 corrected = calculateGRM.getCorrectedDosagesNoRegression(response,None)
 expected = [-0.8600733,2.3253833,-0.8600733,-0.8600733,-0.8600733,0.7326550,0.7326550,0.7326550,-0.8600733,0.7326550,
              0.7326550,0.7326550,2.3253833,0.7326550,0.7326550,0.7326550,0.7326550,-0.8600733,0.7326550,-0.8600733,
             -0.8600733,-0.8600733,0.7326550,-0.8600733,0.7326550,-0.8600733,-0.8600733,0.7326550,-0.8600733,-0.8600733,
              0.7326550,-0.8600733,0.7326550,-0.8600733,0.7326550,0.7326550,0.7326550,-0.8600733,0.7326550,0.7326550,
              0.7326550,-0.8600733,-0.8600733,0.7326550,-0.8600733,-0.8600733,-0.8600733,-0.8600733,-0.8600733,-0.8600733]
 corDos = corrected[1]
 for idx in range(len(corDos)):
  if ( abs(expected[idx] - corDos[idx]) > 1e-7 ):
   raise TestException("Unexpected dosage correction. Expected: %e  Observed: %e" % (expected[idx],corDos[idx]))
 responseWithMissing = provider.getIndependentDosages("testCorrectionSansRegression",50,1,1,0.1)[0]
 correctedWithMissing = calculateGRM.getCorrectedDosagesNoRegression(responseWithMissing,None)
 corWMDos = correctedWithMissing[1]
 expected = [-2.4563392,-0.8798827,-0.8798827,0.6965738,-0.8798827,0.6965738,2.2730303,0.6965738,0.6965738,0.6965738,-0.8798827,0.6965738,0.6965738,-0.8798827,-0.8798827,0.6965738,0.6965738,0.6965738,0.6965738,0.6965738,0.6965738,-0.8798827,-0.8798827,-0.8798827,-0.8798827,-0.8798827,-2.4563392,-0.8798827,-0.8798827,-0.8798827,0.6965738,0.6965738,0.6965738,0.6965738,2.2730303,0.6965738,-2.4563392,0.6965738,-0.8798827,-2.4563392,0.6965738,-2.4563392,-2.4563392,-0.8798827,-0.8798827,-0.8798827,-0.8798827,-2.4563392,-0.8798827,-0.8798827] 
 for idx in range(len(corWMDos)):
  if ( abs(expected[idx]-corWMDos[idx]) > 1e-7 ):
   raise TestException("Unexpected dosage correction. Expected %e  Observed: %e" % (expected[idx],corWMDos[idx]))

def intervalsSubsetUnitTest():
 """ Tests the accuracy of the removeStaleIntervals, findCompletedVariants, and removeStaleVariants functions
 """
 args = provider.getArgs(500)
 intervals = provider.getStandardIntervals()
 genotypes = provider.getIndependentDosages("intervalsSubsetUnitTest",50,15,250,-1.0)
 regressorVariants = list()
 variantsToCorrect = list()
 dosages = dict()
 individualMissingCalculate = dict()
 individualMissingCorrect = dict()
 accum = lambda x: calculateGRM.accumulateVariant(x,regressorVariants,variantsToCorrect,intervals,args.localCorrectionBP,args,individualMissingCalculate,individualMissingCorrect,len(dosages)) 
 for i in range(14):
  accum(genotypes[i])
 done = calculateGRM.removeStaleIntervals(genotypes[11][0],regressorVariants,intervals,args.localCorrectionBP) 
 if ( done != [("1",11000,12000)] ):
  raise TestException("Error removing stale interval. Expected ('1',11000,12000), observed %s" % (repr(done)))
 completedVariants = calculateGRM.findCompletedVariants(done,regressorVariants,variantsToCorrect,args.localCorrectionBP)
 if ( len(completedVariants) != 5 ):
  raise TestException("Error aggregating variants interval. Expected: 5, observed: %d" %(len(completedVariants))) 
 posList = list(map(lambda t: t[0][0].pos, completedVariants))
 if ( posList != [11000,11250,11500,11750,12000] ):
  raise TestException("Incorrect variants aggregated into interval. Expected [11000,11250,11500,11750,12000], observed %s" % (repr(posList)))
 if ( completedVariants[0][1] != completedVariants[4][1] ):
  raise TestException("Mismatch between variants used for correction at start and end of interval. Start: %s, End: %s" % (repr(completedVariants[0][1]),repr(completedVariants[4][1])))
 corPos = list(map(lambda t: t[0].pos,completedVariants[0][1]))
 if ( corPos != [10500,10750,12250,12500] ):
  raise TestException("Unexpected corrector variant positions. Expected: [10500,10750,12250,12500], observed: %s" % (repr(corPos)))
 calculateGRM.removeStaleVariants(done,intervals,regressorVariants,variantsToCorrect,args.localCorrectionBP)
 regressorPos = list(map(lambda t: t[0].pos,regressorVariants)) 
 correctPos = list(map(lambda t: t[0].pos,variantsToCorrect))
 if ( regressorPos != [12500,12750] ):
  raise TestException("Unexpected regressor variant positions. Expected: [12500,12750], observed: %s" %(repr(regressorPos)))
 if ( correctPos != [13000,13250] ):
  raise TestException("Unexpected variantsToCorrect positions. Expected: [13000,13250], observed: %s" %(repr(correctPos)))

def accumulateVariantsUnitTest():
 singleVariantAccumulationTest()
 missingnessConsistencyTest()

def missingnessConsistencyTest():
 """ The missingness indeces placed into the dictionary by accumulate missing variants should be consistent regardless of the interval used.
     In particular, after removing stale variants, the index counter should keep increasing, and not diverge from the actual variant count.
 """
 args = provider.getArgs(50) # lots more variants here
 genotypes = provider.getIndependentDosages("missingnessConsistencyTest",100,100,args.localCorrectionBP,0.2) # high missingness
 intervals = provider.getStandardIntervals()
 regressorVariants = list()
 variantsToCorrect = list()
 dosages = dict()
 individualMissingCalculate = dict()
 individualMissingCorrect = dict()
 accum = lambda x: calculateGRM.accumulateVariant(x,regressorVariants,variantsToCorrect,intervals,args.localCorrectionBP,args,individualMissingCalculate,individualMissingCorrect,len(dosages))
 removeStaleInterval = lambda x: calculateGRM.removeStaleIntervals(x[0],regressorVariants,intervals,args.localCorrectionBP)
 getCompleted = lambda x: calculateGRM.findCompletedVariants(x,regressorVariants,variantsToCorrect,args.localCorrectionBP)
 removeStaleVariant = lambda x: calculateGRM.removeStaleVariants(x,intervals,regressorVariants,variantsToCorrect,args.localCorrectionBP)
 for geno in genotypes[0:41]:
  accum(geno)
  doneInt = removeStaleInterval(geno)
  cmp = getCompleted(doneInt)
  assert cmp == []
  removeStaleVariant(doneInt)
 # variantsToCorrect should consist of all variants from 11000 to 12000 by 50s
 if ( len(variantsToCorrect) != 21 ):
  raise TestException("Unexpected number of variants to correct. Expected: 21, observed: %d" %(len(variantsToCorrect)))
 if ( individualMissingCalculate[0] != [1,12,16] ):
  raise TestException("Individual number 1 has mismatching variant missing indeces. Expected: [1,12,16], observed: %s" % (repr(individualMissingCalculate[0])))
 if ( individualMissingCalculate[9] != [1,5,17] ):
  raise TestException("Individual number 9 has mismatching variant missing indeces. Expected: [1,5,17], observed: %s" % (repr(individualMissingCalculate[9])))
 for geno in genotypes[41:]:
  accum(geno)
  doneInt = removeStaleInterval(geno)
  cmp = getCompleted(doneInt)
  for c in cmp:
   dosages[c[0][0]]=c[0][1] 
  removeStaleVariant(doneInt)
 if ( individualMissingCalculate[4] != [0,2,4,11,17,21,26,27,28,32,33,34,41]): 
  raise TestException("Individual number 4 has mismatching variant missing indeces. Expected: [0,2,4,11,17,21,26,27,28,32,33,34,41], observed: %s" % (repr(individualMissingCalculate[4])))
 if ( individualMissingCalculate[20] != [0,17,21,25,32,33] ):
  raise TestException("Individual number 20 has mismatching variant missing indeces. Expected: [0,17,21,25,32,33], observed: %s" % (repr(individualMissingCalculate[20])))
 if ( individualMissingCalculate[37] != [7,10,13,20,35,38,40,41] ):
  raise TestException("Individual number 37 has mismatching variant missing indeces. Expected: [7,10,13,20,35,38,40,41], observed: %s" % (repr(individualMissingCalculate[37])))

def singleVariantAccumulationTest():
 """ Tests the accuracy of the accumulateVariant method
 """
 # test where a variant gets put wrt the intervals
 args = provider.getArgs(500)
 intervals = provider.getStandardIntervals()
 genotypes = provider.getIndependentDosages("singleVariantAccumulationTest",50,20,250,-1.0)
 regressorVariants = list()
 variantsToCorrect = list() 
 dosages = dict()
 individualMissingCalculate = dict()
 individualMissingCorrect = dict()
 accum = lambda x: calculateGRM.accumulateVariant(x,regressorVariants,variantsToCorrect,intervals,args.localCorrectionBP,args,individualMissingCalculate,individualMissingCorrect,len(dosages))
 # first variant should be ignored. It's only at 10000
 geno = genotypes[0]
 assert geno[0].pos == 10000
 accum(geno)
 if ( regressorVariants != [] or variantsToCorrect != [] ):
  raise TestException("Unexpected change in lists. Expected: len(regressor): 0  len(toCorrect): 0, observed: %d  %d" %(len(regressorVariants),len(variantsToCorrect)))
 # this should also be ignored
 accum(genotypes[1])
 # this goes into the regressor variants
 accum(genotypes[2])
 accum(genotypes[3])
 distance = calculateGRM.calcDistanceToInterval(genotypes[2][0],intervals[0])
 if ( distance != 500 ):
  raise TestException("Distance to beginning of interval miscalculated. Expected: 500, observed: %d" %(distance))
 if ( len(regressorVariants) != 2 ):
  raise TestException("Failure to aggregate nearby variants prior to interval. Exptected: 2, observed: %d" % (len(regressorVariants)))
 if ( len(variantsToCorrect) != 0 ):
  raise TestException("Nearby out-of-interval variants classified as within interval. Expected: 0, observed: %d" %(len(variantsToCorrect)))
 # next five 11(000,250,500,750,000) should go into the interval itself (variants to correct)
 accum(genotypes[4])
 accum(genotypes[5])
 accum(genotypes[6])
 accum(genotypes[7])
 accum(genotypes[8])
 if ( len(regressorVariants) != 2 ):
  raise TestException("Inside-interval variants classified as near interval. Expected: 2, observed: %d" % (len(regressorVariants)))
 if ( len(variantsToCorrect) != 5 ):
  raise TestException("Unexpected number of aggregated within-interval variants. Expected: 4, observed: %d" %(len(variantsToCorrect)))
 assert genotypes[8][0].pos == 12000
 # next two should accumulate into the nearby variants
 accum(genotypes[9])
 accum(genotypes[10])
 if ( len(regressorVariants) != 4 ):
  raise TestException("Failure to aggregate nearby variants posterior to interval. Exptected: 4, observed: %d" % (len(regressorVariants)))
 if ( len(variantsToCorrect) != 5 ):
  raise TestException("Nearby out-of-interval variants classified as within interval. Expected: 5, observed: %d" %(len(variantsToCorrect)))
 # this one is nearby - now to the next interval rather than the top one
 accum(genotypes[11])
 if ( len(regressorVariants) != 5 ):
  raise TestException("Failure to aggregate nearby variants posterior to interval. Exptected: 5, observed: %d" % (len(regressorVariants)))
 if ( len(variantsToCorrect) != 5 ):
  raise TestException("Nearby out-of-interval variants classified as within interval. Expected: 5, observed: %d" %(len(variantsToCorrect)))
 # this one enters the next interval
 accum(genotypes[12])
 if ( len(regressorVariants) != 5 ):
  raise TestException("Inside-interval variants classified as near interval. Expected: 5, observed: %d" % (len(regressorVariants)))
 if ( len(variantsToCorrect) != 6 ):
  raise TestException("Unexpected number of aggregated within-interval variants. Expected: 6, observed: %d" %(len(variantsToCorrect)))
 # finally this one should have no effect
 accum(genotypes[19])
 if ( len(regressorVariants) != 5 ):
  raise TestException("Downstream variant accumulated into nearby variants")
 if ( len(variantsToCorrect) != 6 ):
  raise TestException("Downstream variant accumulated into interval variants")
 # since that's ignored we can just test this. The first interval should be popped.
 doneInt = calculateGRM.removeStaleIntervals(genotypes[10][0],regressorVariants,intervals,args.localCorrectionBP)
 if ( doneInt != [] ):
  raise TestException("Interval being improperly removed. Interval: %s, distance: %s, variant: %s" %(repr(doneInt),repr(args.localCorrectionBP),repr(genotypes[10][0])))
 doneInt = calculateGRM.removeStaleIntervals(genotypes[11][0],regressorVariants,intervals,args.localCorrectionBP)
 if ( doneInt != [('1', 11000, 12000)] ):
  raise TestException("Error popping out completed interval. Expected [('1', 11000, 12000)], observed %s" % (repr(doneInt)))

def testBedReading():
 singleVariantBedTest()
 multiVariantBedTest()
 optimizedBedTest()

def optimizedBedTest():
 args = provider.getArgs(0,"./test/1000G_subset.chr20.79234")
 reader = PlinkReader.SiteOptimizedPlinkBinaryReader(args.bedBase)
 if ( not reader.numGenotypesPerMajor == 12 ):
   raise TestException("Mismatching samples. Expected: %d, observed: %d" %(12,reader.numGenotypesPerMajor))
 genotypeDosages = calculateGRM.getNextVariantOptimized(reader,args)
 reader2 = PlinkReader.getReader(args.bedBase)
 dosages2 = calculateGRM.getNextVariant(reader2,args)
 if ( str(genotypeDosages) != str(dosages2) ):
   raise TestException("Mismatching dosage items. Expected: %s, observed: %s" % (str(genotypeDosages),str(dosages2)))

def singleVariantBedTest():
 args = provider.getArgs(0,"./test/1000G_subset.chr20.79234")
 reader = PlinkReader.getReader(args.bedBase)
 if ( not reader.snpMajor() ):
  raise TestException("The test bed file ./test/1000G_subset.chr20.79234 is being read in individual mode rather than snp major mode.")
 if ( not reader.numGenotypesPerMajor == 12 ):
  raise TestException("Mismatching samples. Expected: %d, observed: %d" % (12,reader.numGenotypesPerMajor))
 genotypeDosages = calculateGRM.getNextVariant(reader,args)
 if ( not genotypeDosages[0].pos == 79234 ):
  raise TestException("Mismatching variant positions. Expected: %d, observed: %d" % (79234,genotypeDosages[0].pos))
 if ( not genotypeDosages[0].frequency == (11.0)/24 ):
  raise TestException("Mismatching variant frequencies. Expected: %f, observed: %f" % (11.0/24,genotypeDosages[0].pos))
 dosages = genotypeDosages[1]
 if ( not numpy.sum(dosages) == 11 ):
  raise TestException("Mismatching dosage sums. Expected: %f, observed: %f" % (11,numpy.sum(dosages)))
 if ( not numpy.linalg.norm(dosages)**2 == 17 ):
  raise TestException("Mismatching dosage norms. Expected: %f, observed: %f" % (17,numpy.linalg.norm(dosages)**2))
 nextDosages = calculateGRM.getNextVariant(reader,args)
 if ( not nextDosages == None ):
  raise TestException("Reading beyond end of single-variant bed file. Expected None, observed: "+str(nextDosages)) 
 assert str(genotypeDosages) == "(Variant: (20,79234), array([ 0.,  0.,  2.,  0.,  2.,  2.,  1.,  0.,  1.,  1.,  1.,  1.]), [])"  

def multiVariantBedTest():
 args = provider.getArgs(0,"./test/1000G_subset.chr20")
 reader = PlinkReader.getReader(args.bedBase)
 allDosages = list()
 genotypeDosages = calculateGRM.getNextVariant(reader,args)
 while ( genotypeDosages != None ):
  allDosages.append(genotypeDosages)
  genotypeDosages = calculateGRM.getNextVariant(reader,args)
 if ( len(allDosages) != 23919 ):
  raise TestException("Mismatching number of variants. Expected: %d, observed: %d" %(23919,len(allDosages))) 
 if ( allDosages[len(allDosages)-1][0].pos != 19999982 ):
  raise TestException("Mismatching position. Expected: %d, observed: %d" % (19999982,allDosages[len(allDosages)-1][0].pos))
 freqSum = sum(map(lambda t: t[0].frequency, allDosages))
 if ( abs(freqSum - 7854.007576) > 1e-6 ):
  raise TestException("Mismatching frequencies. Expected: %f, observed: %f, diff %e" % (7854.007576,freqSum,7854.007576-freqSum))

if ( __name__ == '__main__' ):
 unitTest() 
