import cProfile
import PlinkBedReader as PlinkReader
import argparse
import numpy
import math
import linear_old as linear
import functools

reduce = functools.reduce

def parseArgs():
 parser = argparse.ArgumentParser(description='Parse out arguments for the generalized burden test')
 parser.add_argument("-grm",action='store',default=None,dest="grm",
               required=True,help="The file into which to write the GRM values.")
 parser.add_argument("-bedBase",action='store',default=None,dest="bedBase",
               required=True,help="The base of the bed (so 'xxxx' if file named 'xxxx.bed'). Requires a bim and fam with the same base name.")
 parser.add_argument("-intervals",action='store',default=None,dest="intervals",required=False,help="A file listing intervals over which the GRM should be calculated")
 parser.add_argument("-localCorrectionBP",action='store',default=-1,type=int,dest="localCorrectionBP",required=False,
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
 parser.add_argument("-opt",action="store_true",dest="opt",required=False,
         help="Use optimized bed reader. Experimental.")
 return parser.parse_args()

## stream in the SNPs, break them into bins of "independent correctors" and "dependent genotypes", and run a correction
def streamAndLocalCorrect(args,reader):
 """ Stream in genotypes from the reader and obtain dosages that have been locally corrected by regions outside of the GRM intervals.
     Variants are broken into two bins: "independent correctors" and "dependent genotypes".
     @args: An argument parser object. We will use the following arguments it contains
           + args.intervals - an interval list which contains those intervals the GRM should be calculated from 
           + args.localCorrectionBP - the number of BP such that, if a variant falls within X of an interval, it is used for local correction
     @reader: a PlinkBedReader that streams in genotypes from the bed file passed in
     @return: a list of corrected dosages. This is a list of tuples: List[(PlinkVariant,List[Float])]
 """
 ## todo -- enable the expansion factor to change based on allele frequency
 dependentGenotypeIntervals = getIntervals(args.intervals) 
 regressorVariants = list()
 variantsToCorrect = list() 
 dosages = dict()
 individualMissingCalculate = dict()
 individualMissingCorrect = dict()
 if ( args.opt ):
  nextVar = getNextVariantOptimized
 else:
  nextVar = getNextVariant
 genotypes = nextVar(reader,args)
 # this is a 3-ple: (PlinkVariant,List[Float],List[Int])
 siteNo = 0
 while ( genotypes != None and len(dependentGenotypeIntervals) > 0):
  if ( siteNo % 50 == 1 ):
   print("read : "+str(siteNo) + " at : " + str(genotypes[0].pos) + " accum : "+str(len(regressorVariants))+ " compute : "+str(len(variantsToCorrect)))
  accumulateVariant(genotypes,regressorVariants,variantsToCorrect,dependentGenotypeIntervals,args.localCorrectionBP,args,individualMissingCalculate,individualMissingCorrect,len(dosages)) 
  doneIntervals = removeStaleIntervals(genotypes[0],regressorVariants,dependentGenotypeIntervals,args.localCorrectionBP)
  completedVariants = findCompletedVariants(doneIntervals,regressorVariants,variantsToCorrect,args.localCorrectionBP)
  removeStaleVariants(doneIntervals,dependentGenotypeIntervals,regressorVariants,variantsToCorrect,args.localCorrectionBP)
  for regressionVariants in completedVariants:
   varToCorrect = regressionVariants[0]
   varsToCorrectWith = regressionVariants[1]
   correctedVar,correctedDosages = getCorrectedDosages(varToCorrect,varsToCorrectWith,args) 
   dosages[correctedVar] = correctedDosages
  genotypes = nextVar(reader,args)
  siteNo += 1
 # there may yet be intervals not totally complete
 for completedVariant in findCompletedVariants(dependentGenotypeIntervals,regressorVariants,variantsToCorrect,args.localCorrectionBP):
  varToCorrect = completedVariant[0]
  varsToUse = completedVariant[1]
  correctedVar,correctedDosages = getCorrectedDosages(varToCorrect,varsToUse,args)
  dosages[correctedVar]=correctedDosages
  # why is this a map?  
 return (dosages,individualMissingCalculate)

def getCorrectedDosagesNoRegression(responseVar,predictorVar):
 """ A regression-free version of getCorrectedDosages that does not do regression, but only corrects by the mean and variance of the response.
     No-calls are not dealt with here. In the non-debug version, they are replaced with mean values (2*variant.frequency) prior to regression.
     Missing-value indeces are not propagated forward, as they are maintained (in a more useful orientation) by individual missing dictionaries.
     @responseVar - a 3ple (Variant,numpy.array[Double],List[Int]) of the GRM-relevant variant site, its dosages, and the list of individuals that are missing values
     @predictorVar - a List[(Variant,numpy.array[Double],List[Int])] of the LD-correction predictor variant sites, their dosages, and thier missing individuals lists
     @return - a 2ple (Variant,numpy.array[Double]) where, if the frequency of the site is f, the dosages have been corrected by (x-f)/sqrt(f(1-f))
 """
 # just pass the response through normally, subtracting the mean
 mean = responseVar[0].frequency
 if ( mean < 1e-9 or mean > 1-1e-9 ):
  std = 1
 else:
  std = math.sqrt(2*mean*(1.0-mean))
 mean = 2*responseVar[0].frequency
 respArray = (responseVar[1] - mean)/std
 return (responseVar[0], respArray)

def formatResponse(respVar):
 """ A helper function that takes a list of doubles and reformats it as a numpy.array of doubles. This list is housed inside of a 2-ple.
     @respVar - a 2ple (PlinkVariant,List[Float]) of a variant site and the per-sample dosages.
     @return - a 2ple (PlinkVariant,numpy.array[Float]) of a variant site and the per-sample dosages.
 """
 return (respVar[0],numpy.array(respVar[1]))

def getCorrectedDosages(responseVar,predictorVar,args):
 """ The workhorse method for LD correction. Given dosages (in @responseVar), calculate a set of corrected dosages based on those
     in @predictorVar. If @predictorVar is empty, this is simply normalizing by mean and variance. If it has dosages in it, then
     this is a bit more complicted: call into the linear library to run logistic regression of @responseVar on @predictorVar,
     and take the residuals as dosages. These residuals will have a very different variance from @responseVar, which has var = p(1-p).
     instead if you work out the variance of a bernoulli mix, it's sum p_i(1-p_i), and thus here to correct for "frequency" it should be
     through the "average" variance: 1/n sum p_i|G(1-p_i|G). However, for pairwise comparisons, the variance is sample-dependent (depending, of course
     on the surrounding haplotype), and so it should be (r1*r2)/(sqrt(v1)*sqrt(v2)) for residuals and variance. Thus normalize each residual by the
     variance of the predicted value. For r^2 = 0 or r^2 = 1, this reduces to sqrt(p(1-p)).
     @responseVar - a 3ple (PlinkVariant, List[Float],List[Int]) of the variant site and per-sample dosages to be corrected
     @predictorVar - a list of 3ples List[(PlinkVariant,List[Float],List[Int])] of variant sites and per-sample dosages to use for correction
     @args - the command line arguments, to check if missing values have been fixed or whether they should be excluded
     @return - A (PlinkVariant,numpy.array([Float])) pair that gives the variant site, and the corrected frequencies.
 """ 
 if ( len(predictorVar) == 0 ):
  return getCorrectedDosagesNoRegression(responseVar,predictorVar)
 # 1: stick the response into a numpy vector array
 if ( args.naCorrect ):
  respArray = numpy.array(responseVar[1])
  predArray = numpy.matrix(list(map(lambda u: u[1], predictorVar)))
 else:
  respMissing = responseVar[2]
  respArray = list(map(lambda x: responseVar[1][x], filter(lambda u: u not in respMissing, range(len(responseVar[1])))))
  respArray = numpy.array(respArray) 
  # 2: stick the predictors into a numpy matrix array
  predArray = list()
  for pv in predictorVar:
   # missing values in predictors always set to twice frequency (i.e. mean dosage)
   for mv in pv[2]:
    pv[1][mv] = pv[0].frequency*2
   predArray.append(list(map(lambda i: pv[1][i], filter(lambda u: u not in respMissing, range(len(pv[1]))))))
  predArray = numpy.matrix(predArray)
 # add ones to the predictor array
 predArray = numpy.append(predArray,numpy.matrix(numpy.ones(predArray.shape[1])),0)
 # 3: run a regression on this
 if ( respArray.shape[0] - len(predictorVar) < 1 ):
  raise BaseError("Too many variants for the sample size. Try decreasing localCorrectBP.")
 elif ( respArray.shape[0] - len(predictorVar) > 50 ):
  regression = linear.GLM.Logistic.Fit.newton(respArray,predArray.T,2)
 else:
  regression = linear.GLM.Logistic.Fit.gradientDescent(respArray,predArray.T,2)
 # 4: if missing values have been removed, put the dosages back in
 predP = regression.predictedValues/2
 normalizedResids = regression.residuals/numpy.sqrt(2*predP*(1-predP))
 if ( not args.naCorrect ):
  fullIdx = 0
  residIdx = 0
  finalDosages = list()
  while ( fullIdx < len(responseVar[1]) ):
   if ( fullIdx in responseVar[2] ):
    finalDosages.append(-1.0)
   else:
    finalDosages.append(normalizedResids[residIdx])
    residIdx += 1
   fullIdx += 1
  return (responseVar[0],finalDosages)
 return (responseVar[0],normalizedResids)


def globalCorrect(dosage,globalGenotypes):
 # another workhorse method: dosage contains residuals from earlier regression. Fit these with the global genotypes.
 predArray = numpy.array(list(map(lambda u: numpy.array(list(map(lambda v: v.getDosage(),u[1]))),globalGenotypes)))
 regression = linear.GLM.Logistic.simpleNewton(dosage,predArray,N)
 return regression.residuals

def globalCorrectDebug(dosage,globalGenotypes):
 # just a verison of the above that ignores regression for debugging purposes. Pass through the dosages.
 return [dosage[t] for t in sorted(dosage)]

def findCompletedVariants(completedIntervals,regressorVariants,variantsToCorrect,distance):
 """ for intervals that are complete, combine the nearby (within @distance bp) regressor variants and those variants to correct within the interval.
     This generates a list of (variantToCorrect,List[RegressorVariant]) pairs, one for each completed interval.
     @completedIntervals - a list of intervals that were determined complete. Will extract regressor/regressand variables based on this list. A list: List[(String,Int,Int)]
     @regressorVariants - a list of variants falling close to intervals to be used for local correction. A list: List[(PlinkVariant,List[PlinkGenotype])]
                           This algorithm will pair variants in this list with variants in @variantsToCorrect based on proximity to the same interval
     @variantsToCorrect - a list of variants falling inside intervals from which the corrected dosage is to be extracted. A list: List[(PlinkVariant,List[PlinkGenotype])]
                           This method will pair variants in this list with variants in @regressorVariants based on proximity to the same interval
     @returns - a list of pairs: (PlinkVariant,List[PlinkVariant]) extacted from the finished intervals representing a dosage-correction regression. E.g. for (Y,X): Y ~ X.
 """
 marshalledRegressions = list()
 for completedInterval in completedIntervals:
  # assemble those regressors that are nearby. Perhaps there's a faster way to do this (like a cached map from intervals to the first and last regressor variants close to it)
  intervalRegressors = list()
  for regressorGenotype in regressorVariants:
   regressorVariant = regressorGenotype[0]
   if ( calcDistanceToInterval(regressorVariant,completedInterval) < -distance ):
    break
   elif ( calcDistanceToInterval(regressorVariant,completedInterval) <= distance ):
    intervalRegressors.append(regressorGenotype)
  # now pair the variants to correct in the interval with these regressors
  for correctGenotype in variantsToCorrect:
   if ( calcDistanceToInterval(correctGenotype[0],completedInterval) < 0 ):
    break
   elif ( calcDistanceToInterval(correctGenotype[0],completedInterval) == 0 ):
    marshalledRegressions.append((correctGenotype,intervalRegressors))
 return marshalledRegressions

def getRegressionVariants(variant,accumulator,distance):
 ## there is a faster way to do this by keeping track of the index of the variant in findCompletedVariants. If too slow can fix here.
 nearVariants = list()
 for v in accumulator:
  if ( abs(variant[0].distanceTo(v[0])) <= distance ):
   nearVariants.append(v)
 return (variant,nearVariants) 
 
def accumulateVariant(genotypes,correctorVariantsList,variantsToCorrectList,intervalsToCorrectList,distance,args,individualMissingCalculate,individualMissingCorrect,calcVarComplete):
 """ Takes a variant genotype and identifies the correct bin in which to put it. This identification is solely dependent on
     the distance between the variant and its closest interval. In particular:
      + > @distance BP == ignore
      + <= @distance BP, but > 0bp == it is a corrector variant
      + =0 BP == it is a variant to correct
     @genotypes - The variant site and genotypes to place into bin. A 3-ple (PlinkVariant,List[Double],List[Int]).
     @correctorVariantsList - the list of variants to use for downstream genotype correction. A list: List[(PlinkVariant,List[Double])]
     @variantsToCorrectList - the list of variants to correct and get dosages for during downstream correction. A list: List[(PlinkVariant,List[Double])]
     @intervalsToCorrectList - the list of intervals that determine which variants to correct. A list: List[(String,Int,Int)]
     @distance - the number of BP from the start or end of an interval for which variants should be used for local correction. An Int.
     @args - the command line arguments (for minimum and maximum frequency cutoffs)
     @individualMissingCalculate - a dict { Int -> List[Int] } mapping individual indeces to indeces of variants in variantsToCorrectList that are missing genotypes
     @individualMissingCorrect - a dict { Int -> List[Int] } mapping individual indeces to indeces of variants in correctorVariantsList that are missing genotypes
     @calcVarComplete - an Int - the number of variants used for the final calculation that have been corrected to this point. Used to get the proper index for missingness.
     @return - none
     @sideEffects - Update one of (@correctorVariantsList,@individualMissingCorrect) or (@variantsToCorrectList,@individualMissingCalculate)
 """
 freq = genotypes[0].frequency 
 inInterval = False
 nearInterval = False
 variant = genotypes[0] # get the PlinkVariant
 for correctInterval in intervalsToCorrectList:
  distToInterval = calcDistanceToInterval(variant,correctInterval)
  if ( distToInterval > distance ):
   # more than @distance to the interval. Break.
   break
  elif ( distToInterval > 0 ):
   # within distance of the interval, but still before
   nearInterval = True
  elif ( distToInterval == 0 ):
   # inside of the interval
   inInterval = True
  elif ( distToInterval >= -distance ):
   # after the interval, but still within distance
   nearInterval = True
  else:
   # way after the interval, continue
   continue
 if ( inInterval ):
  if( freq < args.maxFrequency and freq > args.minFrequency ):
   variantsToCorrectList.append(genotypes)
   updateMissing(individualMissingCalculate,genotypes[2],calcVarComplete + len(variantsToCorrectList)-1)
 elif ( nearInterval and freq < args.maxFrequencyCorrect and freq > args.minFrequencyCorrect ):
  correctorVariantsList.append(genotypes)
  updateMissing(individualMissingCorrect,genotypes[2],calcVarComplete + len(correctorVariantsList)-1)

def updateMissing(individualMissingMap,individualIndeces,variantIdx):
 """ For each individual in @individualIndeces, appends @variantIdx to the list of missing genotypes (by site) in the missing map
     @individualMissingMap - a dict { Int -> List[Int] } mapping individual indeces to indeces of variants that are missing genotypes
     @individualIndeces - a List[Int] of individual indeces that are missing this locus
     @variantIdx - the index of this locus
     @return - nothing
     @sideEffect - individualMissingMap is updated
 """
 for iIdx in individualIndeces:
  if ( iIdx not in individualMissingMap ):
   individualMissingMap[iIdx] = list()
  individualMissingMap[iIdx].append(variantIdx)

def removeStaleIntervals(variant,correctorVariantsList,intervalsToCorrectList,distance):
 """ Examine the lists to identify, remove, and return intervals that are stale (e.g. intervals far enough from varaints as to never include a variant again)
     @variant - the most recent variant read from the bedFile. Used to cull intervals out.
     @correctorVariantsList - the list of variants to use for downstream genotype correction. A list: List[(PlinkVariant,List[PlinkGenotype])]. Remove variants if too far before first interval.
     @intervalsToCorrectList - the list of intervals that determine which variants to correct. A list: List[(String,Int,Int)]. Remove intervals if too far before current variant.
     @distance - the number of BP from the start or end of an interval for which variants should be used for local correction. An Int. 
     @return - a list of intervals removed
 """
 intervalsRemoved = list()
 while ( len(intervalsToCorrectList) > 0 and calcDistanceToInterval(variant,intervalsToCorrectList[0]) < -distance ):
  intervalsRemoved.append(intervalsToCorrectList.pop(0))
 return intervalsRemoved

def removeStaleVariants(completedIntervalsList,remainingIntervalsList,regressorVariantsList,variantsToCorrectList,distance):
 """ The completed intervals describe a set of variants that are done: namely all of the variants to correct within the
     interval, and those coming before it. Need to be careful about not removing anything involved in remaining intervals, and re-mapping missing indeces
     @completedIntervalsList - a list of intervals for which we have seen a variant more than @distance after its end position. A list: List[(String,Int,Int)]
     @remainingIntervalsList - a list of intervals for which we have *not yet* seen a variant more than @istance after their end position. A list: List[(String,Int,Int)]
     @regressorVariantsList - a list of variants to use as independent variables for predicting variants in @variantsToCorrect. A list: List[(PlinkVariant,List[PlinkGenotype])]
                              This method will remove from this list all variants falling near the completed intervals and not within @distance of the remaining intervals
     @variantsToCorrectList - a list of variants from which to obtain corrected dosages via regression on @regressorVariantsList. A list: List[(PlinkVariant,List[PlinkGenotype])]
                              This method will remove from this list all variants falling inside of completed intervals.
     @distance - The number of BP from the start or end of an interval for which variants should be used for local correction. An Int.
     @return - nothing. But modify the lists above in the ways mentioned.
 """
 if ( len(remainingIntervalsList) == 0 ):
  # there are no more intervals. Clear everything.
  regressorVariantsList = list()
  variantsToCorrectList = list()
 else:
  firstRemainingInterval = remainingIntervalsList[0]
  for completedInterval in completedIntervalsList:
   # remove the variants near the interval
   while ( len(regressorVariantsList) > 0 and calcDistanceToInterval(regressorVariantsList[0][0],completedInterval) > -distance and
           calcDistanceToInterval(regressorVariantsList[0][0],firstRemainingInterval) > distance ):
    regressorVariantsList.pop(0)
   # remove the variants in the interval
   while ( len(variantsToCorrectList) > 0 and calcDistanceToInterval(variantsToCorrectList[0][0],completedInterval) == 0 ):
    variantsToCorrectList.pop(0)

def calcDistanceToInterval(var,interval):
 if ( var.chrStr != interval[0] ):
  return 500000000
 varpos = var.pos
 if ( varpos > interval[2] ): # after the interval
  return -(varpos - interval[2])
 elif ( varpos < interval[1] ): # before the interval
  return interval[1]-varpos
 else:
  # must be in the interval
  return 0

def getIntervals(intervalFile):
 return list(map(lambda y: (y[0],int(y[1]),int(y[2])), map(lambda x: x.strip().split(), open(intervalFile).readlines())))

def getNextVariant(reader,args):
 """ Given a plink binary bed file reader (in SNP major mode), aggregates genotypes across individuals until the whole variant has been accumulated.
     Genotypes are not necessary at this point: only dosages. Thus the genotype objects are discarded in favor of doubles.
     @reader - a plink binary reader
     @args - the command line arguments. For deciding whether to correct missing dosages, or keep track of them.
     @return - a 3ple: (Variant, List[Double], List[Int]) of the Variant object representing the site information, a list of genotype dosages (converted to a numpy array),
                  and a list of individuals missing genotypes (by index)
 """
 try:
  # accumulate all samples into a single tuple of (common variant, list<genotype>)
  genotypes = list() 
  sIdx = 0
  sumDos = 0.0
  an = 0
  variant = None
  missing = list()
  while ( sIdx < reader.numGenotypesPerMajor ):
   geno = reader.next()
   if ( variant == None ):
    variant = geno.variant
   dos = geno.getDosage()
   if ( dos > -1 ):
    sumDos += dos
    an += 2
   else:
    missing.append(sIdx)
   genotypes.append(float(dos))
   sIdx += 1
  freq = sumDos/an
  variant.frequency = freq
  if ( args.naCorrect ):
   for t in missing:
    genotypes[t] = 2*freq
  return (variant,numpy.array(genotypes),missing)
 except StopIteration:
  return None

def psum(x,y):
 """ A version of summing x = x + y where values of y < 0 are ignored
 """
 return x if y < 0 else y+x

def getNextVariantOptimized(reader,args):
 """ An optimized version of getNextVariant. Assumes the reader is the optimized reader, and calls the function to directly get the dosages.
     Missingness and frequency still need to be calculated for the site.
     @reader - a plink binary reader
     @args - the command line arguments. For deciding whether to correct missing dosages, or keep track of them.
     @return - a 3ple: (Variant, List[Double],List[Int]) of the Variant object representing the site information, a list of genotype dosages (as a numpy array),
               and a list of individuals missing gneotypes (by index)
 """
 try:
   var,dosage = reader.nextDosageFloat()
   missing = [i for i in range(len(dosage)) if dosage[i] < 0] 
   freq = reduce(psum,dosage)/(2*(len(dosage)-len(missing)))
   var.frequency = freq
   if ( args.naCorrect ):
     for t in missing:
       dosage[t] = freq
   return (var,numpy.array(dosage),missing)
 except StopIteration:
   return None

## given the dosages, we want to further correct by some global set of SNPs (e.g. known causal/associated loci)
def globallyCorrectDosages(args,dosages,individualMissingDosages):
 # need to read in the global variants. This means re-reading, but whatever.
 bedReader = PlinkReader.getReader(args.bedBase)
 # get the IDs of the global correctors
 globalIDs = set(map(lambda x: x.strip(),open(args.globalCorrection).readlines()))
 # do some matching
 globalGenotypes = list()
 var = getNextVariant(reader)
 while ( var != None ):
  if ( var[0].id in globalIDs ):
   globalGenotypes.append(var)
 # do some corrections
 correctedDosages = list()
 vars = set(dosages.keys())
 vars.sort()
 for v in vars:
  dosage = dosages[v]
  corrected = globalCorrectDebug(dosage,globalGenotypes) 
  correctedDosages.append(corrected)
 return correctedDosages

def printGRMFromDosages(args,correctedDosageMap,individualMissingDosages):
 """ Given the dosages corrected for local haplotype structures, calculate the pairwise sample distances (including inbreeding, e.g. self-self),
     and output the resulting matrix in linearized form.
     @args - the command line argument structure.
     @correctedDosageMap - a list of corrected dosages. This is a list of tuples: List[(PlinkVariant,List[Float])]
     @individualMissingDosages -  a dict {Int -> List[Int]} mapping individual index to the list of variant indeces that are missing genotypes for that individual
 """
 out = open(args.grm,'w')
 # get the dosages
 correctedDosages = numpy.matrix([correctedDosageMap[t] for t in correctedDosageMap])
 # the dosages come in as a list of lists, where correctedDosages[i] are the dosages for snp [i]. So for dosages[i][j] we need to iterate over [j] first.
 for samIdx1 in range(correctedDosages[0,:].size):
  for samIdx2 in range(samIdx1+1):
   nVariants,relatedness = calcDistance(samIdx1,samIdx2,correctedDosages,individualMissingDosages,args)
   out.write("%d\t%d\t%d\t%.8e\n" % (1+samIdx1,1+samIdx2,nVariants,relatedness))

def calcDistance(idx1,idx2,dosages,missing,args):
 """ compute the relatedness between sample 1 and sample 2, given the dosages
     dosages come in as a matrix dosage[snp][sample] = value
      note on this calculation:
      GCTA uses the formula Ajk = 1/N sum_i (Xij - 2pi)(Xik-2pi)/2pi(1-pi)
      the equivalent here would be to do the following:
       > numVar = len(dosages[0])
       > meanDos = map(lambda x: mean(dosages[x]),range(numVar))
       > sam1dos = map(lambda x: dosages[x][idx1]-2*meanDos[x],range(numVar))
       > sam2dos = map(lambda x: dosages[x][idx2]-2*meanDos[x],range(numVar))
       > components = map(lambda x: sam1dos[x]*sam2dos[x]/(2*meanDos[x]*(1-meanDos[x])),range(numVar))
       > relatedness = sum(components)/numVar
      however, it is not entirely clear that the assumptions built into this model
      namely either that effect size ~ sqrt(p(1-p)) or that variance = p(1-p) --
      should hold with LD-corrected dosages.
      for now - don't bother normalizing. The mean dosage should already be 0, from the regression.
   @idx1 - the index for sample 1
   @idx2 - the index for sample 2
   @dosages - a numpy.matrix[Float] of normalized dosages: rows SNP, columns samples
   @missing - a dict {Int -> List[Int]} mapping sample indexes to the list of SNP indexes that are missing for those samples
   @args - the command line arguments
 """
 # todo -- make the behavior togglable via a command line argument.
 # todo -- this is just an inner product. numpy should be able to speed it up.
 sam1Dos = dosages[:,idx1]
 sam2Dos = dosages[:,idx2]
 # get the missing indeces
 sam12missing = set()
 if ( idx1 in missing ):
  sam12missing.update(missing[idx1])
 if ( idx2 in missing ):
  sam12missing.update(missing[idx2])
 if ( not args.naCorrect ):
  sam12missing = list(sam12missing)
  # extract the values
  sam1val = sam1Dos[sam12missing]
  sam2val = sam2Dos[sam12missing]
  # set values to 0
  sam1Dos[sam12missing] = 0
  sam2Dos[sam12missing] = 0
 # inner product
 dist = numpy.dot(sam1Dos.T,sam2Dos)[0]
 # normalize
 nVar = sam1Dos.size - len(sam12missing)
 dist /= nVar
 if ( not args.naCorrect ):
  # revert the values
  sam1Dos[sam12missing] = sam1val
  sam2Dos[sam12missing] = sam2val
 # return
 return (nVar,dist) 

def assertArgsAreGood(args):
 """ Assert that the command line arguments are within expected parameters.
 """
 if ( args.localCorrectionBP < 0 ):
  raise BaseError("The number of BP for local correction is 0 at minimum.")

 if ( args.minFrequency < 1e-9 or args.minFrequencyCorrect < 1e-9):
  raise BaseError("The minimum frequency must be >0")

 if ( args.maxFrequency > 1.0 - 1e-9 or args.maxFrequencyCorrect > 1.0-1e-9 ):
  raise BaseError("The maximum frequency must be <1")

def writeGRMID(reader,args):
 """ Write out the mapping from individual index to sample ID.
 """
 grmIDOut = open(args.grm+".id",'w')
 idx = 0
 for sample in PlinkReader.getSamples(reader): 
  grmIDOut.write("%s\t%s\n" % (sample.family_id,sample.individual_id))
  idx += 1

def runMain():
 # get the arguments
 args = parseArgs()
 assertArgsAreGood(args)
 # get the reader
 if ( args.opt ):
  bedReader = PlinkReader.SiteOptimizedPlinkBinaryReader(args.bedBase)
 else:
  bedReader = PlinkReader.getReader(args.bedBase)
 if ( not bedReader.snpMajor() ):
  raise BaseError("This tool currently only works with bed files in SNP-major mode. Please transpose the bed file.")
 writeGRMID(bedReader,args)
 # step 1: generate a locally-corrected per-sample dosages
 dosages,missingByIndividual = streamAndLocalCorrect(args,bedReader)
 # step 2: generate globally-corrected per-sample dosages
 if ( args.globalCorrection != None ):
  raise NotImplementedError("Global correction not yet fully tested.")
  dosages = globallyCorrectDosages(args,dosages,missingByIndividual)
 # step 3: calculate the GRM and print it out
 #printGRMFromDosages(args,dosages,missingByIndividual)


if __name__ == "__main__":
 #runMain()
 cProfile.run('runMain()')

