require("ggplot2")
require("randomForest")

args = commandArgs(TRUE);
onCommandLine <- ! is.na(args[1])

if ( onCommandLine ) pdf(height=8.5, width=11)

# variables controlling the system behavior:
BIG_MEMORY <- grepl("gsa", Sys.info()[4])
CHR20_ONLY = F
ENABLE_CACHING = F
EXPLORE_EXAMPLES = F
ANALYZE_BY_N_TRAINING_SITES <- F
ANALYZE_BY_N_ANNOTATIONS <- F
ANALYZE_BY_N_TREES <- T
ANALYZE_ROBUSTNESS_TO_NOISE_ANNOTATIONS <- F
ANALYZE_SENSITIVITY_TO_NEGATIVE_TRAINING_SET <- F
ANALYZE_NA_TREATMENT <- F
MAX_POLY_SITES_TO_EVAL = 1000000
DEFAULT_MAX_NEG_TRAINING_FRACTION = 0.25
nRocPoints = 1000
FORCE_RELOAD = F
MAX_ROWS_TO_READ = 1000000

N_TREES <- 5000
#N_TREES <- 100

dataDir = "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/VQSRv2"
setwd(dataDir)

USE_EXACT_FALSE_POSITIVES = T
indelAll <- read.table("/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/falsePositivesFromOmni/indel.chr20.pos.txt", col.names=c("POS"))
indel5 <- read.table("/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/falsePositivesFromOmni/indel5.chr20.pos.txt", col.names=c("POS"))
mono <- read.table("/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/falsePositivesFromOmni/NOT_POLY.chr20.pos.txt", col.names=c("POS"))
phaseIValidationFailure <- read.table("/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/falsePositivesFromOmni/phaseI.validation.failed.sites", header=T)
pilotValidationFailure <- read.table("/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/falsePositivesFromOmni/pilot.validation.failed.sites", col.names=c("POS"))

LABEL = "TrainingLabel"
LABEL_POSITIVE = "positive"
POLY_LABELS = c("POLY", "OmniPoly", "MillsPoly", "Phase1Poly", "DoubleHitPoly", "GoldStandardPoly")

moreLabels <- c(LABEL, "status", "FILTER", "VQSLOD")

addDerivedColumn <- function(data, field, columns) {
  values = rep("unknown", dim(data)[1])
  for ( col in columns ) {
    indices = !is.na(data[[col]])
    values[indices] = col
  }
  data[[field]] <- factor(values)
  data
}

readVQSRData <- function(variable, file, snps) {
  if ( T | ! exists(variable) ) {
    data = read.table(file, header=T, nrows=MAX_ROWS_TO_READ)
#    data <- addDerivedColumns(data)
    if ( snps ) 
      data = addDerivedColumn(data, "status", c("OmniPoly", "OmniMono", "Phase1Mono", "PilotMono", "BadTDT"))
    else
      #data = addDerivedColumn(data, "status", c("DoubleHitPoly", "MillsMono", "MillsPoly", "Phase1Mono", "Phase1Poly", "BadTDT", "GoldStandardPoly"))
      data = addDerivedColumn(data, "status", c("Phase1Mono", "Phase1Poly", "BadTDT", "GoldStandardPoly"))
  }
}

goNL.snps <- readVQSRData("goNL.snps", "/humgen/gsa-hpprojects/GATK/data/vqsrGrandChallenge/GoNL.recal.chr20.snp.dat", T)
goNL.indels <- readVQSRData("goNL.indels", "/humgen/gsa-hpprojects/GATK/data/vqsrGrandChallenge/GoNL.recal.chr20.indels.dat", F)
Autism.snps <- readVQSRData("Autism.snps", "/humgen/gsa-hpprojects/GATK/data/vqsrGrandChallenge/Autism.recal.exome.snp.dat", T)
Autism.indels <- readVQSRData("Autism.indels", "/humgen/gsa-hpprojects/GATK/data/vqsrGrandChallenge/Autism.recal.exome.indels.dat", F)

# FIXME -- GLOBAL TRAINING AND EVAL SITES
DATA.SET = "goNL.snps"
#DATA.SET = "goNL.indels"
#DATA.SET = "Autism.snps"
#DATA.SET = "Autism.indels"
ALL_SITES <- get(DATA.SET)
TRAINING_SITES = ALL_SITES
EVAL_SITES = subset(ALL_SITES, ! is.na(status) & status != "unknown")

originalAnn <- c("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "FS", "MQ","InbreedingCoeff", "DP")

if (grepl("snps", DATA.SET)) { 
  trainingAnn <- c("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "FS", "MQ", 
  	                    "InbreedingCoeff", "DP", "AC", "QUAL", "MQ0", "BaseQRankSum")
  #trainingAnn <- c("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "FS", "MQ", 
  # 	            "InbreedingCoeff", "DP", "AC", "QUAL", "MQ0", "BaseQRankSum", "SB")
} else {
  trainingAnn <- c("QD", "ReadPosRankSum", "FS", "InbreedingCoeff",
  "AC", "QUAL")
}
print(paste("Training annotations", trainingAnn))

roc <- function(nChunks, data, scoreLabel, name) {
  df <- data.frame()

  nAllPoly <- sum(data$status %in% POLY_LABELS)
  nAllMono <- sum(!(data$status %in% POLY_LABELS))
  
  scores <- data[[scoreLabel]]
  sorted <- data[order(scores, decreasing=T),]
  
  n <- dim(data)[1]
  chunks = round(c(seq(1, n, n / nChunks), n))
  
  lastChunk = 0
  lastPoly = 0
  lastMono = 0
  for ( chunk in chunks ) {
    #print(list(last=lastChunk, c=chunk))
    sub = sorted[(lastChunk+1):chunk,]

    nPoly <- lastPoly + sum(sub$status %in% POLY_LABELS)
    nMono <- lastMono + sum(!(sub$status %in% POLY_LABELS))
    
    sensitivity = nPoly / nAllPoly
    specificity = 1 - nMono / nAllMono
    one <- data.frame(sensitivity = sensitivity, specificity = specificity, name = name, nTrees = NA, nTrainingSites=NA, nAnn = 8)
    df <- rbind(df, one)
    
    # update
    lastChunk = chunk
    lastPoly = nPoly
    lastMono= nMono
  }
  
  df
}

if ( ! exists("vqsr.roc") & exists("EVAL_SITES") ) 
  vqsr.roc <- roc(nRocPoints, EVAL_SITES, "VQSLOD", "VQSR")

trainTree <- function(name, nTrees, nTrainingSites, trainingAnn, maxNegTrainingFraction = DEFAULT_MAX_NEG_TRAINING_FRACTION, na.func=na.roughfix) {
  if ( LABEL == "TrainingLabel" ) {
    trainingData <- subset(TRAINING_SITES, TrainingLabel != "neutral")
  } else {
    trainingData <- subset(TRAINING_SITES, status != "UNKNOWN")
  }  
  
  trainingData <- trainingData[, c(trainingAnn, moreLabels)]
  trainingData <- na.func(trainingData)

  # if required, take only the worst maxNegTrainingFraction of the negative training data 
  if ( maxNegTrainingFraction != -1 ) {
    print(table(trainingData$TrainingLabel))
    pos = subset(trainingData, TrainingLabel == "positive")
    neg = subset(trainingData, TrainingLabel == "negative")
    nNeg = dim(neg)[1]
    nToKeep = round(nNeg * maxNegTrainingFraction)
    keep = order(neg$VQSLOD)[1:nToKeep]
    negToKeep = neg[keep,]
    trainingData = rbind(pos, negToKeep)
    print(table(trainingData$TrainingLabel))
  }
  
  # If requested, sample only nTrainingSites from the total training set
  if ( nTrainingSites != -1 & nTrainingSites < dim(trainingData)[1]) {
    indices = sample(1:dim(trainingData)[1], nTrainingSites)
    trainingData <- trainingData[indices,]
    print(table(trainingData$TrainingLabel))
  }
  trainingData[[LABEL]] <- factor(trainingData[[LABEL]]) # refactor the label, so that missing levels are eliminated

  # create the name of the tree, and if it already exists on disk read it in and use it instead
  nActuaTrainingSites = dim(trainingData)[1]
  annotationString = do.call("paste", c(sep="_", as.list(trainingAnn)))
  name = paste(sep=".", name, paste(sep="_", "nTrees", nTrees, "nActualTrainingSites", nActuaTrainingSites, "annotations", annotationString), "tree")
  print(name)
  
  if ( ENABLE_CACHING & file.exists(name) ) {
    print(paste("Tree cached on disk, loading", name))
    loadTree(name)
  } else {
    print(paste("Building tree", name))
    # actually do the work to build the tree, and write it out to disk
    training.urf <- randomForest(x=trainingData[,trainingAnn], 
                                 y=trainingData[[LABEL]], 
                                 importance=TRUE, 
                                 proximity=F,
                                 ntree=nTrees,
                                 keep.forest=TRUE,
                                 do.trace=100)
    l = list(name=name, rf=training.urf, trainingAnn = trainingAnn, nTrees = nTrees, nTrainingSites = nTrainingSites)
    save(l, file=name)
    loadTree(name)
  }
}

loadTree <- function(filename) {
  env = new.env()
  load(file=filename, env=env)
  get(ls(env), env)
}

rocForTree <- function(tree) {
  noNA <- na.omit(EVAL_SITES[,c(tree$trainingAnn, "VQSLOD", "FILTER", "TrainingLabel", "status")]) #, "EVENTLENGTH")])
  #omni.pred <- predict(tree$rf, noNA)
  #print(table(observed = noNA$status, predicted = omni.pred))
  omni.prob <- predict(tree$rf, noNA, type="prob")
  omni.with.prob <- cbind(noNA, omni.prob)
  
  myRoc = roc(nRocPoints, omni.with.prob, LABEL_POSITIVE, tree$name)
  myRoc$nTrees <- tree$nTrees
  #annotationString = do.call("paste", c(sep=",", as.list(tree$trainingAnn)))
  #roc$annotationString <- annotationString
  myRoc$nTrainingSites <- tree$nTrainingSites
  myRoc$nAnn <- length(tree$trainingAnn)
  myRoc
}

importancePlot <- function(rf) {
  df <- melt(tree.small$rf$importance) 
  names(df) <- c("RF.variable", "Importance.measure", "value")
  p <- ggplot(data=df, aes(y=reorder(RF.variable, value), x=value, color=Importance.measure))
  p <- p + facet_grid(. ~ Importance.measure, scale="free")
  p <- p + geom_point(size=5) #  + geom_linerange(aes(ymin=0, ymax=value), size=2)
#  p <- p + coord_flip()
  p <- p + xlab("Score") + ylab("Random forest training variable")
  print(p)
}

plotRocs <- function(rocs) {
  #print(rocs)
  rocs$shortName = substr(rocs$name, 1, 20)
  rocs$isVQSR = rocs$shortName == "VQSR"
  p <- ggplot(data=rocs, aes(y=sensitivity, x=1-specificity, group=shortName, color=shortName, size=isVQSR))
  p <- p + opts(title=paste("ROC", DATA.SET))  
  #p <- p + geom_point()
  p <- p + geom_line() # + scale_x_log10() + scale_y_log10()

  print(p + xlim(0.0, 1.0) + ylim(0.0, 1.01))
  print(p + xlim(0.0, 0.5) + ylim(0.9, 1.01))
  print(p + xlim(0.0, 0.1) + ylim(0.9, 1.01))

  p
}

if ( EXPLORE_EXAMPLES ) {
  #N_TREES = 1000
  tree.small <- trainTree("tree.small", N_TREES, 500, trainingAnn)
  tree.medium <- trainTree("tree.medium", N_TREES, 5000, trainingAnn)
  tree.big <- trainTree("tree.big", N_TREES, 500000, trainingAnn)
  tree.big.original <- trainTree("tree.vqsr.ann", N_TREES, 500000, originalAnn)
  
  trees <- list(tree.small, tree.medium, tree.big, tree.big.original)#, tree.size)

  rocs = do.call("rbind", lapply(trees, rocForTree))
  plotRocs(rbind(vqsr.roc, rocs))  
}

standardSensSpecTargets <- function(rocs) {
  SENSITIVITY_TARGETS = c(0.97, 0.98, 0.99)
  SPECIFICITY_TARGETS = c(0.95, 0.98, 0.99)

  at.sensitivities = data.frame()
  for ( target in SENSITIVITY_TARGETS ) {
    tmp = subset(rocs, sensitivity >= target) # ddply(rocs, .variables=c("name"), subset, specificity >= target)
    at.sensitivity = ddply(tmp, .variables=c("name"), subset, order(sensitivity) == 1)
    at.sensitivity$target = target
    at.sensitivities = rbind(at.sensitivity, at.sensitivities)
  }
  at.sensitivities$target <- factor(at.sensitivities$target)
  
  at.fdrs = data.frame()
  for ( target in SPECIFICITY_TARGETS ) {
    tmp = subset(rocs, specificity >= target) # ddply(rocs, .variables=c("name"), subset, specificity >= target)
    at.fdr = ddply(tmp, .variables=c("name"), subset, order(sensitivity, decreasing=T) == 1)
    at.fdr$target = target
    at.fdrs = rbind(at.fdr, at.fdrs)
  }
  at.fdrs$target <- factor(at.fdrs$target)
  
  list(at.sensitivity=at.sensitivities, at.fdrs=at.fdrs)
}
#standardSensSpecTargets(rocs)


# ------------------------------------------------------------------------------------------
#
# Performance of classifier by N training sites
#
# ------------------------------------------------------------------------------------------
byNTrainingSites <- function(arg = 1) {
  rocs = data.frame()
  for ( nTrainingSites in c(100, 500, 1000, 2500, 5000, 10000, 50000, 100000, 250000, 500000, 750000, 1000000) ) {
  #for ( nTrainingSites in c(10, 100, 500, 1000, 2000, 3000, 5000, 7500, 10000, 20000) ) {
    name = paste("tree.", nTrainingSites, sep="")
    # with 1M training sites we can only build 100 trees
    tree = trainTree(name, 100, nTrainingSites, trainingAnn)
    trees = c(list(tree), trees)
    rocs = rbind(rocs, rocForTree(tree))
  }
  plotRocCuts("nTrainingSites", rocs)
  rocs
}

distributeGraphRows <- function(graphs, heights = c()) {
  if (length(heights) == 0) {
    heights <- rep.int(1, length(graphs))
  }
  heights <- heights[!is.na(graphs)]
  graphs <- graphs[!is.na(graphs)]
  numGraphs <- length(graphs)
  Layout <- grid.layout(nrow = numGraphs, ncol = 1, heights=heights)
  grid.newpage()
  pushViewport(viewport(layout = Layout))
  subplot <- function(x) viewport(layout.pos.row = x, layout.pos.col = 1)
  for (i in 1:numGraphs) {
    print(graphs[[i]], vp = subplot(i))
  }
}

plotRocCuts <- function(xfield, rocs, log=T) {
  addInfo = function(title, p) {
    p = p + geom_point(size=3) + geom_line() + geom_smooth()
    p = p + opts(title=title)
    if ( log) p = p + scale_x_log10()
    p = p + xlab(xfield)
    p
  }
  
  x = standardSensSpecTargets(rocs)
  t = x$at.sensitivity
  p = ggplot(data=data.frame(x=t[[xfield]], specificity=t$specificity, target=t$target), 
             aes(x=x, y=specificity, group=target, color=target))
  p.sensitivity = addInfo("Specificity at given sensitivity", p)

  t = x$at.fdrs
  p = ggplot(data=data.frame(x=t[[xfield]], sensitivity=t$sensitivity, target=t$target), 
             aes(x=x, y=sensitivity, group=target, color=target))
  p.specificity =addInfo("Sensitivity at given specificity", p)
  
  distributeGraphRows(list(p.sensitivity, p.specificity), c(1,1))
}

if ( ANALYZE_BY_N_TRAINING_SITES ) {
  x <- byNTrainingSites()
  plotRocs(rbind(vqsr.roc, x))
}


# ------------------------------------------------------------------------------------------
#
# assessing sensitivity to no. good and noise annotations
#
# ------------------------------------------------------------------------------------------

# in order of mean decrease in accuracy
good.anns.in.order <- c("MQ0", 
                       "QD", 
                       "MQ", 
                       "InbreedingCoeff",
                       "QUAL",
                       "SB",
                       "MQRankSum",
                       "HaplotypeScore", 
                       "AC",
                       "DP", 
                       "ReadPosRankSum")
noise.anns <- c("BaseQRankSum", 
               "Dels", 
               "HRun",
               "TRANSITION", 
               "AN", 
               "POS", 
               "REF", 
               "ALT")

byGoodAnnotations <- function(nTrainingSites, goodAnnotationsInOrder) {
  trees = list()
  for ( stopI in 2:length(goodAnnotationsInOrder) ) {
    anns = goodAnnotationsInOrder[1:stopI]
    name = do.call("paste", c(sep="+", "tree", as.list(anns)))
    #print(anns)
    print(name)
    tree = trainTree(name, N_TREES, nTrainingSites, anns)
    trees = c(list(tree), trees)
  }
  rocs = do.call("rbind", lapply(trees, rocForTree))
  plotRocCuts("nAnn", rocs, log=F)
  list(trees=trees, rocs=rocs)
}

if ( ANALYZE_BY_N_ANNOTATIONS ) {
  by.good.ann <- byGoodAnnotations(10000, good.anns.in.order)
}

byNoiseAnnotations <- function(nTrainingSites, noiseAnnotations) {
  trees = list()
  for ( stopI in 0:length(noiseAnnotations) ) {
    noiseAnns = noiseAnnotations[0:stopI]
    nNoise = length(noiseAnns)
    name = paste("tree.with.noise.n=", nNoise, sep="")
    print(name)
    tree = trainTree(name, N_TREES, nTrainingSites, c(good.anns.in.order, noiseAnns))
    tree$nNoiseAnn = nNoise
    trees = c(list(tree), trees)
  }
  
  rocForTreeWithNoiseCount <- function(tree) {
    roc = rocForTree(tree)
    roc$nNoiseAnn = tree$nNoiseAnn
    roc
  }
  
  rocs = do.call("rbind", lapply(trees, rocForTreeWithNoiseCount))
  plotRocCuts("nNoiseAnn", rocs, log=F)
  list(trees=trees, rocs=rocs)
}

if ( ANALYZE_ROBUSTNESS_TO_NOISE_ANNOTATIONS ) {
  by.noise.ann <- byNoiseAnnotations(10000, noise.anns)
}

# ------------------------------------------------------------------------------------------
#
# By severity of negative training sites
#
# ------------------------------------------------------------------------------------------

byNegTrainingFraction <- function() {
  trees = list()
  #for ( percentNegative in c(0.05, 0.1, 0.25, 0.5, 0.75, 1.00) ) {
  for ( percentNegative in c(0.05, 0.25, 1.00) ) {
    name = paste("tree.", percentNegative, sep="")
    tree = trainTree(name, 100, -1, trainingAnn, maxNegTrainingFraction=percentNegative)
    trees = c(list(tree), trees)
  }
  print(trees)
  rocs = do.call("rbind", lapply(trees, rocForTree))
  plotRocs(rbind(vqsr.roc, rocs))
  #plotRocCuts("nTrainingSites", rocs)
  list(trees=trees, rocs=rocs)
}

if ( ANALYZE_SENSITIVITY_TO_NEGATIVE_TRAINING_SET ) {
  by.neg.training.fraction <- byNegTrainingFraction()
}

#
# performance of new roc calculation
# 
# MAG = 10
# big <- data.frame()
# for ( i in 1:MAG ) big <- rbind(big, EVAL_SITES)
# print(dim(big))
# vqsr.mar.roc <- rocMarginal(nRocPoints, big, "VQSLOD", "VQSR.mar")
# vqsr.brute.roc <- rocOriginal(nRocPoints, big, "VQSLOD", "VQSR.brute")
# plotRocs(rbind(vqsr.brute.roc, vqsr.mar.roc))

#tree1 = trainTree("foo", 500, 100000, trainingAnn)
#tree2 = trainTree("foo", 1000, 10000, trainingAnn)

na.fix_with_imputation <- function(df) {
  # shockingly expensive
  rfImpute(x=df[,trainingAnn], y=df[[LABEL]])
}

byNATreatment <- function() {
  trees = list()
  for ( na.func in c("na.omit", "na.roughfix") ) {
    name = paste("tree.", na.func, sep="")
    tree = trainTree(name, N_TREES, -1, trainingAnn, na.func=get(na.func)) # , maxNegTrainingFraction=1)
    trees = c(list(tree), trees)
  }

  rocs = do.call("rbind", lapply(trees, rocForTree))
  plotRocs(rbind(vqsr.roc, rocs))
  list(trees=trees, rocs=rocs)
}

if ( ANALYZE_NA_TREATMENT ) {
  by.na.treatment <- byNATreatment()
}

#
# Performance vs. N_TREES
#
byNTrees <- function() {
  trees = list()
  ns <- c(50, 10)
  ns <- c(2000, 1000, 500, 100, ns)
  #ns <- c(8000, 4000, ns)
  #ns <- c(32000, 24000, 16000, ns)
  for ( nTrees in ns ) { 
    name = paste("tree.", nTrees, sep="")
    tree = trainTree(name, nTrees, -1, trainingAnn) # , maxNegTrainingFraction=1)
    trees = c(list(tree), trees)
  }

  rocs = do.call("rbind", lapply(trees, rocForTree))
  plotRocs(rbind(vqsr.roc, rocs))
  list(trees=trees, rocs=rocs)
}

if ( ANALYZE_BY_N_TREES ) {
  by.n.trees <- byNTrees()
}


# must be at the end
if ( onCommandLine ) dev.off()

