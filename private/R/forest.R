require("ggplot2")
require("randomForest")

args = commandArgs(TRUE);
onCommandLine <- ! is.na(args[1])

# variables controlling the system behavior:
BIG_MEMORY <- grepl("gsa", Sys.info()[4])
CHR20_ONLY = F
EXPLORE_EXAMPLES = F
ANALYZE_BY_N_TRAINING_SITES <- T
ANALYZE_BY_N_ANNOTATIONS <- F
ANALYZE_ROBUSTNESS_TO_NOISE_ANNOTATIONS <- F
ANALYZE_SENSITIVITY_TO_NEGATIVE_TRAINING_SET <- F
MAX_POLY_SITES_TO_EVAL = 1000000
DEFAULT_MAX_NEG_TRAINING_FRACTION = 0.25
nRocPoints = 1000
FORCE_RELOAD = F
MAX_ROWS_TO_READ = -1 # 100000

N_TREES <- 1000
#N_TREES <- 100

dataDir = "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/VQSRv2"
setwd(dataDir)

USE_EXACT_FALSE_POSITIVES = T
indelAll <- read.table("/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/falsePositivesFromOmni/indel.chr20.pos.txt", col.names=c("POS"))
indel5 <- read.table("/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/falsePositivesFromOmni/indel5.chr20.pos.txt", col.names=c("POS"))
mono <- read.table("/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/falsePositivesFromOmni/NOT_POLY.chr20.pos.txt", col.names=c("POS"))
phaseIValidationFailure <- read.table("/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/falsePositivesFromOmni/phaseI.validation.failed.sites", header=T)
pilotValidationFailure <- read.table("/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/falsePositivesFromOmni/pilot.validation.failed.sites", col.names=c("POS"))

originalAnn <- c("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "FS", "MQ","InbreedingCoeff", "DP")
standardAnn <- c("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "FS", "MQ", 
                 "InbreedingCoeff", "DP", "AC", "QUAL", "MQ0", "BaseQRankSum", "SB")
# standardAnn <- c("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "FS", "MQ", 
#                  "InbreedingCoeff", "DP", "AC", "QUAL", "MQ0", "BaseQRankSum", "Dels", "HRun",
#                  "SB", "pop", "TRANSITION", "AN", "POS", "known", "REF", "ALT")
trainingAnn <- standardAnn
#trainingAnn <- c("QD", "FS", "InbreedingCoeff")

addDerivedColumns <- function(data) {
  data$known <- factor(data$ID == ".", labels=c("known", "novel"))
  
  data$status <- "UNKNOWN"
  ONLY_EXTRAS = F
  if ( ! ONLY_EXTRAS )
    data$status[! is.na(data$OmniMono)] <- "MONO" 
  data$status[! is.na(data$OmniPoly)] <- "POLY"
  
  if ( USE_EXACT_FALSE_POSITIVES ) {
    #data$status[data$POS %in% indel5$POS] <- "INDEL"
    #data$status[data$POS %in% mono$POS] <- "MONO"
    data$status[data$POS %in% phaseIValidationFailure$POS] <- "PhaseIValidationMono"
    data$status[data$POS %in% pilotValidationFailure$POS] <- "PilotValidationMono"
  }
  
  data$status <- factor(data$status, unique(data$status))
  data
}

if ( ! exists("allSites") ) {
  print("Loading allSites")
  if ( CHR20_ONLY ) {
    allSites <- read.table("combined.phase1.wgs.recal.snps.20.dat", header=T, nrows=MAX_ROWS_TO_READ)
  } else {
    allSites <- read.table("combined.phase1.wgs.recal.snps.training.dat", header=T, nrows=MAX_ROWS_TO_READ)
  }
  allSites <- addDerivedColumns(allSites)
}

if ( ! exists("omniSites") ) {
  print("Loading omniSites")
  omniSites <- read.table("combined.phase1.wgs.recal.snps.omni.dat", header=T, nrows=MAX_ROWS_TO_READ)
  omniSites <- addDerivedColumns(omniSites)

  nonOmniMono = subset(allSites, status != "POLY" & status != "UNKNOWN" & status != "MONO")
  omniSites <- rbind(omniSites, nonOmniMono)

  omniSites <- na.omit(omniSites[,c(trainingAnn, "VQSLOD", "FILTER", "TrainingLabel", "status")])
  
  nPoly = sum(omniSites$status == "POLY")
  if ( nPoly > MAX_POLY_SITES_TO_EVAL ) {
    poly = subset(omniSites, status == "POLY")
    polysub = sample(1:nPoly, MAX_POLY_SITES_TO_EVAL)
    omniSites = rbind(poly[polysub,], subset(omniSites, status != "POLY" ))
  }
}

#LABEL = "status"
#LABEL_POSITIVE = "POLY"

LABEL = "TrainingLabel"
LABEL_POSITIVE = "positive"

#omniSites <- subset(all, status != "UNKNOWN")
moreLabels <- c(LABEL, "status", "known", "FILTER", "VQSLOD")

roc <- function(nChunks, data, scoreLabel, name) {
  df <- data.frame()

  nAllPoly <- sum(data$status == "POLY")
  nAllMono <- sum(data$status != "POLY")
  
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

    nPoly <- lastPoly + sum(sub$status == "POLY")
    nMono <- lastMono + sum(sub$status != "POLY")
    
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

trainTree <- function(name, nTrees, nTrainingSites, trainingAnn, maxNegTrainingFraction = DEFAULT_MAX_NEG_TRAINING_FRACTION) {
  if ( LABEL == "TrainingLabel" ) {
    trainingData <- subset(allSites, TrainingLabel != "neutral")
  } else {
    trainingData <- subset(allSites, status != "UNKNOWN")
  }  
  
  trainingData <- trainingData[, c(trainingAnn, moreLabels)]
  trainingData <- na.omit(trainingData)

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
  
  if ( file.exists(name) ) {
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
                                 do.trace=10)
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
  omni.pred <- predict(tree$rf, omniSites)
  print(table(observed = omniSites$status, predicted = omni.pred))
  omni.prob <- predict(tree$rf, omniSites, type="prob")
  omni.with.prob <- cbind(omniSites, omni.prob)
  
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

if ( EXPLORE_EXAMPLES ) {
  #tree.all <- trainTree("tree.all", N_TREES, -1, trainingAnn)
  N_TREES = 20
  tree.small <- trainTree("tree.500", N_TREES, 500, trainingAnn)
  tree.medium <- trainTree("tree.5000", N_TREES, 5000, trainingAnn)
  tree.big <- trainTree("tree.-1", N_TREES, 500000, trainingAnn)
  tree.big.original <- trainTree("tree.-1.vqsr.ann", N_TREES, 500000, originalAnn)
  
  trees <- list(tree.small, tree.medium, tree.big, tree.big.original) #, tree.big) # , tree.all)

  #importancePlot(tree.big$rf)

  trees.roc <- lapply(trees, rocForTree)
  rocs <- roc(nRocPoints, omniSites, "VQSLOD", "VQSR")
  for ( myRoc in trees.roc )
    rocs <- rbind(rocs, myRoc)
  plotRocs(rocs)
}

plotRocs <- function(rocs) {
  #print(rocs)
  p <- ggplot(data=rocs, aes(y=sensitivity, x=1-specificity, group=name, color=name))
  p <- p + geom_point() + geom_line() # + scale_x_log10() + scale_y_log10()
  p <- p + xlim(0.0, 0.5) + ylim(0.9, 1.01)
  print(p)
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

if ( ! exists("vqsr.roc") ) 
  vqsr.roc <- roc(nRocPoints, omniSites, "VQSLOD", "VQSR")

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
# for ( i in 1:MAG ) big <- rbind(big, omniSites)
# print(dim(big))
# vqsr.mar.roc <- rocMarginal(nRocPoints, big, "VQSLOD", "VQSR.mar")
# vqsr.brute.roc <- rocOriginal(nRocPoints, big, "VQSLOD", "VQSR.brute")
# plotRocs(rbind(vqsr.brute.roc, vqsr.mar.roc))

#tree1 = trainTree("foo", 500, 100000, trainingAnn)
#tree2 = trainTree("foo", 1000, 10000, trainingAnn)
