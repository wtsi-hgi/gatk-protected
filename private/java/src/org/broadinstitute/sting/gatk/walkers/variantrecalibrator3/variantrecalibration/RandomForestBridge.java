package org.broadinstitute.sting.gatk.walkers.variantrecalibrator3.variantrecalibration;

import org.broadinstitute.sting.utils.R.RScriptExecutor;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 3/2/12
 * Time: 12:12 PM
 * To change this template use File | Settings | File Templates.
 */
public class RandomForestBridge extends RecalibrationModel {

    protected String intermediateFile;
    protected List<String> annotations;

    //private static final String[] HEADER_START = {"CHR","POS","TYP","NAL"};
    private static final String[] HEADER_START = {"CHR","POS","TYP"};
    private static final String[] HEADER_END = {"TRAINING","KNOWN","status","vqslod"};

    public RandomForestBridge(boolean isFinal, VariantDataManager manager) {
        super(isFinal);
        annotations = isFinal ? manager.finalKeys : manager.initialKeys;
    }

    public void initialize(List<VariantDatum> data, VariantRecalibratorArgumentCollection args) {
        intermediateFile = args.ROutputFile;
        File fileForROutput = write(data,args.RInputFile, intermediateFile, args.RScriptFile, args.numTree);
        RScriptExecutor executor = new RScriptExecutor();
        executor.setExceptOnError(true);
        executor.addScript(new File(args.RScriptFile));
        logger.info("Executing: "+executor.getApproximateCommandLine());
        executor.exec();
        // output lines up one-to-one
        XReadLines xrl;
        try {
            xrl = new XReadLines(fileForROutput);
        } catch (FileNotFoundException e) {
            throw new StingException("file not found", e);
        }
        String line = xrl.next();
        for ( VariantDatum d : data ) {
            if ( useFinal )
                d.setFinalLod(Double.parseDouble(line));
            else
                d.setInitialLod(Double.parseDouble(line));
            line = xrl.hasNext() ? xrl.next() : "NO DATA";
        }
    }


    private File write(List<VariantDatum> data, String oFileStr, String rOutput, String scriptFileString, String num_tree) {
        File oFile = new File(oFileStr);
        try {
            PrintStream out = new PrintStream(oFile);
            List<String> header = new ArrayList<String>(annotations.size()+HEADER_END.length+HEADER_START.length);
            header.addAll(Arrays.asList(HEADER_START));
            header.addAll(annotations);
            header.addAll(Arrays.asList(HEADER_END));
            out.printf("%s%n", Utils.join("\t",header));
            for ( VariantDatum d : data ) {
                StringBuilder bd = new StringBuilder();
                bd.append(d.contig);
                bd.append("\t");
                bd.append(d.start);
                bd.append("\t");
                bd.append(d.isSNP ? "SNP" : "INDEL");
                bd.append("\t");
                double[] annot = useFinal ? d.finalAnnotations : d.initialAnnotations;
                for ( double dd : annot ) {
                    bd.append(dd);
                    bd.append("\t");
                }
                bd.append(d.atPositiveTrainingSite ? "T" : ( d.atNegativeTrainingSite ? "F" : "U" ) );
                bd.append("\t");
                bd.append(d.isKnown ? "K" : "N");
                bd.append("\t");
                bd.append( d.atTruthSite ? "T" : ( d.atMonomorphicSite ? "F" : "U"));
                bd.append("\t");
                bd.append(d.getLod());
                out.printf("%s%n",bd.toString());
            }
        } catch (FileNotFoundException e) {
            throw new StingException("File not found",e);
        }

        // write the script
        try {
            PrintStream script = new PrintStream(new File(scriptFileString));
            script.printf("require(\"ggplot2\")\n" +
                    "require(\"randomForest\")\n" +
                    "library(foreach)\n" +
                    "library(doMC)\n" +
                    "\n" +
                    "inputFile = \""+oFileStr+"\"\n" +
                    "outputFile = \""+rOutput+"\"\n" +
                    "\n" +
                    "N_TREE = "+num_tree+"\n" +
                    "MAX_LOD = 6.04\n" +
                    "\n" +
                    "roc <- function(nChunks,classifiedData,scoreLabel,name) {\n" +
                    "\t# extract the positive and the orthogonal negative points\n" +
                    "\t#print(\"DEBUG 55\")\n" +
                    "\tdata <- subset(classifiedData, ! is.na(status) & status !=\"U\")\n" +
                    "\t#print(data)\n" +
                    "\tscores <- factor(data[[scoreLabel]])\n" +
                    "\tsorted <- data[order(scores,decreasing=T),]\n" +
                    "\tn <- dim(data)[1]\n" +
                    "\tchunks = round(c(seq(1,n,n/nChunks),n))\n" +
                    "\tlastChunk = 0\n" +
                    "\tlastPos = 0\n" +
                    "\tlastNeg = 0\n" +
                    "\tnAllPos = sum(data$status==\"T\")\n" +
                    "\t#print(\"NALLPOS\")\n" +
                    "\t#print(nAllPos)\n" +
                    "\tnAllNeg = sum(data$status==\"F\")\n" +
                    "\t#print(\"NALLNEG\")\n" +
                    "\t#print(nAllNeg)\n" +
                    "\tdf <- data.frame()\n" +
                    "\tfor ( chunk in chunks ) {\n" +
                    "\t\tsub = sorted[(lastChunk+1):chunk,]\n" +
                    "\t\tnPos <- lastPos + sum(sub$status==\"T\")\n" +
                    "\t\tnNeg <- lastNeg + sum(sub$status==\"F\")\n" +
                    "\t\tsensitivity = nPos/nAllPos\n" +
                    "\t\tspecificity = 1 - nNeg/nAllNeg\n" +
                    "\t\tone <- data.frame(sensitivity = sensitivity, specificity = specificity, name = name, nTrees = NA, nTrainingSites=NA, nAnn = 8)\n" +
                    "\t\tdf <- rbind(df,one)\n" +
                    "\t\tlastChunk = chunk\n" +
                    "\t\tlastPos = nPos\n" +
                    "\t\tlastNeg = nNeg\t\n" +
                    "\t}\n" +
                    "\t\n" +
                    "\tdf\t\n" +
                    "}\n" +
                    "\n" +
                    "plotRocs <- function(rocs, title) {\n" +
                    "  #print(rocs)\n" +
                    "  rocs$shortName = substr(rocs$name, 1, 20)\n" +
                    "  rocs$isVQSR = grepl(\"VQSR\", rocs$shortName)\n" +
                    "  p <- ggplot(data=rocs, aes(y=sensitivity, x=1-specificity, group=shortName, color=shortName, size=isVQSR))\n" +
                    "  p <- p + opts(title=paste(\"ROC\", title))  \n" +
                    "  #p <- p + geom_point()\n" +
                    "  p <- p + geom_line() # + scale_x_log10() + scale_y_log10()\n" +
                    "\n" +
                    "  print(p + xlim(0.0, 1.0) + ylim(0.0, 1.01))\n" +
                    "  print(p + xlim(0.0, 0.5) + ylim(0.9, 1.01))\n" +
                    "  print(p + xlim(0.0, 0.1) + ylim(0.9, 1.01))\n" +
                    "\n" +
                    "  p\n" +
                    "}\n" +
                    "\n" +
                    "standardSensSpecTargets <- function(rocs) {\n" +
                    "  SENSITIVITY_TARGETS = c(0.97, 0.98, 0.99)\n" +
                    "  SPECIFICITY_TARGETS = c(0.95, 0.98, 0.99)\n" +
                    "\n" +
                    "  at.sensitivities = data.frame()\n" +
                    "  for ( target in SENSITIVITY_TARGETS ) {\n" +
                    "    tmp = subset(rocs, sensitivity >= target) # ddply(rocs, .variables=c(\"name\"), subset, specificity >= target)\n" +
                    "    at.sensitivity = ddply(tmp, .variables=c(\"name\"), subset, order(sensitivity) == 1)\n" +
                    "    at.sensitivity$target = target\n" +
                    "    at.sensitivities = rbind(at.sensitivity, at.sensitivities)\n" +
                    "  }\n" +
                    "  at.sensitivities$target <- factor(at.sensitivities$target)\n" +
                    "  \n" +
                    "  at.fdrs = data.frame()\n" +
                    "  for ( target in SPECIFICITY_TARGETS ) {\n" +
                    "    tmp = subset(rocs, specificity >= target) # ddply(rocs, .variables=c(\"name\"), subset, specificity >= target)\n" +
                    "    at.fdr = ddply(tmp, .variables=c(\"name\"), subset, order(sensitivity, decreasing=T) == 1)\n" +
                    "    at.fdr$target = target\n" +
                    "    at.fdrs = rbind(at.fdr, at.fdrs)\n" +
                    "  }\n" +
                    "  at.fdrs$target <- factor(at.fdrs$target)\n" +
                    "  \n" +
                    "  list(at.sensitivity=at.sensitivities, at.fdrs=at.fdrs)\n" +
                    "}\n" +
                    "\n" +
                    "distributeGraphRows <- function(graphs, heights = c()) {\n" +
                    "  if (length(heights) == 0) {\n" +
                    "    heights <- rep.int(1, length(graphs))\n" +
                    "  }\n" +
                    "  heights <- heights[!is.na(graphs)]\n" +
                    "  graphs <- graphs[!is.na(graphs)]\n" +
                    "  numGraphs <- length(graphs)\n" +
                    "  Layout <- grid.layout(nrow = numGraphs, ncol = 1, heights=heights)\n" +
                    "  grid.newpage()\n" +
                    "  pushViewport(viewport(layout = Layout))\n" +
                    "  subplot <- function(x) viewport(layout.pos.row = x, layout.pos.col = 1)\n" +
                    "  for (i in 1:numGraphs) {\n" +
                    "    print(graphs[[i]], vp = subplot(i))\n" +
                    "  }\n" +
                    "}\n" +
                    "\n" +
                    "plotRocCuts <- function(xfield, rocs, log=T) {\n" +
                    "  addInfo = function(title, p) {\n" +
                    "    p = p + geom_point(size=3) + geom_line() + geom_smooth()\n" +
                    "    p = p + opts(title=title)\n" +
                    "    if ( log) p = p + scale_x_log10()\n" +
                    "    p = p + xlab(xfield)\n" +
                    "    p\n" +
                    "  }\n" +
                    "  \n" +
                    "  x = standardSensSpecTargets(rocs)\n" +
                    "  t = x$at.sensitivity\n" +
                    "  p = ggplot(data=data.frame(x=t[[xfield]], specificity=t$specificity, target=t$target), \n" +
                    "             aes(x=x, y=specificity, group=target, color=target))\n" +
                    "  p.sensitivity = addInfo(\"Specificity at given sensitivity\", p)\n" +
                    "\n" +
                    "  t = x$at.fdrs\n" +
                    "  p = ggplot(data=data.frame(x=t[[xfield]], sensitivity=t$sensitivity, target=t$target), \n" +
                    "             aes(x=x, y=sensitivity, group=target, color=target))\n" +
                    "  p.specificity =addInfo(\"Sensitivity at given specificity\", p)\n" +
                    "  \n" +
                    "  distributeGraphRows(list(p.sensitivity, p.specificity), c(1,1))\n" +
                    "}\n" +
                    "\n" +
                    "variantDataTable <- read.table(inputFile,header=TRUE)\n" +
                    "colNames <- colnames(variantDataTable)\n" +
                    "trainAnnot <- colNames[3:(length(colNames)-4)] # (chr,pos,type,alleles,{ATTRIBUTES},training,known,status,vqslod)\n" +
                    "# get only the labeled points; drop the known rod attribute\n" +
                    "trainingData <- subset(variantDataTable,TRAINING!=\"U\")\n" +
                    "# remove the empty \"U\" level\n" +
                    "response <- factor(trainingData[[\"TRAINING\"]])\n" +
                    "# drop out training and known columns\n" +
                    "trainingData <- na.roughfix(trainingData[,c(trainAnnot)])\n" +
                    "# now grab all the data\n" +
                    "classificationData <- na.roughfix(variantDataTable[,c(trainAnnot)])\n" +
                    "# train a random forest on this training data\n" +
                    "trainedRF <- randomForest(x = trainingData,\n" +
                    "                          y = response,\n" +
                    "                          importance = TRUE,\n" +
                    "                          proximity = F,\n" +
                    "                          ntree = N_TREE,\n" +
                    "                          keep.forest = FALSE,\n" +
                    "                          xtest = classificationData,\n" +
                    "                          norm.votes=TRUE,\n" +
                    "                          do.trace=100)\n" +
                    "classProbs <- trainedRF$test$votes\n" +
                    "d1 <- (classProbs[,c(\"T\")]>(1.0-10^(-MAX_LOD)))\n" +
                    "d2 <- (classProbs[,c(\"T\")]<10^(-MAX_LOD)) \n" +
                    "classProbs[d1,c(\"T\")] <- 1.0-10^(-MAX_LOD)\n" +
                    "classProbs[d1,c(\"F\")] <- 10^(-MAX_LOD)\n" +
                    "classProbs[d2,c(\"T\")] <- 10^(-MAX_LOD)\n" +
                    "classProbs[d2,c(\"F\")] <- 1.0-10^(-MAX_LOD)\n" +
                    "classLod <- log10(classProbs[,c(\"T\")]) - log10(classProbs[,c(\"F\")])\n" +
                    "classificationData$lod <- classLod\n" +
                    "classificationData$status <- variantDataTable$status\n" +
                    "classificationData$vqslod <- variantDataTable$vqslod\n" +
                    "rocData <- classificationData\n" +
                    "treeRoc <- roc(30,rocData,\"lod\",paste(\"Tree\",N_TREE,sep=\"\"))\n" +
                    "vqsrRoc <- roc(30,rocData,\"vqslod\",\"GMM\")\n" +
                    "rocs <- rbind(treeRoc,vqsrRoc)\n" +
                    "plotRocs(rocs,\"Filtering\")\n" +
                    "write(classLod,file=outputFile,ncolumns=1)\n" +
                    "dev.off()");
            script.close();

        } catch (FileNotFoundException e) {
            throw new StingException("File not found",e);
        }

        return new File(rOutput);
    }

    private File writeOld(List<VariantDatum> data, String oFileStr, String rOutput, String scriptFileString, String num_tree) {
        File oFile = new File(oFileStr);
        try {
            PrintStream out = new PrintStream(oFile);
            List<String> header = new ArrayList<String>(annotations.size()+HEADER_END.length+HEADER_START.length);
            header.addAll(Arrays.asList(HEADER_START));
            header.addAll(annotations);
            header.addAll(Arrays.asList(HEADER_END));
            out.printf("%s%n", Utils.join("\t",header));
            for ( VariantDatum d : data ) {
                StringBuilder bd = new StringBuilder();
                bd.append(d.contig);
                bd.append("\t");
                bd.append(d.start);
                bd.append("\t");
                bd.append(d.isSNP ? "SNP" : "INDEL");
                bd.append("\t");
                double[] annot = useFinal ? d.finalAnnotations : d.initialAnnotations;
                for ( double dd : annot ) {
                    bd.append(dd);
                    bd.append("\t");
                }
                bd.append(d.atPositiveTrainingSite ? "T" : ( d.atNegativeTrainingSite ? "F" : "U" ) );
                bd.append("\t");
                bd.append(d.isKnown ? "K" : "N");
                bd.append("\t");
                bd.append( d.atTruthSite ? "T" : ( d.atMonomorphicSite ? "F" : "U"));
                bd.append("\t");
                bd.append(d.getLod());
                out.printf("%s%n",bd.toString());
            }
        } catch (FileNotFoundException e) {
            throw new StingException("File not found",e);
        }

        // write the script
        try {
            PrintStream script = new PrintStream(new File(scriptFileString));
            script.printf("%s%n","require(\"ggplot2\")\n" +
                          "require(\"randomForest\")");
            script.printf("%s <- \"%s\"%n","inputFile",oFileStr);
            script.printf("%s <- \"%s\"%n","outputFile",rOutput);
            script.print("MAX_LOD = 6.04\n" +
                    "\n" +
                    "variantDataTable <- read.table(inputFile,header=TRUE)\n" +
                    "colNames <- colnames(variantDataTable)\n" +
                    "trainAnnot <- colNames[3:(length(colNames)-2)] # (chr,pos,type,alleles,{ATTRIBUTES},training,known)\n" +
                    "# get only the labeled points; drop the known rod attribute\n" +
                    "trainingData <- subset(variantDataTable,TRAINING!=\"U\")\n" +
                    "# remove the empty \"U\" level\n" +
                    "response <- factor(trainingData[[\"TRAINING\"]])\n" +
                    "# drop out training and known columns\n" +
                    "trainingData <- na.roughfix(trainingData[,c(trainAnnot)])\n" +
                    "# now grab all the data\n" +
                    "classificationData <- na.roughfix(variantDataTable[,c(trainAnnot)])\n" +
                    "# train a random forest on this training data\n" +
                    "trainedRF <- randomForest(x = trainingData,\n" +
                    "                          y = response,\n" +
                    "                          importance = TRUE,\n" +
                    "                          proximity = F,\n" +
                    "                          ntree = "+num_tree+",\n" +
                    "                          keep.forest = FALSE,\n" +
                    "                          xtest = classificationData,\n" +
                    "                          norm.votes = TRUE,\n" +
                    "                          do.trace=100)\n" +
                    "classProbs <- trainedRF$test$votes\n" +
                    "d1 <- (classProbs[,c(\"T\")]>(1.0-10^(-MAX_LOD)))\n" +
                    "d2 <- (classProbs[,c(\"T\")]<10^(-MAX_LOD)) \n" +
                    "classProbs[d1,c(\"T\")] <- 1.0-10^(-MAX_LOD)\n" +
                    "classProbs[d1,c(\"F\")] <- 10^(-MAX_LOD)\n" +
                    "classProbs[d2,c(\"T\")] <- 10^(-MAX_LOD)\n" +
                    "classProbs[d2,c(\"F\")] <- 1.0-10^(-MAX_LOD)\n" +
                    "classLod <- log10(classProbs[,c(\"T\")]) - log10(classProbs[,c(\"F\")])\n" +
                    "write(classLod,file=outputFile,ncolumns=1)");
            script.close();

        } catch (FileNotFoundException e) {
            throw new StingException("File not found",e);
        }

        return new File(rOutput);
    }

    public double evaluateDatum(VariantDatum d) {
        // note: the lod gets set to its own value.
        return d.getLod();
    }
}
