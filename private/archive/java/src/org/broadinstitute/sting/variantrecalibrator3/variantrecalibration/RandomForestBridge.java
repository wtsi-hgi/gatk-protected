/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

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
        File fileForROutput = write(data,args.RInputFile, intermediateFile, args.RScriptFile, args.numTree,args.pdf);
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


    private File write(List<VariantDatum> data, String oFileStr, String rOutput, String scriptFileString, String num_tree, File pdf) {
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
                bd.append(d.loc.getContig());
                bd.append("\t");
                bd.append(d.loc.getStart());
                bd.append("\t");
                bd.append(d.isSNP ? "0" : "1");
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
                    "library(reshape)\n" +
                    "\n" +
                    //"pdf(\"" + pdf.getAbsolutePath() + "\",8,10)\n" +
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
                    "trainAnnot <- colNames[2:(length(colNames)-4)] # (chr,pos,type,{ATTRIBUTES},training,known,status,vqslod)\n" +
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
                    "# take a little bit of extra time to generate importance plots\n" +
                    "\n" +
                    "importancePlot <- function(rfs, title) {\n" +
                    "    # Combined importance plot\n" +
                    "    df = data.frame()\n" +
                    "    df <- melt(rfs$importance)\n" +
                    "    print(df)\n" +
                    "    names(df) <- c(\"RF.variable\", \"Importance.measure\", \"value\")\n" +
                    "    p <- ggplot(data=df, aes(y=reorder(RF.variable, value), x=value))\n" +
                    "    p <- p + facet_grid(. ~ Importance.measure, scale=\"free\")\n" +
                    "    p <- p + geom_point(size=3) #  + geom_linerange(aes(ymin=0, ymax=value), size=2)\n" +
                    "    p <- p + xlab(\"Score\") + ylab(\"Random forest training variable\")\n" +
                    "    print(p)\n" +
                    "}\n" +
                    "smallRF <- randomForest(x = trainingData,\n" +
                    "                          y = response,\n" +
                    "                          importance = TRUE,\n" +
                    "                          proximity = F,\n" +
                    "                          ntree = 120,\n" +
                    "                          keep.forest = TRUE,\n" +
                    "                          xtest = na.roughfix(variantDataTable[,c(trainAnnot)]),\n" +
                    "                          norm.votes=TRUE,\n" +
                    "                          do.trace=100)\n" +
                    "importancePlot(smallRF,\"ImportancePlot\")\n" +
                    "dev.off()");
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
