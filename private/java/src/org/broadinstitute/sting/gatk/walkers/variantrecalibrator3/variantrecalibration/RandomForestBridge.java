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
    private static final String[] HEADER_END = {"TRAINING","KNOWN"};

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
                bd.append(d.atTruthSite ? "K" : "N");
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
