package org.broadinstitute.sting.gatk.walkers.misc;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.R.RScriptExecutor;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.codecs.table.TableFeature;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 5/1/12
 * Time: 1:35 PM
 * To change this template use File | Settings | File Templates.
 */
public class ScoreSeqWalker extends RodWalker<Integer,Integer> {

    @Input(doc="A VCF containing the SNPs you want to test",required=true,shortName="V",fullName="variants")
    RodBinding<VariantContext> variantContext;

    @Input(doc="A tabular ROD (table or bedTable) listing, for each variant, the weight(s) associated with the variant",
            shortName="weightTable",fullName="weightTable",required=true)
    public RodBinding<TableFeature> weightTable = null;

    @Argument(doc="Delta value by which to modify weights. Defaults to 32/sqrt(nSamples)",required=false,shortName="delta",fullName="delta")
    public double delta = -1;

    @Argument(doc="Use soft-calls for dosages?",shortName="softCall",fullName="softCall",required=false)
    public boolean softCall = false;

    @Argument(doc="Re-weigh the dosage vector at no-calls, rather than rough-fixing the dosage",shortName="noroughfix",fullName="noRoughFix",required=false)
    public boolean noRoughFix = false;

    @Input(doc="An input file of sample covariates, tab-delimited. Random effects currently not supported, so please \n"
            +"encode cohorts or groups as dummy {0,1} variables",required=true,shortName="covFile",fullName="covariatesFile")
    public File covariateFile = null;

    @Argument(doc="The covariates to use in the model (column names of covariate file)",required=true,shortName="cov",fullName="covariates")
    public List<String> covariates = new ArrayList<String>();

    @Argument(doc="The name/header of the response variable to use in the model",required=true,shortName="resp",fullName="response")
    public String respName = null;

    @Argument(doc="Is the response vector binary (alters the link function)",required=false,shortName="b",fullName="binaryResponse")
    public boolean isBinaryResponse = false;

    @Argument(doc="The Sample ID (header of the ID column), defaults to \"IID\"",required=false,fullName="sampleIDKey")
    public String sampleIDKey = "IID";

    @Argument(doc="The header of the ID column for the variant weights, defaults to \"weight\"",required=false,shortName="wk",fullName="weightKey")
    public String weightKey = "weight";

    @Output(doc="The name of the covariate-and-dosage file to write",required=true,fullName="phenoOut",shortName="PO")
    public File outPheno;

    //@ArgumentCollection
    //ScoreSeqArgumentCollection ssCollection;

    @Output
    File out;

    private Map<String,Double> dosages;

    public void initialize() {
        // get number of samples in the VCF
        String vname = variantContext.getName();
        Set<String> sampleNames = SampleUtils.getSampleListWithVCFHeader(getToolkit(),Arrays.asList(vname));
        // sanity check the bound samples
        try {
            int numSamples = (new XReadLines(covariateFile)).readLines().size()-1;
            if ( numSamples < sampleNames.size() ) {
                throw new UserException("More samples in the VCF than provided in the covariates file. Either covariate file is empty, or many samples are missing");
            }
        } catch (FileNotFoundException e) {
            throw new StingException("Covariate file not found: "+covariateFile.getAbsolutePath(),e);
        }
        if ( delta < 0) {
            delta = 32.0/sampleNames.size();
        }

        dosages = new HashMap<String,Double>(sampleNames.size());
        for ( String s : sampleNames ) {
            dosages.put(s,0.0);
        }

    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer a, Integer b) {
        return a + b;
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null || ! tracker.hasValues(variantContext) )
            return 0;

        VariantContext vc = tracker.getFirstValue(variantContext);

        if ( vc.isFiltered() || ! vc.isBiallelic() )
            return 0;

        if ( ! tracker.hasValues(weightTable) || tracker.getValues(weightTable,ref.getLocus()).size() == 0 )
            throw new UserException("Variant Context at "+ref.getLocus().toString()+" does not have a weight associated with it");

        TableFeature weightFeature = tracker.getValues(weightTable,ref.getLocus()).get(0);

        incorporate(weightFeature,vc);

        return 1;
    }

    protected int incorporate(TableFeature weightFeature, VariantContext vc) {
        double weight = Double.parseDouble(weightFeature.get(weightKey));
        double modWeight = weight;
        if ( weight != 0.0 )
            modWeight = weight < 0 ? weight - delta : weight + delta;

        if ( softCall ) {
            for (Genotype g : vc.getGenotypes()) {
                double gDosage = 2.0*Math.pow(10, -0.1 * g.getLikelihoods().getAsMap(false).get(Genotype.Type.HOM_VAR));
                gDosage += Math.pow(10,-0.1*g.getLikelihoods().getAsMap(false).get(Genotype.Type.HET));
                dosages.put(g.getSampleName(),dosages.get(g.getSampleName())+gDosage*modWeight);
            }
        } else {
            for ( Genotype g : vc.getGenotypes() ) {
                if ( g.isHomRef() )
                    continue;
                double dosage = dosages.get(g.getSampleName());
                if ( g.isHet() )
                    dosage += modWeight;
                else if ( g.isHomVar() )
                    dosage += 2*modWeight;
                else // nocall
                    if ( noRoughFix )
                        continue;
                    else
                        dosage += modWeight*Double.parseDouble(vc.getAttribute("AF",0.0).toString());
                dosages.put(g.getSampleName(),dosage);
            }
        }

        return 1;
    }

    public void onTraversalDone(Integer records) {
        StringBuilder headerBuilder = new StringBuilder();
        headerBuilder.append(sampleIDKey);
        headerBuilder.append("\t"); // for the ease of reading table into R as a data frame
        headerBuilder.append(respName);
        headerBuilder.append("\t");
        for ( String s : covariates ) {
            headerBuilder.append(s);
            headerBuilder.append("\t");
        }
        headerBuilder.append("DOSAGE");
        XReadLines covariateFileLines;
        try {
            covariateFileLines = new XReadLines(covariateFile);
        } catch ( FileNotFoundException e) {
            throw new StingException("Covariate file was found earlier but now is not found. Did the filesystem monster consume it in its vengeful wrath?",e);
        }
        PrintStream phenoStream;
        try {
            phenoStream = new PrintStream(outPheno);
        } catch ( FileNotFoundException e ) {
            throw new UserException("Error opening output file for writing: "+outPheno.getAbsolutePath());
        }

        phenoStream.printf("%s%n",headerBuilder.toString());
        // get the indeces
        String[] covFileHeader = covariateFileLines.next().split("\t");
        Map<String,Integer> columnOffsets = new HashMap<String,Integer>(covFileHeader.length);
        for ( int i = 0; i < covFileHeader.length ; i++ ) {
            columnOffsets.put(covFileHeader[i],i);
        }

        for ( String samLine : covariateFileLines ) {
            String[] samData = samLine.split("\t");
            StringBuilder samBuilder = new StringBuilder();
            // sample name
            String samID = samData[columnOffsets.get(sampleIDKey)];
            samBuilder.append(samID);
            samBuilder.append("\t");
            samBuilder.append(samData[columnOffsets.get(respName)]);
            for ( String cov : covariates ) {
                samBuilder.append("\t");
                samBuilder.append(samData[columnOffsets.get(cov)]);
            }
            samBuilder.append("\t");
            samBuilder.append(String.format("%f",dosages.get(samID)));
            phenoStream.printf("%s%n",samBuilder.toString());
        }

        phenoStream.close();
        runAssociation();
    }

    protected void runAssociation() {
        File tmpRScript;
        PrintStream rScriptStream;
        try {
            tmpRScript = File.createTempFile("ScoreSeqAssoc",".R");
            rScriptStream = new PrintStream(tmpRScript);
        } catch (IOException e) {
            throw new UserException("Unable to open temporary file for writing R script. Make sure your temporary directory is writable.",e);
        }

        StringBuilder nullModelEq = new StringBuilder();
        nullModelEq.append(respName);
        nullModelEq.append(" ~ ");
        nullModelEq.append(Utils.join(" + ",covariates));
        nullModelEq.append(" + 1");
        StringBuilder altModelEq = new StringBuilder(nullModelEq.toString());
        altModelEq.append(" + ");
        altModelEq.append("DOSAGE");
        String family = isBinaryResponse ? "binomial" : "gaussian";

        rScriptStream.printf("%s", "library(\"nlme\")\n" +
                "library(\"MASS\")\n" +
                "covAndDosage <- read.table(\""+outPheno.getAbsolutePath()+"\",header=T)\n" +
                "covOnlyPQL <- glm("+nullModelEq.toString()+", family="+family+", data=covAndDosage)\n" +
                "withDosagePQL <- glm("+altModelEq.toString()+", family="+family+",data=covAndDosage)\n" +
                "# compute p-values from likelihood ratio test\n" +
                "lrt <- function(o1,o2) {\n" +
                "\tl0 <- logLik(o1)\n" +
                "\tl1 <- logLik(o2)\n" +
                "\tdif = as.vector(-2*(l0-l1))\n" +
                "\tdf <- attr(l1,\"df\")-attr(l0,\"df\")\n" +
                "\tpval <- pchisq(dif,df,lower.tail=F)\n" +
                "\tpval\n" +
                "}\n" +
                "pChiSqLogLik <- lrt(covOnlyPQL,withDosagePQL)\n" +
                "# calculate p-values via the scoreseq method\n" +
                "U <- sum(covOnlyPQL$residuals*covAndDosage$DOSAGE)\n" +
                "# vi = exp(gamma*Z)/((1+exp(gamma*Z))^2\n" +
                "r <- t(as.matrix(covOnlyPQL$coefficients))%*%t(covOnlyPQL$model)\n" +
                "VI <- exp(r)/(1+exp(r))^2\n" +
                "V_1 <- sum(VI*covAndDosage$DOSAGE*covAndDosage$DOSAGE)\n" +
                "V_2 <- colSums(t(VI*covAndDosage$DOSAGE)*covOnlyPQL$model)\n" +
                "ZZ <- as.matrix((as.matrix(covOnlyPQL$model)[1,]))\n" +
                "ZZ <- ZZ %*% t(ZZ)\n" +
                "ZZ <- VI[1] * ZZ\n" +
                "for ( j in 2:length(r) ) {\n" +
                "\tm <- as.matrix((as.matrix(covOnlyPQL$model)[1,]))\n" +
                "\tZZ <- ZZ + VI[j] * (m %*% t(m))\n" +
                "}\n" +
                "V <- V_1 - t(as.matrix(V_2))%*%ginv(ZZ)%*%as.matrix(V_2)\n" +
                "T <- U/sqrt(V)\n" +
                "pNorm = pnorm(T,lower.tail=F)\n" +
                "withDosagePQL$pChiSq <- pChiSqLogLik\n" +
                "withDosagePQL$U <- U\n" +
                "withDosagePQL$V <- V\n" +
                "withDosagePQL$Tval <- T\n" +
                "withDosagePQL$pNormOfTval <- pNorm\n" +
                "withDosagePQL$logPNormOfTval <- log(pNorm)\n" +
                "write(str(withDosagePQL),file = \""+out.getAbsolutePath()+"\")");
        rScriptStream.close();
        RScriptExecutor executor = new RScriptExecutor();
        executor.setExceptOnError(true);
        executor.addScript(tmpRScript);
        executor.exec();
    }

    protected void clearDosages() {
        // for subclasses that want to do multiple testing
        Map<String,Double> newDosage = new HashMap<String,Double>(dosages.size());
        for ( Map.Entry<String,Double> etry : dosages.entrySet() ) {
            newDosage.put(etry.getKey(),0.0);
        }
        dosages = newDosage;
    }

}
