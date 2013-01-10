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

package org.broadinstitute.sting.gatk.walkers.misc;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
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
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypeType;
import org.broadinstitute.variant.variantcontext.VariantContext;

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
public class ScoreSeq extends RodWalker<Integer,Integer> {

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
                double gDosage = 2.0*Math.pow(10, -0.1 * g.getLikelihoods().getAsMap(false).get(GenotypeType.HOM_VAR));
                gDosage += Math.pow(10,-0.1*g.getLikelihoods().getAsMap(false).get(GenotypeType.HET));
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
