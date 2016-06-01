/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").
* 
* WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
* WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
* NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
* 
* 1. DEFINITIONS
* 1.1 PROGRAM shall mean copyright in the object code and source code known as GATK3 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/gatk on the EFFECTIVE DATE.
* 
* 2. LICENSE
* 2.1 Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
* The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only. For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
* 2.2 No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD. LICENSEE shall ensure that all of its users agree to the terms of this Agreement. LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
* 2.3 License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
* 
* 3. PHONE-HOME FEATURE
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system ("PHONE-HOME") which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE’S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2016 Broad Institute, Inc.
* Notice of attribution: The GATK3 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
* LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
* 
* 5. INDEMNIFICATION
* LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
* 
* 6. NO REPRESENTATIONS OR WARRANTIES
* THE PROGRAM IS DELIVERED AS IS. BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
* IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
* 
* 7. ASSIGNMENT
* This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
* 
* 8. MISCELLANEOUS
* 8.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
* 8.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
* 8.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
* 8.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested. All notices under this Agreement shall be deemed effective upon receipt.
* 8.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
* 8.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
* 8.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.gatk.tools.walkers.rnaseq;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.tools.walkers.annotator.StrandBiasTableUtils;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Created by ami on 8/6/14.
 */
public class ASECaller extends RodWalker<Integer, Integer> {


    @Input(fullName="dnaVariant", shortName = "dna", doc="Input DNA VCF file", required=true)
    public RodBinding<VariantContext> dnaVariants;

    @Input(fullName="rnaVariant", shortName = "rna", doc="Input RNA VCF file", required=true)
    public RodBinding<VariantContext> rnaVariants;

    @Output(doc="File to which all variants with allele specific expression should be written")
    protected VariantContextWriter vcfWriter = null;

    @Argument(fullName = "sampleName", shortName = "sn", doc = "sample of interest", required = false)
    protected String sampleName = "NA12878";

    @Argument(fullName = "emitGQdifferances", shortName = "GQdiff", doc = "only emit sites where the different between the QC is higher then this threshold; default 0 ", required = false)
    protected int gqDiff = 0;

    @Argument(fullName = "emitQDdifferances", shortName = "GDdiff", doc = "only emit sites where QD is not more then X times the other QD this threshold; default 2  ", required = false)
    protected int dpDiff = 2;

    @Argument(fullName = "minGQ", shortName = "minGQ", doc = "only emit sites where both GQ is at least this threshold; default 20  ", required = false)
    protected int minGQ = 20;

    @Argument(fullName = "onlyHet", shortName = "het", doc = "only emit sites where the genotype is het; default true  ", required = false)
    protected boolean onlyHet = true;

    @Argument(fullName = "differentGenotypes", shortName = "diffGenotype", doc = "emit also sites where the genotypes are different; default true  ", required = false)
    protected boolean diffGenotypes = true;

    private static final double MIN_PVALUE = 1E-320;


    public void initialize() {
       //todo get a header
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        final String contig = context.getLocation().getContig();
        final long position = context.getPosition();

        List<RodBinding<VariantContext>> vatiants = new ArrayList<>();
        vatiants.add(dnaVariants);
        vatiants.add(rnaVariants);
        Collection<VariantContext> vcs = tracker.getValues(vatiants, context.getLocation());
        int numOfVCs = vcs.size();
        if (numOfVCs == 0) {
            logger.warn(contig+ ":"+position+" does not contain any VC");
            return 0;
        }
        if (numOfVCs > 2)
            logger.warn(contig+":"+position+" has more then 2 VC");
        if(numOfVCs == 1) {
            logger.warn(contig + ":" + position + " has only one VC");
            return 0;
        }
        VariantContext dnaVC = null, rnaVC= null;
        int countVC = 1;
        for(final VariantContext vc : vcs){
            if (numOfVCs == 2 && countVC == 1)
                dnaVC = vc;
            else if (numOfVCs == 2 && countVC == 2)
                rnaVC = vc;
            else{
                logger.warn(contig+":"+position+" strange number of VC: "+vc);
            }
            countVC++;
        }
        if(numOfVCs > 2)
            return 0;

        Genotype dnaGenotype, rnaGenotype;
        int dnaDP, rnaDP;
        int[] dnaADs, rnaADs;
        double dnaGQ, rnaGQ;
        if(dnaVC != null && dnaVC.hasGenotype(sampleName)) {
            dnaGenotype = dnaVC.getGenotype(sampleName);
            dnaDP = dnaGenotype.getDP();
            dnaGQ = dnaGenotype.getGQ();
            dnaADs = dnaGenotype.getAD();
        }
        else
            throw new UserException("no dna genotype or no matching sample name for "+sampleName);

        if(rnaVC != null && rnaVC.hasGenotype(sampleName)) {
            rnaGenotype = rnaVC.getGenotype(sampleName);
            rnaDP = rnaGenotype.getDP();
            rnaGQ = rnaGenotype.getGQ();
            rnaADs = rnaGenotype.getAD();
        }
        else
            throw new UserException("no rna genotype or no matching sample name for "+sampleName);

        if(dnaGenotype.sameGenotype(rnaGenotype) ){//&& dnaGQ >= minGQ && rnaGQ >= minGQ && (dnaDP/rnaDP) < dpDiff && (dnaDP/rnaDP) > 1/ dpDiff && (dnaGQ - rnaGQ > gqDiff || dnaGQ - rnaGQ < -gqDiff)) {
            if (!onlyHet || dnaGenotype.isHet()) {
                //if(dnaGQ >= lowDnaGQ){
                //    System.out.println("different genotype: possible ASE: " + contig + ":" + position + "\t" + dnaGenotype.toString() + "\t" + rnaGenotype.toString());
                //}

                int[][] table = new int[2][2];
                table[0][0] = dnaADs[0];
                table[0][1] = dnaADs[1];
                table[1][0] = rnaADs[0];
                table[1][1] = rnaADs[1];
                double pValue = StrandBiasTableUtils.FisherExactPValueForContingencyTable(table);
                final Object value = String.format("%.3f", QualityUtils.phredScaleErrorRate(Math.max(pValue, MIN_PVALUE)));
                double dnaRatio = (double)(dnaADs[0])/(dnaADs[1]+dnaADs[0]);
                double rnaRatio = (double)(rnaADs[0])/(rnaADs[1]+rnaADs[0]);
                System.out.println("same genotype: " + contig + ":" + position + "\t" + dnaGenotype.toString() + "\t" + rnaGenotype.toString() +"\t"+dnaRatio+"\t"+rnaRatio+"\t"+value);
            }
        }

        if(!dnaGenotype.sameGenotype(rnaGenotype) && diffGenotypes && dnaGQ >= minGQ && rnaGQ >= minGQ) {
            if(rnaGenotype.isHomVar() && dnaGenotype.isHet())
                System.out.println("different genotype: possible ASE: " + contig + ":" + position + "\t" + dnaGenotype.toString() + "\t" + rnaGenotype.toString());
            else if(dnaGenotype.isHomRef() && !rnaGenotype.isHomRef())
                System.out.println("different genotype: possible RNA editing: " + contig + ":" + position + "\t" + dnaGenotype.toString() + "\t" + rnaGenotype.toString());
            else
                System.out.println("different genotype: unknown case: " + contig + ":" + position + "\t" + dnaGenotype.toString() + "\t" + rnaGenotype.toString());
        }
        return 1;


    }

    @Override
    public Integer reduceInit() {
        return null;
    }

    @Override
    public Integer reduce(Integer value, Integer sum) {
        return null;
    }
}
