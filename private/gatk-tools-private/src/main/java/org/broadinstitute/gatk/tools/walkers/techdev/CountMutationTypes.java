/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 (“BROAD”) and the LICENSEE and is effective at the date the downloading is completed (“EFFECTIVE DATE”).
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
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system (“PHONE-HOME”) which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE’S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.techdev;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.report.GATKReport;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFUtils;

import java.io.PrintStream;
import java.util.*;

/**
 * Simple walker to plot the percentages of each type of SNP.
 *
 * <p>
 *  Features of this walker:
 *  <li>allow filtering of sites based on JEXL expressions </li>
 *  <li>allow filtering based on variant Quals (default 20)</li>
 *  <li>produce a sorted list of the SNP type</li>
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * The VCF file and an optional interval list, JEXL expressions and variant qual threshold for filtering
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A GATK Report with the percentages of each type of SNP (sorted by the SNP type: A -> C ... T->G)
 *
 * <p/>
 * <h3>Examples</h3>
 * <pre>
 * java -Xmx4g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T CountMutationTypes \
 *   -V myData.vcf \
 *   -L interesting.intervals \
 *   -o report.grp
 * </pre>
 *
 * <pre>
 * java -Xmx4g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T CountMutationTypes \
 *   -V myData.vcf \
 *   -L interesting.intervals \
 *   -select "QD > 2"
 *   -minQual 500
 *   -o report.grp
 * </pre>
 *
 * <pre>
 * java -Xmx4g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T CountMutationTypes \
 *   -V myData.vcf \
 *   -L interesting.intervals \
 *   -select "QD > 2"
 *   --printPositions
 *   -o report.grp
 * </pre>
 *
 * @author ami
 * @since 4/9/14.
 */

@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARMANIP, extraDocs = {CommandLineGATK.class} )
public class CountMutationTypes extends RodWalker<Integer, Integer> {

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Output(fullName = "out", shortName = "o", doc = "output file, if not provided will print to stdout", required = false)
    protected PrintStream out = null;

    /**
     *  CountMutationTypes accepts any number of JEXL expressions (so you can have two filters by using
     *  --selectExpression "X < 1" --selectExpression "X > 2").
     */
    @Argument(fullName="selectExpression", shortName="select", doc="One or more expression used with INFO fields to filter", required=false)
    protected ArrayList<String> SELECT_EXPRESSIONS = new ArrayList<>();

    @Argument(fullName = "minQual",shortName = "minQual", doc = "min qual to include in the analysis", required = false)
    protected double minQual = 20;

    @Argument(fullName = "printPositions", shortName = "positions", doc = "print VCF files with the mutation sites that were counted (one file per mutation type, called <input file name>.printMutation.<mutation type, i.g. A_G>.vcf")
    protected boolean printPositions = false;

    // JEXL expressions for the filters
    private List<VariantContextUtils.JexlVCMatchExp> jexls = null;
    private final ArrayList<String> selectNames = new ArrayList<>();
    private int countMutations = 0;
    private Map<AllelePair,MutableInt> mutationCounter;
    private Map<AllelePair,List<VariantContext>> mutationVC;
    private VCFHeader vcfHeader;

    public void initialize() {

        for (int i = 0; i < SELECT_EXPRESSIONS.size(); i++) {
            selectNames.add(String.format("select-%d", i));
        }
        jexls = VariantContextUtils.initializeMatchExps(selectNames, SELECT_EXPRESSIONS);
        mutationCounter = new TreeMap<>(); //use TreeMap to get a sorted output
        mutationVC = new TreeMap<>(); //hold the VC in order to be able to print the sites if needed.
        if(printPositions){
            final Map<String, VCFHeader> vcfRods = GATKVCFUtils.getVCFHeadersFromRods(getToolkit());
            final Set<String> samples = SampleUtils.getSampleList(vcfRods, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);
            final Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), true);
            vcfHeader = new VCFHeader(headerLines, samples);
        }
    }

    /**
     *
     * @param tracker  the meta-data tracker
     * @param ref      the reference base
     * @param context  the context for the given locus
     * @return 1 if the locus was successfully processed, 0 if otherwise
     */
    @Override
    public Integer map(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context) {
        if ( tracker == null )
            return 0;

        final Collection<VariantContext> VCs = tracker.getValues(variantCollection.variants, context.getLocation());

        if ( VCs == null || VCs.size() == 0) {
            return 0;
        }

        for ( final VariantContext vc : VCs ) {
            if(countableVariant(vc)){
                addMutationToStat(vc);
            }
        }
        return 1;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(final Integer value, final Integer sum) {
        return sum + value;
    }

    /**
     * Tell the user the number of loci processed and close out the new variants file.
     * Create a GATK report with the percentages of each mutation type
     * (optional) print VCF files (one for each mutation type) with the records of the counted mutation
     *
     * @param result  the number of loci seen.
     */
    @Override
    public void onTraversalDone(final Integer result) {
        logger.info(result + " records processed.");
        final GATKReport report = GATKReport.newSimpleReport("MutationTypeDistribution","SNP","counts","percentage");

        for (final Map.Entry<AllelePair,MutableInt> entry : mutationCounter.entrySet()){
            report.addRow("<" + entry.getKey().first + "," + entry.getKey().second + ">", entry.getValue().getValue(), (double) 100 * entry.getValue().getValue() / countMutations);
        }
        report.print(out);

        if(printPositions){
            for (final Map.Entry<AllelePair,List<VariantContext>> entry: mutationVC.entrySet()){
                final String mutationType = entry.getKey().first.getBaseString() +"_"+entry.getKey().second.getBaseString();
                final VariantContextWriter vcfWriter = new VariantContextWriterBuilder()
                        .setOutputFile("mutationFile."+mutationType+".vcf")
                        .setReferenceDictionary(getToolkit().getMasterSequenceDictionary())
                        .build();
                vcfWriter.writeHeader(vcfHeader);
                for(final VariantContext vc : entry.getValue())
                    vcfWriter.add(vc);
                vcfWriter.close();

            }

        }
    }

    private boolean matchAllJexlExpressions (final VariantContext vc) {
        for ( final VariantContextUtils.JexlVCMatchExp jexl : jexls ) {
            if ( !VariantContextUtils.match(vc, jexl) )
                return false;
        }
        return true;
    }

    /**
    * a variant is countable if it is bi-allelic SNP that pass the user provided variant qual threshold and match the Jexl expressions
    **/
    private boolean countableVariant(final VariantContext vc){
        return vc.isSNP() && vc.isBiallelic() && vc.getPhredScaledQual() >= minQual && matchAllJexlExpressions(vc);
    }

    private void addMutationToStat(final VariantContext vc){
        countMutations++;
        final AllelePair pair = new AllelePair(vc.getReference(),vc.getAlternateAllele(0)); //we know it is a bi-allelic
        final MutableInt count = mutationCounter.get(pair);
        if(count == null){
            mutationCounter.put(pair,new MutableInt());
            if(printPositions){
                final List<VariantContext> vcList = new LinkedList<>();
                vcList.add(vc);
                mutationVC.put(pair,vcList);
            }
        }
        else{
            count.increment();
            if(printPositions)
                mutationVC.get(pair).add(vc);
        }
    }


    //--------------------- private classes -------------

    private final class MutableInt {
        private int value = 1; // we start at 1 since we're counting
        protected void increment() {++value;}
        protected int getValue() {return value;}
    }

    private final class AllelePair implements Comparable {

        Allele first;
        Allele second;

        public AllelePair(final Allele x, final Allele y){
            first = x;
            second = y;
        }

        public int compareTo(final Object other ){
            if(other == null)
                throw new UserException.BadInput("try to compare a AllelePair to null");
            final int compareFirst = first.compareTo(((AllelePair)other).first);
            return (compareFirst != 0) ?  compareFirst  : second.compareTo(((AllelePair)other).second);

        }
    }

}
