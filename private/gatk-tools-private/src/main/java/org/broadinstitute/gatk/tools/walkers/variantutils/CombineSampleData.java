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

package org.broadinstitute.gatk.tools.walkers.variantutils;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.arguments.GenotypeCalculationArgumentCollection;
import org.broadinstitute.gatk.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.gatk.tools.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.gatk.tools.walkers.genotyper.UnifiedGenotypingEngine;
import org.broadinstitute.gatk.tools.walkers.genotyper.afcalc.GeneralPloidyFailOverAFCalculatorProvider;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.genotyper.IndexedSampleList;
import org.broadinstitute.gatk.utils.genotyper.SampleList;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;

import java.io.File;
import java.util.*;

/**
 * Combines sequencing data from whole-genome and exome studies of the same sample
 *
 *  <p>
 * CombineSampleData merges genotyped VCF records that were produced as part of the reference model-based variant discovery pipeline
 * (see documentation for more details) using the '-ERC GVCF' or '-ERC BP_RESOLUTION' mode of the HaplotypeCaller and
 * --uniquifySamples in GenotypeGVCFs.  This tool combines uniqufied samples with the same base name and combines the PLs
 * in a mathematically sound manner.
 *
 *
 * <h3>Input</h3>
 * <p>
 * A VCF containing pairs of samples, as uniquified by GenotypeGVCFs. The set of sample calls in a pair should be derived from WGS and WEx data for the same sample.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A combined VCF with combined calls for each pair of samples specified and de-uniquified sample names.
 * </p>
 *
 * <h3>Examples</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T CombineSampleData \
 *   --variant vcf1.vcf \
 *   -o output.vcf
 * </pre>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T CombineSampleData \
 *   --variant vcf1.vcf \
 *   --uniquified_sample_name NA12878.variant \
 *   --uniquified_sample_name NA12878.variant2
 *   -o output.vcf
 * </pre>
 *
 * <h3> Caveats </h3>
 *
 * This tool assumes a ploidy of 2.
 * When combining genotypes, only PLs, ADs, and DPs will be copied. Other attributes like phasing or genotype filters
 * should be reapplied after using this tool.
 *
 */

/*TODO: when this tool is moved into protected the following will have to be addressed:
    * Do more robust error checking on sample name de-uniquification -- right now checks for pairs of <sampleName>.variantX and <sampleName>.variantY but should be extended to allow tagged VCF input into GenotypeGVCFs, which will produce names like <sampleName>.RODtagName
    * Move sample name uniqufication/de-uniquification to SampleListUtils.java
    * Check to make sure all genotype attributes are preserved after merge, e.g. allele phasing and genotype filters
    * Generalize for all ploidies?
    * Change GenotypeGVCFs --uniquifySamples argument from hidden (maybe still keep @advanced?)
    *
 */

@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARMANIP, extraDocs = {CommandLineGATK.class} )
public class CombineSampleData extends RodWalker<Integer,Integer> {

    /**
     * The gVCF files to merge together
     */
    @Input(fullName="variant", shortName = "V", doc="One or more input gVCF files", required=true)
    public List<RodBindingCollection<VariantContext>> variantCollections;
    final private List<RodBinding<VariantContext>> variants = new ArrayList<>();

    /**
     * Optionally, specify the names of samples to merge.  By default, merges will be done by de-uniquifying samples.
     */
    @Input(fullName="uniquified_sample_name", shortName="usn", doc="Uniquified names of samples to merge (should come in pairs) -- unlisted samples will be omitted", required = false)
    protected Set<String> selectUniquifiedSampleNameSet = new HashSet<String>();


    @Output(fullName="out", shortName = "o", doc="File to which variants should be written")
    protected VariantContextWriter vcfWriter = null;

    @ArgumentCollection
    public GenotypeCalculationArgumentCollection genotypeArgs = new GenotypeCalculationArgumentCollection();


    // the genotyping engine
    private UnifiedGenotypingEngine genotypingEngine;

    // the annotation engine
    private VariantAnnotatorEngine annotationEngine;


    private List<String> sampleNameList;


    public void initialize() {
        final List<String> annotationsToUse = new ArrayList<>();

        // collect the actual rod bindings into a list for use later
        for ( final RodBindingCollection<VariantContext> variantCollection : variantCollections )  {
            variants.addAll(variantCollection.getRodBindings());
            //annotationsToUse.addAll(variantCollection.getRodBindings().)
        }

        final GenomeAnalysisEngine toolkit = getToolkit();
        final Map<String, VCFHeader> vcfRods = GATKVCFUtils.getVCFHeadersFromRods(toolkit, variants);

        final SampleList samples = new IndexedSampleList(SampleUtils.getSampleList(vcfRods, GATKVariantContextUtils.GenotypeMergeType.UNIQUIFY));

        // create the genotyping engine
        genotypingEngine = new UnifiedGenotypingEngine(createUAC(), samples, toolkit.getGenomeLocParser(), GeneralPloidyFailOverAFCalculatorProvider.createThreadSafeProvider(toolkit, genotypeArgs, logger),
                toolkit.getArguments().BAQMode);



        // create the annotation engine
        //annotationEngine = new VariantAnnotatorEngine(Arrays.asList("none"), annotationsToUse, Collections.<String>emptyList(), this, toolkit);


        // Initialize VCF header
        final Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), true);
        headerLines.add(new VCFHeaderLine("source", "CombineSampleData"));

        if (selectUniquifiedSampleNameSet.isEmpty())
            selectUniquifiedSampleNameSet = SampleUtils.getSampleList(vcfRods, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

        //check for proper uniquified sample name format
        boolean improperSampleNameFormat = false;
        for (final String sampleName : selectUniquifiedSampleNameSet) {
            if (!sampleName.contains("variant")) improperSampleNameFormat = true;
        }
        if (improperSampleNameFormat)
            throw new UserException.BadInput("Uniqufied sample names are improperly formatted -- should be <sampleName>.variantX. Tagged input variants may have been used in GenotypeGVCFs.");

        sampleNameList = deuniquify(selectUniquifiedSampleNameSet);
        final VCFHeader vcfHeader = new VCFHeader(headerLines, sampleNameList);
        vcfWriter.writeHeader(vcfHeader);
    }

    public Integer reduceInit() { return 0; }

    public Integer map(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context) {
        if ( tracker == null || context == null || ref == null ) {
            return 0;
        }

        final Collection<VariantContext> vcs = tracker.getValues(variants, ref.getLocus());

        if (vcs.isEmpty())
            return 0;

        final VariantContextBuilder builder = new VariantContextBuilder(vcs.iterator().next());
        ArrayList<Genotype> storeGenotypes = new ArrayList<>();
        ArrayList<Genotype> mergedGenotypes = new ArrayList<>();


        for (final VariantContext vc : vcs) {
            String baseSampleName = sampleNameList.get(0).split(".variant")[0];
            //sample name is alphabetized, so uniquified samples will occur together
            for (final String sampleName : selectUniquifiedSampleNameSet) {
                final String currentSampleBase = sampleName.split(".variant")[0];

                //group together uniquified samples based on base name
                if(currentSampleBase.equals(baseSampleName)) {
                    storeGenotypes.add(vc.getGenotype(sampleName));
                }
                //after run of sequential same samples ends, process
                else {
                    //there should always be at least one sample in storeGenotypes, but we'll be rigorous
                    if (!storeGenotypes.isEmpty())
                        mergedGenotypes.add(combineSamplePLs(storeGenotypes, vc, baseSampleName));
                    storeGenotypes.clear();
                    baseSampleName = sampleName.split(".variant")[0];
                    storeGenotypes.add(vc.getGenotype(sampleName));
                }

            }
            //do the last sample
            if (!storeGenotypes.isEmpty())
                mergedGenotypes.add(combineSamplePLs(storeGenotypes, vc, baseSampleName));
        }
        builder.genotypes(mergedGenotypes);

        final VariantContext mergedVC = builder.make();

        final VariantContext regenotypedVC;
        if (selectUniquifiedSampleNameSet.isEmpty()) {
            //despite the name, calculate genotypes doesn't change PLs, but will recalculate qual
            regenotypedVC = genotypingEngine.calculateGenotypes(mergedVC);
            //calculateGenotypes can return null
            if (regenotypedVC == null)
                return 0;
        }
        else {
            //don't calculate new QUAL if we're subsetting samples because we can lose calls
            regenotypedVC = mergedVC;
        }

        final VariantContextBuilder builder2 = new VariantContextBuilder(regenotypedVC);
        //calculateGenotypes removes pretty much all of the attributes from mergedVC, so copy them back
        builder2.filters(mergedVC.getFilters());
        builder2.attributes(mergedVC.getAttributes());
        builder2.id(mergedVC.getID());
        VariantContextUtils.calculateChromosomeCounts(builder2, true);
        vcfWriter.add(builder2.make());

        return 1;
    }

    public Integer reduce(final Integer l, final Integer r) { return r + l; }

    /**
     * Creates a UnifiedArgumentCollection with appropriate values filled in from the arguments in this walker
     * @return a complete UnifiedArgumentCollection
     */
    private UnifiedArgumentCollection createUAC() {
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.genotypeArgs = genotypeArgs.clone();
        return uac;
    }

    private List<String> deuniquify(final Set<String> uniqueSampleNameSet) {
        final List<String> deuniquified = new ArrayList<>();
        for(final String sampleName : uniqueSampleNameSet) {
            final String baseName = sampleName.split(".variant")[0];
            if (!deuniquified.contains(baseName))
                deuniquified.add(baseName);
        }
        return deuniquified;
    }

    @Requires("!storeGenotypes.isEmpty()")
    private Genotype combineSamplePLs(final ArrayList<Genotype> storeGenotypes, final VariantContext vc, final String sampleName){
        //When combining WGS and WEx data, WEx may be no-call. In an attempt to preserve any extra genotype attributes,
        // build the combined genotype from whichever has the most genotype data or is "fullest"
        Genotype fullestGenotype = null;

        for (final Genotype g : storeGenotypes) {
            if (g.hasPL()) {
                fullestGenotype = g;
                break;
            }
        }
        //if neither genotype has PLs, just use the first one
        //precondition guarantees storeGenotypes isn't empty
        if (fullestGenotype == null) fullestGenotype = storeGenotypes.get(0);

        final GenotypeBuilder gBuilder = new GenotypeBuilder(fullestGenotype);
        final int N = vc.getNAlleles();
        //assume ploidy = 2 to generalize to multiallelic PLs
        final int[] PLsum = new int[(N*N-N)/2+N];
        final int[] ADsum = new int[N];
        int DPsum = 0;

        for (final Genotype g : storeGenotypes) {
            //add up PLs for data from both samples
            int[] currentPL;
            if (g.hasPL()) {
                currentPL = g.getPL();
                for (int i=0; i< PLsum.length; i++) {
                    PLsum[i] += currentPL[i];
                }
            }
            //add up AD for data from both samples
            int[] currentAD;
            if (g.hasAD()) {
                currentAD = g.getAD();
                for (int i=0; i< ADsum.length; i++) {
                    ADsum[i] += currentAD[i];
                }
            }
            if (g.hasDP())
                DPsum += g.getDP();
        }
        gBuilder.PL(PLsum);
        gBuilder.AD(ADsum);
        gBuilder.DP(DPsum);

        gBuilder.name(sampleName);
        GATKVariantContextUtils.updateGenotypeAfterSubsetting(fullestGenotype.getAlleles(), gBuilder,
                GATKVariantContextUtils.GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, GenotypeLikelihoods.fromPLs(PLsum).getAsVector(), vc.getAlleles());

        return gBuilder.make();
    }

    }
