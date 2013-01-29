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

package org.broadinstitute.sting.gatk.walkers.simulatereads;

import cern.jet.random.Poisson;
import cern.jet.random.engine.MersenneTwister;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.gatk.walkers.annotator.ChromosomeCountConstants;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.variant.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.variant.vcf.*;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.text.TextFormattingUtils;

import java.io.PrintWriter;
import java.util.*;

/**
 * Generates simulated reads for variants
 */
@Reference(window=@Window(start=-200,stop=200))
public class SimulateReadsForVariants extends RefWalker<Integer, Integer> {
    @Argument(fullName = "vcf", shortName = "vcf", doc="Variants underlying the reads",required=true)
    protected VariantContextWriter variantsWriter;

    @Argument(fullName = "sites", shortName = "sites", doc="Variants sites",required=true)
    protected PrintWriter sitesWriter;

    @Output(fullName = "read", shortName = "reads", doc="Reads corresponding to variants",required=true)
    protected StingSAMFileWriter readWriter;

    @Argument(fullName="nSamples", shortName="NS", doc="Number of samples to simulate", required=false)
    public int nSamples = 1;

    @Argument(fullName="readDepth", shortName="DP", doc="Read depths to simulate", required=false)
    public List<Integer> readDepths = Arrays.asList(1);

    @Argument(fullName="errorRate", shortName="ER", doc="Phred-scaled error rate", required=false)
    public List<Integer> errorRates = Arrays.asList(20);

    @Argument(fullName="readLengths", shortName="RL", doc="Read length, in bp", required=false)
    public List<Integer> readLengths = Arrays.asList(3);

    public enum ReadSamplingMode { CONSTANT, POISSON };
    @Argument(fullName="readSamplingMode", shortName="RSM", doc="Sampling mode", required=false)
    public List<ReadSamplingMode> samplingModes = Arrays.asList(ReadSamplingMode.CONSTANT);

    @Argument(fullName="variantsPerBin", shortName="VPB", doc="No. of variants to generate for each bin", required=false)
    public int variantsPerBin = 1;

    @Argument(fullName="verbose", shortName="verbose", doc="Verbose", required=false)
    public boolean verbose = false;

    private class ParameterSet {
        int readDepth, readLength;
        ReadSamplingMode mode;
        byte[] readQuals;
        double errorRate; // in abs fraction (0.01 not Q20)
        int nVariants = 0;
        ParameterSet next = null;
        Poisson poissonRandom = null;
        Iterator<Integer> acs;
        int nSites = 0;

        public ParameterSet(int readDepth, int readLength, ReadSamplingMode mode, int phredErrorRate, ParameterSet next, List<Integer> ACs ) {
            this.readDepth = readDepth;
            this.readLength = readLength;
            this.mode = mode;
            this.readQuals = new byte[readLength];
            Arrays.fill(readQuals, (byte)phredErrorRate);
            this.errorRate = QualityUtils.qualToErrorProb((byte)phredErrorRate);
            this.next = next;
            nSites = ACs.size();
            acs = ACs.iterator();

            if ( mode == ReadSamplingMode.POISSON )
                poissonRandom = new Poisson(readDepth, new MersenneTwister((int)RANDOM_SEED));
        }

        public void incCount() { nVariants++; }
        public boolean done() { return ! acs.hasNext(); }
        public boolean hasNext() { return next != null; }

        public int combinations() {
            return nSites + ( hasNext() ? next.combinations() : 0);
        }
    }

    List<Integer> alleleCounts = new ArrayList<Integer>();

    ParameterSet parameters = null;
    SAMFileHeader header = null;

    private static String SAMPLE_PREFIX = "SAMPLE";
    public static final String PROGRAM_RECORD_NAME = "GATK SimulateReadsForVariants";

    List<String> sampleNames = new ArrayList<String>();
    Map<String, SAMReadGroupRecord> sample2RG = new HashMap<String, SAMReadGroupRecord>();

    private String sampleName(int i) { return sampleNames.get(i); }
    private SAMReadGroupRecord sampleRG(String name) { return sample2RG.get(name); }

    private static final long RANDOM_SEED = 1252863495;
    private static final Random ran = new Random(RANDOM_SEED);

    int SEPARATION_BETWEEN_SITES = 10;

    private SAMReadGroupRecord createRG(String name) {
        SAMReadGroupRecord rg = new SAMReadGroupRecord(name);
        rg.setPlatform("ILLUMINA");
        rg.setSample(name);
        return rg;
    }

	public void initialize() {
        // initialize sample I -> sample info map
        List<SAMReadGroupRecord> sampleRGs = new ArrayList<SAMReadGroupRecord>();

        for ( int i = 0; i < nSamples; i++ ) {
            sampleNames.add(String.format("%s%04d", SAMPLE_PREFIX, i));
            SAMReadGroupRecord rg = createRG(sampleName(i));
            sampleRGs.add(rg);
            sample2RG.put(sampleName(i), rg);
        }

        for ( int i = 0; i <= (2 * nSamples); i++) {
            int nCopies = (int)Math.round((2.0* nSamples) / (Math.max(i, 1)));
            for ( int j = 0; j < (nCopies * variantsPerBin); j++ )
                alleleCounts.add(i);
        }

        // initialize VCF headers
        // todo -- fill out header
        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
        headerLines.add(new VCFHeaderLine("source", "SimulateReadsForVariants"));
        VCFStandardHeaderLines.addStandardInfoLines(headerLines, true,
                VCFConstants.DOWNSAMPLED_KEY,
                VCFConstants.MLE_ALLELE_COUNT_KEY,
                VCFConstants.MLE_ALLELE_FREQUENCY_KEY);

        // FORMAT fields
        VCFStandardHeaderLines.addStandardFormatLines(headerLines, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_QUALITY_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_PL_KEY);
        headerLines.addAll(Arrays.asList(ChromosomeCountConstants.descriptions));
        headerLines.add(new VCFInfoHeaderLine("Q", 1, VCFHeaderLineType.Integer, "Q"));
        headerLines.add(new VCFInfoHeaderLine("MODE", 1, VCFHeaderLineType.Integer, "Mode"));

        variantsWriter.writeHeader(new VCFHeader(headerLines, new HashSet<String>(sampleNames)));

        // initialize BAM headers
        header = new SAMFileHeader();
        header.setSequenceDictionary(getToolkit().getReferenceDataSource().getReference().getSequenceDictionary());
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        header.setReadGroups(sampleRGs);

        final SAMProgramRecord programRecord = new SAMProgramRecord(PROGRAM_RECORD_NAME);
        final ResourceBundle headerInfo = TextFormattingUtils.loadResourceBundle("StingText");
        programRecord.setProgramVersion(headerInfo.getString("org.broadinstitute.sting.gatk.version"));
        programRecord.setCommandLine(getToolkit().createApproximateCommandLineArgumentString(getToolkit(), this));
        header.setProgramRecords(Arrays.asList(programRecord));

        readWriter.writeHeader(header);

        // set up feature sets
        for ( int readLength : readLengths ) {
            if ( readLength % 2 == 0 ) throw new UserException.BadArgumentValue("readLength", "Read lengths must be odd");

            for ( ReadSamplingMode mode : samplingModes ) {
                for ( int errorRate : errorRates ) {
                    for ( int readDepth : readDepths ) {
                        parameters = new ParameterSet(readDepth, readLength, mode, errorRate, parameters, alleleCounts);
                    }
                }
            }
        }
        logger.info("Total number of combinations " + parameters.combinations());
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( parameters.done() ) {
            if ( parameters.hasNext() )
                parameters = parameters.next;
            else
                return 0;   // early abort, we're done generating
        }

        if ( ref.getLocus().getStart() < parameters.readLength || ! BaseUtils.isRegularBase(ref.getBase()) )
            return 0;

        if ( ref.getLocus().getStart() % (parameters.readLength + SEPARATION_BETWEEN_SITES) != 0 )
            return 0;

        byte[] refBases = getBasesForReads(ref, parameters.readLength);

        // at this point, we want to generate variants and reads for the parameters in parameters
        int AC = parameters.acs.next();
        VariantContext vc = generateVariant(context.getLocation(), ref.getBase(), AC, parameters);
        if ( verbose ) logger.info(String.format("Generating reads for %s", vc));
        ReadBackedPileup rbp = generateRBPForVariant(context.getLocation(), vc, refBases, parameters);

        // BED is zero based
        sitesWriter.printf("%s %d %d%n", ref.getLocus().getContig(), ref.getLocus().getStart()-1, ref.getLocus().getStart() );
        variantsWriter.add(vc);
        for ( GATKSAMRecord read : rbp.getReads() ) readWriter.addAlignment(read);

        parameters.incCount();

        return 0;
    }

    private byte[] getBasesForReads(ReferenceContext ref, int readLength) {
        int center = (int)(ref.getLocus().getStart() - ref.getWindow().getStart());
        int start = center - ((readLength - 1) / 2);
        byte[] bases = new byte[readLength];
        System.arraycopy(ref.getBases(), start, bases, 0, readLength);
        return bases;
    }

    private VariantContext generateVariant( GenomeLoc loc, byte refBase, int AC, ParameterSet params ) {
        Allele ref = Allele.create(refBase, true);
        Allele alt = Allele.create(BaseUtils.baseIndexToSimpleBase(BaseUtils.getRandomBaseIndex(BaseUtils.simpleBaseToBaseIndex(refBase))));
        List<Allele> alleles = AC == 0 ? Arrays.asList(ref) : Arrays.asList(ref, alt);

        List<Allele> homRef = Arrays.asList(ref, ref);
        List<Allele> het    = Arrays.asList(ref, alt);
        List<Allele> homAlt = Arrays.asList(alt, alt);

        List<Genotype> genotypes = new ArrayList<Genotype>();
        double p = AC / (2.0 * nSamples);
        //double q = 1 - p;
        int nHomAlt = (int) Math.round(p * p * nSamples);
        int nHet    = AC - nHomAlt * 2;
        //int nHet    = (int) Math.round(2 * p * q * nSamples);
        for ( int i = 0; i < nSamples; i++ ) {
            List<Allele> genotype;

            if ( i < nHomAlt ) { genotype = homAlt; }
            else if ( i < (nHet + nHomAlt) ) { genotype = het; }
            else { genotype = homRef; }

            genotypes.add(GenotypeBuilder.create(sampleName(i), genotype));
        }

        Map<String, Object> attributes = new LinkedHashMap<String, Object>();
        attributes.put(VCFConstants.ALLELE_COUNT_KEY, AC);
    //    attributes.put(VCFConstants.SAMPLE_NUMBER_KEY, nSamples);
        attributes.put(VCFConstants.ALLELE_NUMBER_KEY, 2 * nSamples);
        attributes.put("Q", params.readQuals[0]);
        attributes.put("MODE", params.mode);
        attributes.put("DP", params.readDepth);
        
        return new VariantContextBuilder("anonymous", loc.getContig(), loc.getStart(), loc.getStart(), alleles)
                .genotypes(genotypes).passFilters().attributes(attributes).make();
    }

    private ReadBackedPileup generateRBPForVariant( GenomeLoc loc, VariantContext vc, byte[] refBases, ParameterSet params ) {
        List<GATKSAMRecord> reads = new ArrayList<GATKSAMRecord>();
        int offset = (params.readLength - 1) / 2;

        int start = (int)(loc.getStart() - (params.readLength - 1) / 2);
        byte altBase = vc.isVariant() ? vc.getAlternateAllele(0).getBases()[0] : 0;
        byte[] refHaplotype = Arrays.copyOf(refBases, refBases.length);
        byte[] altHaplotype = Arrays.copyOf(refBases, refBases.length);
        altHaplotype[(params.readLength - 1) / 2] = altBase;

        int gi = 0;
        for ( Genotype g : vc.getGenotypes() ) {
            int myDepth = sampleDepth(params);
            for ( int d = 0; d < myDepth; d++ ) {
                byte[] readBases = trueHaplotype(g, refHaplotype, altHaplotype);
                addMachineErrors(readBases, params.errorRate);

                GATKSAMRecord read = new GATKSAMRecord(header);
                read.setBaseQualities(params.readQuals);
                read.setReadBases(readBases);
                read.setReadName("FOO");
                read.setCigarString(params.readLength + "M");
                read.setReadPairedFlag(false);
                read.setAlignmentStart(start);
                read.setMappingQuality(60);
                read.setReferenceName(loc.getContig());
                read.setReadNegativeStrandFlag(gi++ % 2 == 0);
                read.setAttribute("RG", sampleRG(g.getSampleName()).getReadGroupId());
            
                reads.add(read);
            }
        }

        return new ReadBackedPileupImpl(loc, reads, offset);
    }

    private int sampleDepth(ParameterSet params) {
        switch ( params.mode ) {
            case CONSTANT: return params.readDepth;
            case POISSON: return params.poissonRandom.nextInt();
            default:
                throw new IllegalStateException("Unexpected DepthSamplingType " + params.mode);
        }
    }

    private byte[] trueHaplotype(Genotype g, byte[] refHaplotype, byte[] altHaplotype) {
        double refP = 0.0;

        if ( g.isHomRef() )     refP = 1;
        else if ( g.isHet() )   refP = 0.5;
        else                    refP = 0.0;

        return Arrays.copyOf(ran.nextDouble() < refP ? refHaplotype : altHaplotype, refHaplotype.length);
    }

    private void addMachineErrors(byte[] readBases, double errorRate) {
        for ( int i = 0; i < readBases.length; i++ ) {
            double r = ran.nextDouble();
            if ( r < errorRate ) {
                byte errorBase = BaseUtils.baseIndexToSimpleBase(BaseUtils.getRandomBaseIndex(BaseUtils.simpleBaseToBaseIndex(readBases[i])));
                if ( errorBase == readBases[i] ) throw new IllegalStateException("Read and error bases are the same");
                readBases[i] = errorBase;
            }
        }
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer counter, Integer sum) {
        return counter + sum;
    }

    public void onTraversalDone(Integer sum) {
        //variantsWriter.close();
        sitesWriter.close();
    }
}