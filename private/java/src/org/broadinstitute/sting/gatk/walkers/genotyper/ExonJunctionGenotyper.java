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

package org.broadinstitute.sting.gatk.walkers.genotyper;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.DuplicateReadFilter;
import org.broadinstitute.sting.gatk.filters.FailsVendorQualityCheckFilter;
import org.broadinstitute.sting.gatk.filters.MappingQualityZeroFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportColumn;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.genotyper.afcalc.AFCalc;
import org.broadinstitute.sting.gatk.walkers.genotyper.afcalc.AFCalcFactory;
import org.broadinstitute.sting.gatk.walkers.genotyper.afcalc.AFCalcResult;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.clipping.ReadClipper;
import org.broadinstitute.sting.utils.codecs.refseq.RefSeqFeature;
import org.broadinstitute.sting.utils.codecs.table.TableFeature;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.variant.vcf.VCFConstants;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 11/16/11
 */
@ReadFilters({DuplicateReadFilter.class,FailsVendorQualityCheckFilter.class,MappingQualityZeroFilter.class})
public class ExonJunctionGenotyper extends ReadWalker<ExonJunctionGenotyper.EvaluationContext,ExonJunctionGenotyper.ECLikelihoods> {


    /**
     * A raw, unfiltered, highly sensitive callset in VCF format.
     */
    @Output(doc="File to which variants should be written", required = true)
    protected VariantContextWriter vcfWriterBase = null;

    protected VariantContextWriter vcfWriter;

    @Input(shortName="r",fullName="refSeq",required=true,doc="The RefSeq Gene definition track")
    public RodBinding<RefSeqFeature> refSeqRodBinding;

    @Input(shortName="H",fullName="insertHistogram",required=true,doc="The insert size histogram per read group, either flat file or GATK report formatted. One or more of these can be a .list file containing paths to multiple of such histogram files.")
    public List<File> readGroupInsertHistogram;

    @Input(shortName="Y",fullName="hypothesis",required=true,doc="The set of gene(transcript) hypothesis generated by the Hypothesis Generator that we want to score")
    public RodBinding<TableFeature> hypothesisRodBinding;

    @Argument(shortName="FSW", fullName="fullSmithWaterman",doc="Run SmithWaterman of every exon-intron overlapping read against exon-exon junction. Defaults to doing this only for unmapped pairs of mapped read.",required=false)
    public boolean fullSmithWaterman = false;

    @Hidden
    @Input(fullName="maxInsertSize",required=false)
    public int MAX_INSERT_SIZE = 1000;

    public final static byte MIN_TAIL_QUALITY = 6;

    private Map<String,byte[]> insertQualsByRG = new HashMap<String,byte[]>();

    private boolean initialized = false;
    private IndexedFastaSequenceFile referenceReader;
    private Set<String> samples;
    private List<Allele> noCallAlleles = Arrays.asList(new Allele[]{Allele.NO_CALL,Allele.NO_CALL});
    private Map<String,JunctionHypothesis> activeHypotheses;

    // the likelihoods engine
    IntronLossGenotypeLikelihoodCalculationModel ilglcm= new IntronLossGenotypeLikelihoodCalculationModel(logger,fullSmithWaterman);

    public ECLikelihoods reduceInit() {
        return new ECLikelihoods(ilglcm);
    }

    public ECLikelihoods reduce(EvaluationContext context, ECLikelihoods prevRed) {
        if ( context != null && context.read.getReadGroup() != null ) {
            logger.debug(String.format("%s   %s",prevRed.getLocation(), getToolkit().getGenomeLocParser().createGenomeLoc(context.read)));
            if ( prevRed.getLocation() != null && (context.read.getReferenceIndex() > prevRed.getLocation().getContigIndex() || context.read.getUnclippedStart() > prevRed.getLocation().getStop()) ) {
                printToFile(prevRed);
                prevRed.purge();
                activeHypotheses.clear();
            }
            prevRed.update(context);
        }

        return prevRed;
    }

    public EvaluationContext map(ReferenceContext context, GATKSAMRecord read, RefMetaDataTracker metaDataTracker) {
        if ( read.getReadGroup() == null ) { return null; }

        final HashSet<TableFeature> hypotheses = new HashSet<TableFeature>(metaDataTracker.getValues(hypothesisRodBinding));

        //logger.debug("Tracker Size "+Integer.toString(metaDataTracker.getAllCoveringRods().size())+" Hyp: "+hypotheses.size());

        if ( hypotheses.size() == 0 ) {
            return null;
        }

        //logger.debug(hypotheses);

        GenomeLoc readLoc = null;
        if ( ! read.getReadUnmappedFlag() ) {
            readLoc =  getToolkit().getGenomeLocParser().createGenomeLoc(read);
        } else if ( ! read.getMateUnmappedFlag() ) {
            readLoc = getToolkit().getGenomeLocParser().createGenomeLoc(read.getMateReferenceName(),read.getMateAlignmentStart(),read.getMateAlignmentStart()+read.getReadLength());
        }

        if (readLoc == null ) {
            return null;
        }

        final Map<String,RefSeqFeature> refSeqFeatures = new HashMap<String,RefSeqFeature>(16);
        for (final RefSeqFeature refSeqFeature : metaDataTracker.getValues(refSeqRodBinding) ) {
            refSeqFeatures.put(refSeqFeature.getTranscriptUniqueGeneName(),refSeqFeature);
        }

        EvaluationContext ec = new EvaluationContext();

        if ( ! read.getReadUnmappedFlag() ) {
            int start = read.getUnclippedStart();
            int stop= read.getUnclippedEnd();
            int leftHardClipAlready = read.getCigar().getCigarElement(0).getOperator().equals(CigarOperator.H) ? read.getCigar().getCigarElement(0).getLength() : 0;
            byte[] refBases = getToolkit().getReferenceDataSource().getReference().getSubsequenceAt(read.getReferenceName(),start, stop).getBases();
            GATKSAMRecord initRead = ReadClipper.hardClipLowQualEnds(read, (byte) 12);
            if ( initRead.getCigarString() == null ) {
                return null;
            }
            GATKSAMRecord clippedRead = ReadClipper.revertSoftClippedBases(initRead);

            if( clippedRead.getCigar().isEmpty()) {
                return null;
            }
            int offset = clippedRead.getCigar().getCigarElement(0).getOperator().equals(CigarOperator.H) ? clippedRead.getCigar().getCigarElement(0).getLength() - leftHardClipAlready : 0;
            clippedRead.setAttribute("NM", AlignmentUtils.getMismatchCount(clippedRead, refBases, offset).numMismatches);

            ec.read = clippedRead;
        } else {
            ec.read = ReadClipper.hardClipLowQualEnds(read, (byte) 12);
        }

        ec.refSeqByName = refSeqFeatures;
        // todo -- this is wildly inefficient, why regenerate hypotheses for every read?! just store in a map.
        ec.hypotheses = new TreeSet<JunctionHypothesis>();
        for ( TableFeature f : hypotheses ) {
            if ( activeHypotheses.containsKey(f.getValue(1))) {
                ec.hypotheses.add(activeHypotheses.get(f.getValue(1)));
            } else {
                RefSeqFeature fe = refSeqFeatures.get(f.getValue(1));
                if ( fe == null ) {
                    // cannot form a hypothesis without a ref seq feature
                    continue;
                }
                try {
                    for ( List<Pair<Integer,Integer>> junc : JunctionHypothesis.unwrapPairs(f.getValue(3))) {
                        if ( junc.size() > 0 ) {
                            JunctionHypothesis jc = new JunctionHypothesis(f.getValue(1),fe,junc,referenceReader);
                            ec.hypotheses.add(jc);
                            activeHypotheses.put(f.getValue(1),jc);
                        }
                    }
                } catch (IllegalArgumentException  e ) {
                    // just telling us that there were disparate events that could not be consolidated into a retrogene
                    continue;
                }
            }
        }

        return ec;
    }

    @Override
    public void initialize() {

        Set<String> sampleStr = SampleUtils.getSAMFileSamples(getToolkit());
        ilglcm.setSamples(sampleStr);

        vcfWriter = VariantContextWriterFactory.sortOnTheFly(vcfWriterBase,1500000);
        vcfWriter.writeHeader(new VCFHeader(new HashSet<VCFHeaderLine>(), sampleStr));

        try {
            for ( File readGroupInsertHistogramFile : readGroupInsertHistogram ) {
                logger.debug("Reading: "+readGroupInsertHistogramFile.getAbsolutePath());
                if ( readGroupInsertHistogramFile.getAbsolutePath().endsWith(".list")) {
                    // this is a list of histograms, read each one separately
                    for ( String line : new XReadLines(readGroupInsertHistogramFile) ) {
                        logger.debug("Reading "+line+" within "+readGroupInsertHistogramFile.getAbsolutePath());
                        readHistoFile(new File(line));
                    }
                } else {
                    readHistoFile(readGroupInsertHistogramFile);
                }
            }
        } catch (FileNotFoundException e) {
            throw new UserException("Histogram file not found "+e.getMessage(),e);
        } catch (IOException e) {
            throw new StingException("IO Exception",e);
        }

        StringBuffer debug = new StringBuffer();
        for ( Map.Entry<String,byte[]> etry : insertQualsByRG.entrySet() ) {
            debug.append("  ");
            debug.append(etry.getKey());
            debug.append("->");
            debug.append(Arrays.deepToString(ArrayUtils.toObject(etry.getValue())));
        }

        logger.debug(debug);

        ilglcm.initialize(refSeqRodBinding,getToolkit().getGenomeLocParser(),insertQualsByRG);

        try {
            // fasta reference reader to supplement the edges of the reference sequence
            referenceReader = new CachingIndexedFastaSequenceFile(getToolkit().getArguments().referenceFile);
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(getToolkit().getArguments().referenceFile,ex);
        }

        samples = new HashSet<String>(128);
        for (SAMFileHeader h : getToolkit().getSAMFileHeaders() ) {
            samples.addAll(SampleUtils.getSAMFileSamples(h));
        }

        activeHypotheses = new HashMap<String,JunctionHypothesis>();
    }

    private void readHistoFile(File readGroupInsertHistogramFile ) throws FileNotFoundException, IOException{ // exceptions handled upstream
        XReadLines xrl = new XReadLines(readGroupInsertHistogramFile);
        if ( ! xrl.hasNext() ) {
            logger.warn("The histogram file appears to be empty: "+readGroupInsertHistogramFile.getAbsolutePath());
            return;
        }
        if ( ! xrl.next().startsWith("##:") ) {
            xrl.close();
            for ( String entry : new XReadLines(readGroupInsertHistogramFile) ) {
                String[] split1 = entry.split("\\t");
                String id = split1[0];
                String[] histogram = split1[1].split(";");
                byte[] quals = new byte[histogram.length];
                int idx = 0;
                for ( String histEntry : histogram ) {
                    try {
                        quals[idx++] = Byte.parseByte(histEntry);
                    } catch( NumberFormatException e) {
                        quals[idx-1] = QualityUtils.trueProbToQual(Double.parseDouble(histEntry));
                    }
                }

                insertQualsByRG.put(id,quals);
            }
        } else {
            xrl.close();
            GATKReport report = new GATKReport(readGroupInsertHistogramFile);
            GATKReportTable reportTable = report.getTable("InsertSizeDistributionByReadGroup");
            // rows are insert sizes, columns are read groups
            for (GATKReportColumn reportColumn : reportTable.getColumnInfo() ) {
                // annoyingly, the column has no knowledge of its own rows
                int sum = 0;
                for ( int row = 0; row < reportTable.getNumRows(); row++ ) {
                    sum += Integer.parseInt( (String) reportTable.get(row,reportColumn.getColumnName()));
                }
                byte[] rgHist = new byte[MAX_INSERT_SIZE];
                for ( int row = 0; row < reportTable.getNumRows(); row++ ) {
                    final int insertSize = Integer.parseInt( (String) reportTable.get(row,0));
                    int val = 1;
                    if ( insertSize < rgHist.length ) {
                        val = Integer.parseInt( (String) reportTable.get(row,reportColumn.getColumnName()));
                    }
                    rgHist[row] = QualityUtils.errorProbToQual((((double) val) / sum), QualityUtils.MAX_SAM_QUAL_SCORE);
                }

                insertQualsByRG.put(reportColumn.getColumnName(),rgHist);
            }
        }
    }

    public void onTraversalDone(ECLikelihoods likelihoods) {
        printToFile(likelihoods);
        vcfWriter.close();
    }

    private void printToFile(ECLikelihoods likelihoods) {
        if ( likelihoods.context == null ) {
            return;
        }
        logger.info("HYPO_SIZE: "+Integer.toString(likelihoods.context.hypotheses.size()));
        for ( JunctionHypothesis hypothesis : likelihoods.context.hypotheses ) {
            Map<String,double[]> rawLiks = likelihoods.likelihoods.get(hypothesis);
            if ( rawLiks.size() == 0 ) {
                continue;
            }
            Allele alt = Allele.create(hypothesis.toString(),false);
            GenomeLoc refPos = hypothesis.getLocation().getStartLocation();
            //logger.info(hypothesis.getLocation());
            Allele ref = Allele.create(referenceReader.getSubsequenceAt(refPos.getContig(), refPos.getStart(), refPos.getStop()).getBases(), true);
            Byte paddingBase = referenceReader.getSubsequenceAt(refPos.getContig(),refPos.getStart()-1,refPos.getStart()-1).getBases()[0];
            GenotypesContext GLs = GenotypesContext.create(samples.size());
            final List<Allele> noCall = new ArrayList<Allele>();
            noCall.add(Allele.NO_CALL);
            for ( String s : samples ) {
                //logger.info(gl.getKey() + Arrays.deepToString(ArrayUtils.toObject(gl.getValue())));
                HashMap<String, Object> attributes = new HashMap<String, Object>();
                if ( rawLiks.containsKey(s) ) {
                    attributes.put(VCFConstants.GENOTYPE_PL_KEY, GenotypeLikelihoods.fromLog10Likelihoods(MathUtils.normalizeFromLog10(rawLiks.get(s), false, true)));
                } else {
                    attributes.put(VCFConstants.GENOTYPE_PL_KEY,VCFConstants.MISSING_VALUE_v4);
                }
                GLs.add(GenotypeBuilder.create(s, noCall, attributes));
            }
            Map<String,Object> attributes = new HashMap<String,Object>();
            attributes.put("SNS", likelihoods.getSupNoSupRatio());
            // see if we need to warn on these scores
            if ( hypothesis.getScoreDifference() > 5 ) {
                logger.warn("There is some evidence suggesting an insertion or deletion polymorphism on retrotranscript for " +
                     "hypothesis "+hypothesis.toString()+" . Would suggest running with full smith waterman to avoid misgenotyping");
                attributes.put("FSWF",true);
            }
            double[] prior = computeAlleleFrequencyPriors(GLs.size()*2+1);
            // gls, num alt, priors, result, preserve

            VariantContextBuilder vcb = new VariantContextBuilder("EJG",refPos.getContig(),refPos.getStop(),refPos.getStop(),Arrays.asList(ref,alt));
            vcb.genotypes(GLs);
            List<Allele> alleles = new ArrayList<Allele>(2);
            alleles.add(ref);
            alleles.add(Allele.create(hypothesis.toString()));
            vcb.alleles(alleles);
            VariantContext asCon = vcb.make();
            GenotypesContext genAssigned = GATKVariantContextUtils.assignDiploidGenotypes(asCon);
            vcb.genotypes(genAssigned);

            final AFCalc AFCalculator = AFCalcFactory.createAFCalc(samples.size());
            final AFCalcResult result = AFCalculator.getLog10PNonRef(vcb.make(), prior);

            final double log10POfF0 = result.getLog10PosteriorOfAFGT0();
            logger.debug(log10POfF0);
            vcb.log10PError(log10POfF0);
            attributes.put("MLEAC", result.getAlleleCountsOfMLE()[0]);
            VariantContextUtils.calculateChromosomeCounts(vcb.make(),attributes,false);
            vcb.attributes(attributes);

            vcfWriter.add(vcb.make());
        }

    }

    private double[] computeAlleleFrequencyPriors(int N) {
        // todo -- this stolen from the UG, maybe call into it?
        // calculate the allele frequency priors for 1-N
        double sum = 0.0;
        double heterozygosity = 0.0000001; // ~ 10 differences per person = 10/3gb ~ 1/100 mil
        double[] priors = new double[2*N+1];

        for (int i = 1; i <= 2*N; i++) {
            double value = heterozygosity / (double)i;
            priors[i] = Math.log10(value);
            sum += value;
        }

        // null frequency for AF=0 is (1 - sum(all other frequencies))
        priors[0] = Math.log10(1.0 - sum);
        return priors;
    }

    class EvaluationContext  {
        public GATKSAMRecord read;
        public Map<String,RefSeqFeature> refSeqByName;
        public TreeSet<JunctionHypothesis> hypotheses;

        public void merge(EvaluationContext that) {
            hypotheses.addAll(that.hypotheses);
            refSeqByName.putAll(that.refSeqByName);
            read = that.read;
        }

        public boolean malformed() {
            // reads cam come after the hypotheses they're built from. It's dumb, but a property of the tracker.
            return hypotheses.size() == 0 || ! getToolkit().getGenomeLocParser().createGenomeLoc(read).overlapsP(hypotheses.first().getLocation().getStartLocation().endpointSpan(hypotheses.last().getLocation().getStopLocation()));
        }
    }

    class ECLikelihoods implements HasGenomeLocation {
        public Map<JunctionHypothesis,Map<String,double[]>> likelihoods = new TreeMap<JunctionHypothesis,Map<String,double[]>>();
        public EvaluationContext context;
        private IntronLossGenotypeLikelihoodCalculationModel ilglcm;
        GenomeLoc latestEnd = null;
        private int sup;
        private int noSup;

        public ECLikelihoods(IntronLossGenotypeLikelihoodCalculationModel model) {
            ilglcm = model;
            sup = 0;
            noSup = 0;
        }

        public void purge() {
            likelihoods = new TreeMap<JunctionHypothesis,Map<String,double[]>>();
            context = null;
        }

        public void update(EvaluationContext ec) {
            if ( ! ec.malformed() ) {
                if ( context == null ) {
                    context = ec;
                } else {
                    context.merge(ec);
                }

                if ( context != null && context.hypotheses.size() > likelihoods.size() ) {
                    for ( JunctionHypothesis j : ec.hypotheses ) {
                        GenomeLoc end = j.getLocation().getStopLocation();
                        if ( latestEnd == null || end.isPast(latestEnd) ) {
                            latestEnd = end;
                        }
                    }
                }

                updateLilkelihoods();
            }
        }

        private void updateLilkelihoods() {
            for ( JunctionHypothesis hypothesis : context.hypotheses ) {
                if ( ! likelihoods.containsKey(hypothesis) ) {
                    likelihoods.put(hypothesis,new HashMap<String,double[]>());
                }
                // early return if the read is off the edge of the hypothesized sequence
                boolean offEdge = context.read.getUnclippedStart() < hypothesis.getLocation().getStart() ||
                        context.read.getUnclippedEnd() > hypothesis.getLocation().getStop();
                if ( offEdge ) {
                    continue;
                }
                Pair<String,double[]> lik = ilglcm.getLikelihoods(hypothesis,context.read);
                if ( ! likelihoods.get(hypothesis).containsKey(lik.first) ) {
                    likelihoods.get(hypothesis).put(lik.first,new double[3]);
                }
                //logger.debug("Updating likelihoods for sample "+context.read.getReadGroup().getSample());
                updateLikelihoods(likelihoods.get(hypothesis).get(lik.first), lik.second);
                double[] dlik = likelihoods.get(hypothesis).get(lik.first);
                //logger.debug(String.format("[0] :: %f   [1] :: %f   [2] :: %f",dlik[0],dlik[1],dlik[2]));
            }
            //logger.info(String.format("Context: %d Likelihoods: %d",context.hypotheses.size(),likelihoods.size()));
        }

        private void updateLikelihoods(double[] sum, double[] newLik) {
            double max = Integer.MIN_VALUE;
            int maxIdx = -1;
            for ( int i = 0; i <= 2; i++ ) {
                if ( newLik[i] > max ) {
                    max = newLik[i];
                    maxIdx = i;
                }
                sum[i] += newLik[i];
            }
            if (MathUtils.compareDoubles(newLik[0],newLik[1],1e-5) != 0 ) {
                if ( maxIdx > 0 ) {
                    sup++;
                } else {
                    noSup++;
                }
            }
        }

        public GenomeLoc getLocation() {
            TreeSet<JunctionHypothesis> hypotheses = new TreeSet<JunctionHypothesis>(likelihoods.keySet());
            if ( hypotheses.size() == 0 ) {
                return null;
            }
            return hypotheses.first().getLocation().getStartLocation().endpointSpan(latestEnd);
        }

        public double getSupNoSupRatio() {
            return (0.0+sup)/(noSup);
        }
    }
}
