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

package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.features.annotator.AnnotatorInputTableFeature;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.variant.utils.BaseUtils;
import org.broadinstitute.variant.vcf.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.util.*;

import static org.broadinstitute.sting.utils.variant.GATKVCFUtils.getVCFHeadersFromRods;


/**
 * Walks along all variant ROD loci, and dynamically annotates alleles at MNP records.
 */
@Allows(value = {DataSource.REFERENCE})
@Requires(value = {DataSource.REFERENCE}, referenceMetaData = {@RMD(name = AnnotateMNPsWalker.REFSEQ_ROD_NAME, type = AnnotatorInputTableFeature.class), @RMD(name = AnnotateMNPsWalker.VARIANT_ROD_NAME, type = ReferenceOrderedDatum.class)})

public class AnnotateMNPsWalker extends RodWalker<Integer, Integer> {

    @Output(doc = "File to which variants should be written", required = true)
    protected VCFWriter writer = null;
    private ManualSortingVCFWriter sortingWriter = null;

    @Argument(fullName = "emitOnlyMNPs", shortName = "emitOnlyMNPs", doc = "Only output MNP records; [default:false]", required = false)
    protected boolean emitOnlyMNPs = false;    

    private String rodName = "variant";
    private GenomeLocParser locParser = null;
    private TreeMap<GenomeLoc, Set<GenomeLoc>> MNPstartToStops = null; // Must be TreeMap sorted by START sites!

    public final static String REFSEQ_ROD_NAME = "refseq";
    public final static String VARIANT_ROD_NAME = "variant";

    private LocusToFeatures locusToRefSeqFeatures = null;


    protected final static String MNP_ANNOTATION_KEY_PREFIX = "MNP.refseq.";

    protected final static String REFSEQ_NAME = "name";
    protected final static String REFSEQ_NAME2 = "name2";

    protected final static String REFSEQ_POSITION_TYPE = "positionType";
    protected final static String REFSEQ_CDS = "CDS";

    protected final static String REFSEQ_STRAND = "transcriptStrand";
    protected final static String REFSEQ_POS_STRAND = "+";
    protected final static String REFSEQ_NEG_STRAND = "-";

    protected final static String REFSEQ_CODON_COORD = "codonCoord";
    protected final static String REFSEQ_CODING_FRAME = "frame";

    protected final static String REFSEQ_REF_CODON = "referenceCodon";
    protected final static String REFSEQ_REF_AA = "referenceAA";

    protected final static String REFSEQ_ALT_BASE = "haplotypeAlternate";

    protected final static String REFSEQ_VARIANT_CODON = "variantCodon";
    protected final static String REFSEQ_VARIANT_AA = "variantAA";
    protected final static String REFSEQ_CHANGES_AA = "changesAA";
    protected final static String REFSEQ_FUNCTIONAL_CLASS = "functionalClass";
    protected final static String REFSEQ_PROTEIN_COORD_DESCRIPTION = "proteinCoordStr";

    protected final static String REFSEQ_CODING_ANNOTATIONS = "codingVariants";
    protected final static String REFSEQ_NUM_AA_CHANGES = "numAAchanges";
    protected final static String REFSEQ_HAS_MULT_AA_CHANGES = "alleleHasMultAAchanges";

    public void initialize() {
        locParser = getToolkit().getGenomeLocParser();
        MNPstartToStops = new TreeMap<GenomeLoc, Set<GenomeLoc>>(); // sorted by start sites

        initializeVcfWriter();

        locusToRefSeqFeatures = new LocusToFeatures();
    }

    private void initializeVcfWriter() {
        sortingWriter = new ManualSortingVCFWriter(writer);
        writer = sortingWriter;

        // setup the header fields:
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(GATKVCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));

        Map<String, VCFHeader> rodNameToHeader = getVCFHeadersFromRods(getToolkit(), Arrays.asList(rodName));
        writer.writeHeader(new VCFHeader(hInfo, new TreeSet<String>(rodNameToHeader.get(rodName).getGenotypeSamples())));
    }

    public boolean generateExtendedEvents() {
        return false;
    }

    public Integer reduceInit() {
        return 0;
    }

    /**
     * For each site of interest, annotate it if it's a MNP.
     *
     * @param tracker the meta-data tracker
     * @param ref     the reference base
     * @param context the context for the given locus
     * @return count of MNPs observed
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return null;

        int numMNPsObserved = 0;
        GenomeLoc curLocus = ref.getLocus();
        clearOldLocusFeatures(curLocus);

        for (VariantContext vc : tracker.getValues(VariantContext.class, rodName)) {
            GenomeLoc vcLoc = GATKVariantContextUtils.getLocation(locParser, vc);
            boolean atStartOfVc = curLocus.getStart() == vcLoc.getStart();
            boolean atEndOfVc = curLocus.getStart() == vcLoc.getStop();

            if (vc.isMNP()) {
                logger.debug("Observed MNP at " + vcLoc);

                if (isChrM(vc)) {
                    if (atStartOfVc) {
                        logger.warn("Skipping mitochondrial MNP at " + vcLoc + " due to complexity of coding table [need to know if first codon, etc.]...");
                        writeVCF(vc);
                    }
                    continue;
                }

                GenomeLoc stopLoc = locParser.createGenomeLoc(curLocus.getContig(), vcLoc.getStop());
                final List<Object> refSeqRODs = tracker.getValues(REFSEQ_ROD_NAME);
                for (Object refSeqObject : refSeqRODs) {
                    AnnotatorInputTableFeature refSeqAnnotation = (AnnotatorInputTableFeature) refSeqObject;
                    locusToRefSeqFeatures.putLocusFeatures(curLocus, refSeqAnnotation, stopLoc);
                }

                if (atStartOfVc) { // MNP is starting here, so register that we're waiting for it
                    Set<GenomeLoc> stopLocs = MNPstartToStops.get(curLocus);
                    if (stopLocs == null) {
                        stopLocs = new HashSet<GenomeLoc>();
                        MNPstartToStops.put(curLocus, stopLocs);
                    }
                    stopLocs.add(stopLoc);
                }

                if (atEndOfVc) {
                    numMNPsObserved++; // only count a MNP at its stop site
                    logger.debug("Observed end of MNP at " + curLocus);
                    logger.debug("Current list of per-locus features\n" + locusToRefSeqFeatures);

                    Map<String, Object> MNPannotations = annotateMNP(vc);
                    MNPannotations.putAll(RefSeqDataParser.removeRefSeqAttributes(vc.getAttributes())); // remove any RefSeq INFO, since adding it in more thoroughly here
                    vc = VariantContext.modifyAttributes(vc, MNPannotations);
                    writeVCF(vc);

                    GenomeLoc startLoc = locParser.createGenomeLoc(curLocus.getContig(), vcLoc.getStart());
                    Set<GenomeLoc> stopLocs = MNPstartToStops.get(startLoc);
                    if (stopLocs != null) { // otherwise, just removed stopLocs due to another MNP that has the same (start, stop)
                        stopLocs.remove(stopLoc);
                        if (stopLocs.isEmpty()) // no longer waiting for startLoc
                            MNPstartToStops.remove(startLoc);
                    }
                }
            }
            else if (atStartOfVc && !emitOnlyMNPs) {// only want to write other VariantContexts records once (where they start):
                writeVCF(vc);
            }
        }

        Integer mostUpstreamWritableLoc = null;
        if (!MNPstartToStops.isEmpty()) {
            GenomeLoc waitingForLoc = MNPstartToStops.entrySet().iterator().next().getKey();
            mostUpstreamWritableLoc = waitingForLoc.getStart() - 1;
        }
        sortingWriter.setmostUpstreamWritableLocus(mostUpstreamWritableLoc);

        return numMNPsObserved;
    }

    private static boolean isChrM(final VariantContext vc) {
        return vc.getChr().equals("chrM") || vc.getChr().equals("MT");
    }

    private Map<String, Object> annotateMNP(VariantContext vc) {
        Map<String, Object> annotations = new HashMap<String, Object>();

        RefSeqNameToFeatures nameToPositionalFeatures = new RefSeqNameToFeatures(vc);
        MNPannotationKeyBuilder kb = new MNPannotationKeyBuilder(nameToPositionalFeatures);

        for (Map.Entry<String, RefSeqFeatureList> nameToFeatureEntry : nameToPositionalFeatures.entrySet()) {
            String featureName = nameToFeatureEntry.getKey();
            RefSeqFeatureList feature = nameToFeatureEntry.getValue();
            CodonAnnotationsForAltAlleles codonAnnotationsForAlleles = new CodonAnnotationsForAltAlleles(vc, feature);

            annotations.put(kb.getKey(REFSEQ_CODING_ANNOTATIONS), codonAnnotationsForAlleles.getCodonAnnotationsString());
            annotations.put(kb.getKey(REFSEQ_NUM_AA_CHANGES), codonAnnotationsForAlleles.getNumAAchangesString());
            annotations.put(kb.getKey(REFSEQ_HAS_MULT_AA_CHANGES), codonAnnotationsForAlleles.hasAlleleWithMultipleAAchanges);
            annotations.put(kb.getKey(REFSEQ_NAME), featureName);
            annotations.put(kb.getKey(REFSEQ_NAME2), feature.name2);
            annotations.put(kb.getKey(REFSEQ_POSITION_TYPE), REFSEQ_CDS);
            annotations.put(kb.getKey(REFSEQ_STRAND), (feature.positiveStrand ? REFSEQ_POS_STRAND : REFSEQ_NEG_STRAND));
            annotations.put(kb.getKey(REFSEQ_CODON_COORD), feature.getCodonCoordString());

            kb.incrementFeatureIndex();
        }

        return annotations;
    }

    private static class MNPannotationKeyBuilder {
        private int featureIndex;
        private boolean multipleEntries;

        public MNPannotationKeyBuilder(RefSeqNameToFeatures nameToPositionalFeatures) {
            this.featureIndex = 1;
            this.multipleEntries = nameToPositionalFeatures.nameToFeatures.size() > 1;
        }

        public void incrementFeatureIndex() {
            featureIndex++;
        }

        public String getKey(String type) {
            String annotationKey = MNP_ANNOTATION_KEY_PREFIX + type;
            if (multipleEntries)
                annotationKey += "_" + featureIndex;
            return annotationKey;
        }
    }

    private static byte[] ByteArrayToPrimitive(Byte[] nonNullArray) {
        byte[] primArray = new byte[nonNullArray.length];

        for (int i = 0; i < nonNullArray.length; i++) {
            if (nonNullArray[i] == null)
                throw new ReviewedStingException("nonNullArray[i] == null");
            primArray[i] = nonNullArray[i];
        }

        return primArray;
    }

    private void clearOldLocusFeatures(GenomeLoc curLoc) {
        Iterator<Map.Entry<GenomeLoc, PositionalRefSeqFeatures>> locusFeaturesIt = locusToRefSeqFeatures.entrySet().iterator();
        while (locusFeaturesIt.hasNext()) {
            Map.Entry<GenomeLoc, PositionalRefSeqFeatures> locusFeaturesEntry = locusFeaturesIt.next();
            if (curLoc.isPast(locusFeaturesEntry.getValue().getFurthestLocusUsingFeatures()))
                locusFeaturesIt.remove();
        }
    }

    public Integer reduce(Integer count, Integer total) {
        if (count != null)
            total = total + count;

        return total;
    }

    /**
     * @param result the number of MNPs processed.
     */
    public void onTraversalDone(Integer result) {
        System.out.println("Number of MNPs observed: " + result);
        writer.close();
    }

    private void writeVCF(VariantContext vc) {
        WriteVCF.writeVCF(vc, writer, logger);
    }

    /*
     Inner classes:
     */

    // Maps: RefSeq entry name -> features for ALL positions of a particular VariantContext MNP:

    private class RefSeqNameToFeatures {
        private Map<String, RefSeqFeatureList> nameToFeatures;

        public RefSeqNameToFeatures(VariantContext vc) {
            this.nameToFeatures = new HashMap<String, RefSeqFeatureList>();

            int MNPstart = vc.getStart();
            int MNPstop = vc.getEnd();
            int MNPlength = MNPstop - MNPstart + 1;

            for (int i = 0; i < MNPlength; i++) {
                int genomicPosition = MNPstart + i;
                GenomeLoc posLoc = locParser.createGenomeLoc(vc.getChr(), genomicPosition);

                PositionalRefSeqFeatures locFeatures = locusToRefSeqFeatures.getLocusFeatures(posLoc);
                if (locFeatures == null) // no features for posLoc
                    continue;

                for (Map.Entry<String, PositionalRefSeqFeature> nameToFeatureEntry : locFeatures.entrySet()) {
                    String name = nameToFeatureEntry.getKey();
                    PositionalRefSeqFeature posFeature = nameToFeatureEntry.getValue();

                    RefSeqFeatureList featureList = nameToFeatures.get(name);
                    if (featureList == null) {
                        featureList = new RefSeqFeatureList(MNPlength);
                        nameToFeatures.put(name, featureList);
                    }
                    featureList.updateFeatureAtPosition(i, posFeature);
                }
            }
        }

        public Set<Map.Entry<String, RefSeqFeatureList>> entrySet() {
            return nameToFeatures.entrySet();
        }
    }

    // For a particular RefSeq entry, contains the features for ALL positions of a particular VariantContext MNP

    private static class RefSeqFeatureList {
        private final static String CODON_FRAME_START = "(";
        private final static String CODON_FRAME_END = ")";
        private final static String CODON_DELIM = "|";

        private CodingRefSeqFeature[] refSeqFeatures;
        private String name2;
        private Boolean positiveStrand;

        private Map<Integer, List<Integer>> codonToIndices; // Map of: codon index -> MNP indices that refer to codon

        public RefSeqFeatureList(int MNPlength) {
            this.refSeqFeatures = new CodingRefSeqFeature[MNPlength];
            for (int i = 0; i < MNPlength; i++)
                this.refSeqFeatures[i] = null;

            this.name2 = null;
            this.positiveStrand = null;
            this.codonToIndices = new TreeMap<Integer, List<Integer>>();
        }

        public void updateFeatureAtPosition(int index, PositionalRefSeqFeature feature) {
            if (name2 == null) {
                name2 = feature.name2;
                positiveStrand = feature.positiveStrand;
            }
            else if (!name2.equals(feature.name2) || positiveStrand != feature.positiveStrand) {
                throw new UserException("Inconsistency between previous RefSeq entry and: " + feature);
            }

            CodingRefSeqFeature crsf = new CodingRefSeqFeature(feature);
            refSeqFeatures[index] = crsf;

            List<Integer> indicesWithCodon = codonToIndices.get(crsf.codonCoord);
            if (indicesWithCodon == null) {
                indicesWithCodon = new LinkedList<Integer>();
                codonToIndices.put(crsf.codonCoord, indicesWithCodon);
            }
            indicesWithCodon.add(index);
        }

        public Set<Map.Entry<Integer, List<Integer>>> codonIndicesEntrySet() {
            return codonToIndices.entrySet();
        }

        public String getCodonCoordString() {
            StringBuilder sb = new StringBuilder();

            for (int i = 0; i < refSeqFeatures.length; i++) {
                CodingRefSeqFeature crsf = refSeqFeatures[i];
                if (crsf != null)
                    sb.append(crsf.codonCoord).append(CODON_FRAME_START).append(crsf.codingFrame).append(CODON_FRAME_END);
                if (i < refSeqFeatures.length - 1)
                    sb.append(CODON_DELIM);
            }

            return sb.toString();
        }
    }

    private static class CodingRefSeqFeature {
        protected int codonCoord;
        protected int codingFrame;
        protected String referenceCodon;
        protected String referenceAA;

        public CodingRefSeqFeature(PositionalRefSeqFeature feature) {
            this.codonCoord = feature.codonCoord;
            this.codingFrame = feature.codingFrame;
            this.referenceCodon = feature.referenceCodon.toUpperCase();
            this.referenceAA = feature.referenceAA;
        }
    }

    private static class CodonAnnotationsForAltAlleles {
        protected final static int MIN_CODON_INDEX = 0;
        protected final static int NUM_CODON_INDICES = 3;
        private final static String CODON_ANNOTATION_DELIM = ",";

        private List<SingleCodonAnnotationsForAlleles> alleleAnnotations;
        private int[] alleleToNumAAchanges;
        private boolean hasAlleleWithMultipleAAchanges;

        public CodonAnnotationsForAltAlleles(VariantContext vc, RefSeqFeatureList feature) {
            this.alleleAnnotations = new LinkedList<SingleCodonAnnotationsForAlleles>();

            Set<Allele> altAlleles = vc.getAlternateAlleles();
            int numAltAlleles = altAlleles.size();
            this.alleleToNumAAchanges = new int[numAltAlleles];
            for (int i = 0; i < numAltAlleles; i++)
                this.alleleToNumAAchanges[i] = 0;

            int MNPstart = vc.getStart();
            int MNPstop = vc.getEnd();
            int MNPlength = MNPstop - MNPstart + 1;

            for (Map.Entry<Integer, List<Integer>> codonToIndicesEntry : feature.codonIndicesEntrySet()) {
                int codonIndex = codonToIndicesEntry.getKey();
                List<Integer> indices = codonToIndicesEntry.getValue();
                if (indices.isEmpty())
                    throw new ReviewedStingException("indices should not exist if it's empty!");

                for (int index : indices) {
                    int frame = feature.refSeqFeatures[index].codingFrame;
                    if (feature.refSeqFeatures[index].codonCoord != codonIndex)
                        throw new ReviewedStingException("LOGICAL ERROR: feature.refSeqFeatures[index].codonCoord != codonIndex");
                    if (frame < MIN_CODON_INDEX || frame >= NUM_CODON_INDICES)
                        throw new UserException("RefSeq codon frame not one of {0,1,2}");
                }
                CodingRefSeqFeature firstFeatureForCodon = feature.refSeqFeatures[indices.get(0)];
                String refCodon = firstFeatureForCodon.referenceCodon;

                SingleCodonAnnotationsForAlleles codonAnnotation = new SingleCodonAnnotationsForAlleles(codonIndex, altAlleles, MNPlength, refCodon, firstFeatureForCodon, indices, feature);
                alleleAnnotations.add(codonAnnotation);

                // From a single codon, summarize the data for ALL alleles:
                for (int i = 0; i < numAltAlleles; i++) {
                    if (codonAnnotation.annotationsForAlleles[i].codonFunc.changesAA) {
                        alleleToNumAAchanges[i]++;
                        if (alleleToNumAAchanges[i] > 1)
                            this.hasAlleleWithMultipleAAchanges = true;
                    }
                }
            }
        }

        public String getCodonAnnotationsString() {
            StringBuilder sb = new StringBuilder();

            int index = 0;
            for (SingleCodonAnnotationsForAlleles codonToAlleles : alleleAnnotations) {
                sb.append(codonToAlleles);
                if (index < alleleAnnotations.size() - 1)
                    sb.append(CODON_ANNOTATION_DELIM);
                index++;
            }

            return sb.toString();
        }

        public String getNumAAchangesString() {
            StringBuilder sb = new StringBuilder();

            for (int index = 0; index < alleleToNumAAchanges.length; index++) {
                sb.append(alleleToNumAAchanges[index]);
                if (index < alleleToNumAAchanges.length - 1)
                    sb.append(SingleCodonAnnotationsForAlleles.ALLELE_ANNOTATION_DELIM);
            }

            return sb.toString();
        }
    }

    private static class SingleCodonAnnotationsForAlleles {
        private final static String CODON_MAP_SYMBOL = "->";
        private final static String CODON_ANNOTATION_START = "[";
        private final static String CODON_ANNOTATION_END = "]";
        private final static String REF_CODON_INFO_DELIM = "|";
        private final static String ALLELE_ANNOTATION_DELIM = ",";
        private final static String ASSIGNMENT = ":";

        private int codonIndex;
        private String refCodon;
        private String refAA;

        private SingleCodonAnnotationsForAllele[] annotationsForAlleles;

        public SingleCodonAnnotationsForAlleles(int codonIndex, Collection<Allele> altAlleles, int MNPlength, String refCodon, CodingRefSeqFeature firstFeatureForCodon, List<Integer> indices, RefSeqFeatureList feature) {
            if (refCodon.length() != CodonAnnotationsForAltAlleles.NUM_CODON_INDICES)
                throw new UserException("RefSeq reference codon " + refCodon + " is not of length " + CodonAnnotationsForAltAlleles.NUM_CODON_INDICES);

            AminoAcid refAA = AminoAcidTable.getEukaryoticAA(refCodon);
            if (!refAA.getCode().equals(firstFeatureForCodon.referenceAA))
                throw new UserException("RefSeq: translated reference codon= " + refAA + " != " + firstFeatureForCodon.referenceAA + " = reference AA");

            this.codonIndex = codonIndex;
            this.refCodon = refCodon;
            this.refAA = refAA.getCode();
            this.annotationsForAlleles = new SingleCodonAnnotationsForAllele[altAlleles.size()];

            int altInd = 0;
            for (Allele altAllele : altAlleles) {
                if (altAllele.length() != MNPlength)
                    throw new ReviewedStingException("length(altAllele) != length(MNP)");
                byte[] altBases = altAllele.getBases();

                Byte[] variantCodonArr = new Byte[CodonAnnotationsForAltAlleles.NUM_CODON_INDICES];
                for (int i = CodonAnnotationsForAltAlleles.MIN_CODON_INDEX; i < CodonAnnotationsForAltAlleles.NUM_CODON_INDICES; i++)
                    variantCodonArr[i] = null;

                for (int index : indices) {
                    int frame = feature.refSeqFeatures[index].codingFrame;
                    if (variantCodonArr[frame] != null)
                        throw new UserException("RefSeq assigns codon " + codonIndex + " twice at same frame: " + frame);

                    byte base = altBases[index];
                    if (!feature.positiveStrand) // negative strand codon
                        base = BaseUtils.simpleComplement(base);

                    variantCodonArr[frame] = base;
                }

                /* For missing frames, there MUST exist AT LEAST one index that refers to this codon,
                  so use it to derive the missing bases [ALREADY complemented if on the negative strand]:
                */
                for (int frame = CodonAnnotationsForAltAlleles.MIN_CODON_INDEX; frame < CodonAnnotationsForAltAlleles.NUM_CODON_INDICES; frame++) {
                    if (variantCodonArr[frame] == null)
                        variantCodonArr[frame] = (byte) refCodon.charAt(frame);
                }
                String variantCodon = new String(ByteArrayToPrimitive(variantCodonArr)).toUpperCase();

                SingleCodonAnnotationsForAllele alleleAnnotation = new SingleCodonAnnotationsForAllele(variantCodon, refCodon, refAA, codonIndex);
                annotationsForAlleles[altInd] = alleleAnnotation;
                altInd++;
            }
        }

        public String toString() {
            StringBuilder sb = new StringBuilder();

            sb.append(codonIndex).append(CODON_MAP_SYMBOL).append(CODON_ANNOTATION_START);
            sb.append(REFSEQ_REF_CODON).append(ASSIGNMENT).append(refCodon).append(REF_CODON_INFO_DELIM);
            sb.append(REFSEQ_REF_AA).append(ASSIGNMENT).append(refAA).append(REF_CODON_INFO_DELIM);

            int index = 0;
            for (SingleCodonAnnotationsForAllele annotation : annotationsForAlleles) {
                sb.append(annotation);
                if (index < annotationsForAlleles.length - 1)
                    sb.append(ALLELE_ANNOTATION_DELIM);
                index++;
            }
            sb.append(CODON_ANNOTATION_END);

            return sb.toString();
        }
    }

    private static class SingleCodonAnnotationsForAllele {
        private final static String ALLELE_START = "{";
        private final static String ALLELE_END = "}";
        private final static String CODON_INFO_DELIM = "|";
        private final static String ASSIGNMENT = ":";
        private final static String MNP_DEPENDENT_AA = "MNPdependentAA";

        private CodonFunction codonFunc;
        private String proteinCoordStr;
        private boolean MNPdependentAA;
        private String originalAA;

        public SingleCodonAnnotationsForAllele(String variantCodon, String refCodon, AminoAcid refAA, int codonIndex) {
            this.codonFunc = new CodonFunction(variantCodon, refCodon, refAA);
            this.proteinCoordStr = "p." + refAA.getLetter() + codonIndex + codonFunc.variantAA.getLetter();

            int refCodonLength = refCodon.length();
            if (codonFunc.variantCodon.length() != refCodonLength)
                throw new ReviewedStingException("codonFunc.variantCodon.length() != refCodonLength, but ALREADY checked that they're both 3");

            this.MNPdependentAA = true;
            this.originalAA = "(";
            for (int i = 0; i < refCodonLength; i++) {
                // Take [0,i-1] and [i+1, end] from refCodon, and i from variantCodon:
                String singleBaseChangeCodon = refCodon.substring(0, i) + variantCodon.substring(i, i+1) + refCodon.substring(i+1, refCodonLength);
                CodonFunction singleBaseChangeCodonFunc = new CodonFunction(singleBaseChangeCodon, refCodon, refAA);
                if (singleBaseChangeCodonFunc.variantAA.equals(codonFunc.variantAA)) {
                    this.MNPdependentAA = false;
                    this.originalAA = "";
                    break;
                }

                this.originalAA = this.originalAA + "" + singleBaseChangeCodonFunc.variantAA.getLetter();
                if (i < refCodonLength - 1)
                    this.originalAA = this.originalAA + ",";
            }

            if (this.MNPdependentAA)
                this.originalAA = this.originalAA + ")";
        }

        private static class CodonFunction {
            private String variantCodon;
            private AminoAcid variantAA;
            private boolean changesAA;
            private String functionalClass;

            public CodonFunction(String variantCodon, String refCodon, AminoAcid refAA) {
                this.variantCodon = variantCodon;
                this.variantAA = AminoAcidTable.getEukaryoticAA(this.variantCodon);
                this.changesAA = !refAA.equals(variantAA);

                if (!this.variantCodon.equals(refCodon)) {
                    if (changesAA) {
                        if (variantAA.isStop()) {
                            functionalClass = "nonsense";
                        }
                        else if (refAA.isStop()) {
                            functionalClass = "readthrough";
                        }
                        else {
                            functionalClass = "missense";
                        }
                    }
                    else { // the same aa:
                        functionalClass = "silent";
                    }
                }
                else { // the same codon:
                    functionalClass = "no_change";
                }
            }
        }

        public String toString() {
            StringBuilder sb = new StringBuilder();

            sb.append(ALLELE_START);
            sb.append(REFSEQ_VARIANT_CODON).append(ASSIGNMENT).append(codonFunc.variantCodon).append(CODON_INFO_DELIM);
            sb.append(REFSEQ_VARIANT_AA).append(ASSIGNMENT).append(codonFunc.variantAA.getCode()).append(CODON_INFO_DELIM);
            sb.append(REFSEQ_CHANGES_AA).append(ASSIGNMENT).append(codonFunc.changesAA).append(CODON_INFO_DELIM);
            sb.append(REFSEQ_FUNCTIONAL_CLASS).append(ASSIGNMENT).append(codonFunc.functionalClass).append(CODON_INFO_DELIM);
            sb.append(REFSEQ_PROTEIN_COORD_DESCRIPTION).append(ASSIGNMENT).append(proteinCoordStr).append(CODON_INFO_DELIM);
            sb.append(MNP_DEPENDENT_AA).append(ASSIGNMENT).append(MNPdependentAA).append(originalAA);
            sb.append(ALLELE_END);

            return sb.toString();
        }
    }
}


// External classes:

class LocusToFeatures {
    private Map<GenomeLoc, PositionalRefSeqFeatures> locusToFeatures;

    public LocusToFeatures() {
        this.locusToFeatures = new TreeMap<GenomeLoc, PositionalRefSeqFeatures>();
    }

    public PositionalRefSeqFeatures getLocusFeatures(GenomeLoc loc) {
        return locusToFeatures.get(loc);
    }

    public void putLocusFeatures(GenomeLoc loc, AnnotatorInputTableFeature refSeqAnnotation, GenomeLoc locusUsingThis) {
        PositionalRefSeqFeatures locFeatures = locusToFeatures.get(loc);
        if (locFeatures == null) {
            locFeatures = new PositionalRefSeqFeatures(locusUsingThis);
            locusToFeatures.put(loc, locFeatures);
        }
        locFeatures.putFeature(refSeqAnnotation, locusUsingThis);
    }

    public Set<Map.Entry<GenomeLoc, PositionalRefSeqFeatures>> entrySet() {
        return locusToFeatures.entrySet();
    }

    public String toString() { // INTERNAL use only
        StringBuilder sb = new StringBuilder();

        for (Map.Entry<GenomeLoc, PositionalRefSeqFeatures> locFeatures : entrySet()) {
            GenomeLoc loc = locFeatures.getKey();
            PositionalRefSeqFeatures features = locFeatures.getValue();
            sb.append("Locus: ").append(loc).append("\n").append(features);
        }

        return sb.toString();
    }
}

class PositionalRefSeqFeatures {
    private final static String[] REQUIRE_COLUMNS =
            {AnnotateMNPsWalker.REFSEQ_NAME, AnnotateMNPsWalker.REFSEQ_POSITION_TYPE};

    private Map<String, PositionalRefSeqFeature> nameToFeature;
    private GenomeLoc furthestLocusUsingFeatures;

    public PositionalRefSeqFeatures(GenomeLoc locusUsingThis) {
        this.nameToFeature = new HashMap<String, PositionalRefSeqFeature>();
        this.furthestLocusUsingFeatures = locusUsingThis;
    }

    public void putFeature(AnnotatorInputTableFeature refSeqAnnotation, GenomeLoc locusUsingThis) {
        for (String column : REQUIRE_COLUMNS) {
            if (!refSeqAnnotation.containsColumnName(column))
                throw new UserException("In RefSeq: " + refSeqAnnotation + " Missing column " + column);
        }

        if (locusUsingThis.isPast(furthestLocusUsingFeatures))
            furthestLocusUsingFeatures = locusUsingThis;

        String posType = refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_POSITION_TYPE);
        if (!posType.equals(AnnotateMNPsWalker.REFSEQ_CDS)) // only interested in coding sequence annotations
            return;

        PositionalRefSeqFeature newLocusFeature = new PositionalRefSeqFeature(refSeqAnnotation);

        String refSeqName = refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_NAME);
        PositionalRefSeqFeature locusFeature = nameToFeature.get(refSeqName);
        if (locusFeature == null) {
            locusFeature = newLocusFeature;
            nameToFeature.put(refSeqName, locusFeature);
        }
        else if (!locusFeature.equals(newLocusFeature)) {
            throw new UserException("Inconsistency between previous RefSeq entry and: " + refSeqAnnotation);
        }

        locusFeature.updateFeature(refSeqAnnotation);
    }

    public GenomeLoc getFurthestLocusUsingFeatures() {
        return furthestLocusUsingFeatures;
    }

    public Set<Map.Entry<String, PositionalRefSeqFeature>> entrySet() {
        return nameToFeature.entrySet();
    }

    public String toString() { // INTERNAL use only
        StringBuilder sb = new StringBuilder();

        for (Map.Entry<String, PositionalRefSeqFeature> nameFeatureEntry : entrySet()) {
            String name = nameFeatureEntry.getKey();
            PositionalRefSeqFeature feature = nameFeatureEntry.getValue();
            sb.append(name).append(" -> [").append(feature).append("]\n");
        }

        return sb.toString();
    }
}

class PositionalRefSeqFeature {
    private final static String[] REQUIRE_COLUMNS =
            {AnnotateMNPsWalker.REFSEQ_NAME2, AnnotateMNPsWalker.REFSEQ_STRAND,
                    AnnotateMNPsWalker.REFSEQ_CODON_COORD, AnnotateMNPsWalker.REFSEQ_CODING_FRAME,
                    AnnotateMNPsWalker.REFSEQ_REF_CODON, AnnotateMNPsWalker.REFSEQ_REF_AA};

    protected String name2;
    protected boolean positiveStrand;
    protected int codonCoord;
    protected int codingFrame;
    protected String referenceCodon;
    protected String referenceAA;

    private Map<String, BaseAnnotations> baseToAnnotations;

    public PositionalRefSeqFeature(AnnotatorInputTableFeature refSeqAnnotation) {
        for (String column : REQUIRE_COLUMNS) {
            if (!refSeqAnnotation.containsColumnName(column))
                throw new UserException("In RefSeq: " + refSeqAnnotation + " Missing column " + column);
        }
        this.name2 = refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_NAME2);
        this.positiveStrand = (refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_STRAND).equals(AnnotateMNPsWalker.REFSEQ_POS_STRAND));
        this.codonCoord = Integer.parseInt(refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_CODON_COORD));
        this.codingFrame = Integer.parseInt(refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_CODING_FRAME));
        this.referenceCodon = refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_REF_CODON);
        this.referenceAA = refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_REF_AA);

        this.baseToAnnotations = new HashMap<String, BaseAnnotations>();
    }

    public boolean equals(PositionalRefSeqFeature that) {
        return this.name2.equals(that.name2) && this.positiveStrand == that.positiveStrand && this.codonCoord == that.codonCoord && this.codingFrame == that.codingFrame
                && this.referenceCodon.equals(that.referenceCodon) && this.referenceAA.equals(that.referenceAA);
    }

    public void updateFeature(AnnotatorInputTableFeature refSeqAnnotation) {
        if (!refSeqAnnotation.containsColumnName(AnnotateMNPsWalker.REFSEQ_ALT_BASE))
            throw new UserException("In RefSeq: " + refSeqAnnotation + " Missing column " + AnnotateMNPsWalker.REFSEQ_ALT_BASE);
        String base = refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_ALT_BASE);

        baseToAnnotations.put(base, new BaseAnnotations(refSeqAnnotation));
    }

    public String toString() { // INTERNAL use only
        StringBuilder sb = new StringBuilder();

        sb.append("name2= ").append(name2);
        sb.append(", positiveStrand= ").append(positiveStrand);
        sb.append(", codonCoord= ").append(codonCoord);
        sb.append(", codingFrame= ").append(codingFrame);
        sb.append(", referenceCodon= ").append(referenceCodon);
        sb.append(", referenceAA= ").append(referenceAA);

        sb.append(", baseAnnotations= {");
        for (Map.Entry<String, BaseAnnotations> baseToAnnotationsEntry : baseToAnnotations.entrySet()) {
            String base = baseToAnnotationsEntry.getKey();
            BaseAnnotations annotations = baseToAnnotationsEntry.getValue();
            sb.append(" ").append(base).append(" -> {").append(annotations).append("}");
        }
        sb.append(" }");

        return sb.toString();
    }
}

class BaseAnnotations {
    private final static String[] REQUIRE_COLUMNS =
            {AnnotateMNPsWalker.REFSEQ_VARIANT_CODON, AnnotateMNPsWalker.REFSEQ_VARIANT_AA,
                    AnnotateMNPsWalker.REFSEQ_CHANGES_AA, AnnotateMNPsWalker.REFSEQ_FUNCTIONAL_CLASS,
                    AnnotateMNPsWalker.REFSEQ_PROTEIN_COORD_DESCRIPTION};

    protected String variantCodon;
    protected String variantAA;
    protected boolean changesAA;
    protected String functionalClass;
    protected String proteinCoordStr;

    public BaseAnnotations(AnnotatorInputTableFeature refSeqAnnotation) {
        for (String column : REQUIRE_COLUMNS) {
            if (!refSeqAnnotation.containsColumnName(column))
                throw new UserException("In RefSeq: " + refSeqAnnotation + " Missing column " + column);
        }
        this.variantCodon = refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_VARIANT_CODON);
        this.variantAA = refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_VARIANT_AA);
        this.changesAA = Boolean.parseBoolean(refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_CHANGES_AA));
        this.functionalClass = refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_FUNCTIONAL_CLASS);
        this.proteinCoordStr = refSeqAnnotation.getColumnValue(AnnotateMNPsWalker.REFSEQ_PROTEIN_COORD_DESCRIPTION);
    }


    public String toString() { // INTERNAL use only
        StringBuilder sb = new StringBuilder();

        sb.append("variantCodon= ").append(variantCodon);
        sb.append(", variantAA= ").append(variantAA);
        sb.append(", changesAA= ").append(changesAA);
        sb.append(", functionalClass= ").append(functionalClass);
        sb.append(", proteinCoordStr= ").append(proteinCoordStr);

        return sb.toString();
    }
}