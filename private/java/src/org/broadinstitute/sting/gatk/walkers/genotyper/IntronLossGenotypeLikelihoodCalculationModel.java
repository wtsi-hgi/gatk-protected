package org.broadinstitute.sting.gatk.walkers.genotyper;

import com.google.java.contract.Requires;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.codecs.refseq.RefSeqFeature;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import java.util.*;

/**
 *
 */
public class IntronLossGenotypeLikelihoodCalculationModel {

    public static final int PLOIDY = 2;
    public static final int ALIGNMENT_QUAL_CAP = 500;
    private static final double GARBAGE_COLLECT_THRESHOLD = -20; // Q200
    private GenomeLocParser genomeLocParser;
    private Map<String,byte[]> insertSizeQualHistograms;
    private Set<String> samples;
    private Logger logger;

    public IntronLossGenotypeLikelihoodCalculationModel(Logger logger) {
        this.logger = logger;
    }

    public void initialize(RodBinding<RefSeqFeature> refSeqBinding,
                           GenomeLocParser parser, Map<String,byte[]> histograms) {
        this.genomeLocParser = parser;
        this.insertSizeQualHistograms = histograms;
    }

    public void setSamples(Set<String> s) {
        samples =s;
    }

    public GenomeLoc getGenomeLoc(GATKSAMRecord read) {
        if ( read == null )
            return null;
        int start =  ReadUtils.getRefCoordSoftUnclippedStart(read);
        int end =  Math.max(start, ReadUtils.getRefCoordSoftUnclippedEnd(read));
        return genomeLocParser.createGenomeLoc(read.getReferenceName(),start,end);
    }

    public List<VariantContext> getLikelihoods(IntronLossGenotyperV2.PairedReadBin readBin, RefSeqFeature geneFeature, IndexedFastaSequenceFile refFile) {
        Map<String,double[]> samLikelihoods = new HashMap<String,double[]>(samples.size());
        for ( String s : samples ) {
            samLikelihoods.put(s,new double[1+PLOIDY]);
        }
        GenomeLoc featureStop = genomeLocParser.createGenomeLoc(geneFeature);

        Map<String,IntronLossPostulate> postulates = new HashMap<String,IntronLossPostulate>(36);
        // first postulate events
        for ( Pair<GATKSAMRecord,GATKSAMRecord> pair : readBin ) {
            IntronLossPostulate p = new IntronLossPostulate(geneFeature.getExons(),getGenomeLoc(pair.first),getGenomeLoc(pair.second));
            if ( ! postulates.containsKey(p.toString()) ) {
                postulates.put(p.toString(),p);
            }
            postulates.get(p.toString()).supportingPairs.add(pair);
        }

        List<VariantContext> vcontexts = new ArrayList<VariantContext>(Math.max(0,postulates.size()-1));

        for  ( IntronLossPostulate postulate : postulates.values() ) {
            if ( postulate.clazz != PostulateClass.NONE && postulate.isValid()) {
                GenomeLoc eventLoc = geneFeature.getLocation();
                Allele ref = Allele.create(refFile.getSubsequenceAt(featureStop.getContig(), featureStop.getStop(), featureStop.getStop()).getBases(), true);
                Allele alt = Allele.create(String.format("<:%s:>",postulate),false);
                Map<String,double[]> newLikelihoods = createPostulateLikelihoodsNew(samples,postulate,postulates.get("NONE"),geneFeature,refFile);
                Map<String,MultiallelicGenotypeLikelihoods> GLs = new HashMap<String,MultiallelicGenotypeLikelihoods>(samples.size());
                for ( Map.Entry<String,double[]> gl : samLikelihoods.entrySet() ) {
                    GLs.put(gl.getKey(), new MultiallelicGenotypeLikelihoods(gl.getKey(),new ArrayList<Allele>(Arrays.asList(ref,alt)),gl.getValue(),postulate.supportingPairs.size()));
                }
                Map<String,Object> attributes = new HashMap<String,Object>();
                HashMap<String, Genotype> genotypes = new HashMap<String, Genotype>();
                for ( String s : samples ) {
                    Map<String,Object> genAttribs = new HashMap<String,Object>();
                    GenotypeLikelihoods likelihoods = GenotypeLikelihoods.fromLog10Likelihoods(newLikelihoods.get(s));
                    genAttribs.put(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY,likelihoods);
                    genotypes.put(s,new Genotype(s, Arrays.asList(Allele.NO_CALL), Genotype.NO_NEG_LOG_10PERROR, null, genAttribs, false));
                }

                attributes.put("SR",postulate.supportingPairs.size());
                attributes.put("GN",geneFeature.getGeneName());
                attributes.put("EL",postulate.exonLocs.first.getStopLocation().toString());

                vcontexts.add(new VariantContext("ILGV2",featureStop.getContig(),featureStop.getStop(),featureStop.getStop(),Arrays.asList(ref,alt),genotypes,-1,null,attributes,ref.getBases()[0]));
            }
        }

        return vcontexts;
    }

    public Map<String,double[]> createPostulateLikelihoodsNew(Set<String> samples, IntronLossPostulate myPostulate, IntronLossPostulate nonePostulate, RefSeqFeature geneFeature, IndexedFastaSequenceFile refFile) {
        Map<String,double[]> samLikelihoods = new HashMap<String,double[]>(samples.size());
        for ( String s : samples ) {
            samLikelihoods.put(s,new double[1+PLOIDY]);
        }
        // note: reads at this point have an accurate NM tag

        byte[] junctionSeq = myPostulate.getSequence(refFile);


        for ( Pair<GATKSAMRecord,GATKSAMRecord> supportingPair : myPostulate.supportingPairs ) {
            double lossScore = 0.0;
            double noLossScore = 0.0;
            double realignmentScore = scoreRealignment(supportingPair,junctionSeq);
            double initialScore = scoreAlignment(supportingPair);
            double rawInsertScore = insertSizeProbability(supportingPair.first,0);
            double adjInsertScore = insertSizeProbability(supportingPair.second,myPostulate.getEventSize());
            lossScore += realignmentScore;
            lossScore += adjInsertScore;
            noLossScore += initialScore;
            noLossScore += rawInsertScore;
            for ( int ac = 0; ac <= PLOIDY; ac ++ ) {
                double logAc = Math.log10(ac);
                double logAnMAc = Math.log10(PLOIDY-ac);
                double logPloidy = Math.log10(PLOIDY); // todo -- should be global constant
                // todo -- numerically unstable: delta can be negative infinity
                double delta = MathUtils.log10sumLog10(new double[]{logAc + lossScore - logPloidy, logAnMAc + noLossScore - logPloidy});
                String sam = supportingPair.first.getReadGroup().getSample();
                samLikelihoods.get(sam)[ac] += delta;
            }
        }

        // Note: Early abort if we are more likely to be ref from SUPPORTING READS
        boolean probRef = true;
        for ( double[] lik : samLikelihoods.values() ) {
            probRef &= MathUtils.compareDoubles(MathUtils.arrayMax(lik),lik[0],1e-5) == 0;
        }
        if ( probRef ) {
            return samLikelihoods;
        }

        if ( nonePostulate == null ) {
            return samLikelihoods;
        }

        for ( Pair<GATKSAMRecord,GATKSAMRecord> nonePair : nonePostulate.supportingPairs ) {
            double lossScore = 0.0;
            double noLossScore = 0.0;
            if ( myPostulate.isRelevant(getGenomeLoc(nonePair.first),getGenomeLoc(nonePair.second)) ) {
                double realignmentScore = scoreRealignment(nonePair,junctionSeq);
                double initialScore = scoreAlignment(nonePair);
                lossScore += realignmentScore;
                noLossScore += initialScore;
                if ( lossScore < GARBAGE_COLLECT_THRESHOLD && noLossScore < GARBAGE_COLLECT_THRESHOLD ) {
                    continue;
                }
                for ( int ac = 0; ac <= PLOIDY; ac ++ ) {
                    double logAc = Math.log10(ac);
                    double logAnMAc = Math.log10(PLOIDY-ac);
                    double logPloidy = Math.log10(PLOIDY); // todo -- should be global constant
                    // todo -- numerically unstable: delta can be negative infinity
                    double delta = MathUtils.log10sumLog10(new double[]{logAc + lossScore - logPloidy, logAnMAc + noLossScore - logPloidy});
                    String sam = nonePair.first.getReadGroup().getSample();
                    samLikelihoods.get(sam)[ac] += delta;
                }
            }
        }

        return samLikelihoods;
    }

    private double scoreAlignment(Pair<GATKSAMRecord,GATKSAMRecord> initiallyAlignedPair) {
        // aligned reads are equipped already with necessary metrics
        double logScore = 0.0;
        for ( GATKSAMRecord read : Arrays.asList(initiallyAlignedPair.first,initiallyAlignedPair.second) ) {
            if ( read == null )
                continue;
            logScore += ( (Integer) read.getAttribute("NM") )*(-1.2);
            for ( CigarElement e : read.getCigar().getCigarElements() ) {
                if ( e.getOperator().equals(CigarOperator.D) || e.getOperator().equals(CigarOperator.I) ) {
                    logScore += -1*(2.0-Math.pow(2,1-e.getLength()));
                } else if ( e.getOperator().equals(CigarOperator.S) ) {
                    logScore += -0.50*e.getLength();
                }
            }
        }

        return logScore;
    }

    private double scoreRealignment(Pair<GATKSAMRecord,GATKSAMRecord> inPair, byte[] sequence) {
        double logScore = 0.0;
        for ( GATKSAMRecord rec : Arrays.asList(inPair.first,inPair.second) ) {
            if ( rec == null )
                continue;
            SWPairwiseAlignment alignment = new SWPairwiseAlignment(sequence,rec.getReadBases());
            logScore += scoreRealignment(alignment.getCigar().getCigarElements(),alignment.getAlignmentStart2wrt1(),rec.getReadBases(),sequence);
        }

        return logScore;
    }

    private double scoreRealignment(Collection<CigarElement> elements, int offset, byte[] readBases, byte[] exonBases) {
        int b = 0;
        int phredBasesAreExonBases = 0;
        int indelPen = 15;
        for ( CigarElement elem : elements) {
            if ( elem.getOperator().equals(CigarOperator.M) ) {
                for ( int ct = elem.getLength(); ct > 0; ct-- ) {
                    int pos = b + offset;
                    if (BaseUtils.basesAreEqual(readBases[b],exonBases[pos])) {
                        // bases match, don't do anything
                    } else {
                        // read clipped at Q15, with e/3 this goes to Q12
                        phredBasesAreExonBases += 12;
                    }
                    b++;
                }
            } else if ( elem.getOperator().equals(CigarOperator.D) ) {
                // just increment offset
                offset += elem.getLength();
                phredBasesAreExonBases += indelPen + 15*elem.getLength();
                indelPen += 15;
            } else if (elem.getOperator().equals(CigarOperator.I)) {
                // update the b
                b += elem.getLength();
                offset -= elem.getLength(); // compensate, since pos will move up by elem.getLength()
                phredBasesAreExonBases += indelPen+15*elem.getLength();
                indelPen += 15;
            } else if ( elem.getOperator().equals(CigarOperator.S) ) {
                // just skip the bases and update b; soft clipping occurs at edges, so need to subtract from offset
                b += elem.getLength();
                offset -= elem.getLength();
                phredBasesAreExonBases += 10*elem.getLength();
            } else if ( elem.getOperator().equals(CigarOperator.H) ) {
                b += elem.getLength();
                offset -= elem.getLength();
            } else {
                throw new StingException("Unsupported operator: "+elem.getOperator().toString());
            }
        }

        int cappedPhred = phredBasesAreExonBases > ALIGNMENT_QUAL_CAP  ? ALIGNMENT_QUAL_CAP : phredBasesAreExonBases;

        return  ( (double) cappedPhred)/(-10.0);
    }
    private double insertSizeProbability(GATKSAMRecord read, int adjustment) {
        byte[] quals = insertSizeQualHistograms.get(read.getReadGroup().getId());
        if ( quals == null ) {
            logger.warn("No insert size histogram for RGID: " + read.getReadGroup().getId()+" returning an uninformative probability");
            return QualityUtils.qualToErrorProb((byte) 4);
        }
        int size = Math.max(1,Math.min(quals.length-1,Math.abs(read.getInferredInsertSize())-adjustment));
        return Math.log10(QualityUtils.qualToErrorProb(quals[size]));
    }

    enum PostulateClass {
        LOSS, NONE
    }

    class IntronLossPostulate {

        Pair<GenomeLoc,GenomeLoc> exonLocs;
        Pair<Integer,Integer> exonNums;
        PostulateClass clazz;
        Set<Pair<GATKSAMRecord,GATKSAMRecord>> supportingPairs = new HashSet<Pair<GATKSAMRecord,GATKSAMRecord>>(1000);

        public IntronLossPostulate(List<GenomeLoc> exons, GenomeLoc read, GenomeLoc mate) {
            if ( read == null || mate == null ) {
                clazz = PostulateClass.NONE;
                exonLocs = null;
                exonNums = null;
                return;
            }
            int eReadNum = -1, eMateNum = -1;
            int eNum = 0;
            for ( GenomeLoc e : exons ) {
                if ( e.overlapsP(read) ) {
                    eReadNum = eNum;
                }

                if ( e.overlapsP(mate) ) {
                    eMateNum = eNum;
                }
                eNum++;
            }

            if ( eReadNum != eMateNum && eMateNum > -1 && eReadNum > -1 ) {
                clazz = PostulateClass.LOSS;
                exonLocs = new Pair<GenomeLoc, GenomeLoc>(exons.get(eReadNum),exons.get(eMateNum));
                exonNums = new Pair<Integer, Integer>(eReadNum,eMateNum);
            } else {
                clazz = PostulateClass.NONE;
                exonLocs = null;
                exonNums = null;
            }
        }

        public String testConstructor(List<GenomeLoc> exons, GenomeLoc read, GenomeLoc mate) {
            StringBuffer buf = new StringBuffer();
            int eReadNum = -1, eMateNum = -1;
            int eNum = 0;
            buf.append(read);
            buf.append(" :: ");
            buf.append(mate);
            for ( GenomeLoc e : exons ) {
                buf.append("\t");
                buf.append(e);
                if ( e.overlapsP(read) ) {
                    eReadNum = eNum;
                    buf.append("\tREAD_OVERLAP");
                }

                if ( e.overlapsP(mate) ) {
                    eMateNum = eNum;
                    buf.append("\tMATE_OVERLAP");
                }
            }

            return buf.toString();
        }

        public boolean equals(Object other) {
            if ( ! ( other instanceof IntronLossPostulate ) ) {
                return false;
            }
            IntronLossPostulate otherPostulate = (IntronLossPostulate) other;
            return ( (exonNums.first == otherPostulate.exonNums.first) &&
                    exonNums.second == otherPostulate.exonNums.second ) ||
                    ( ( exonNums.second == otherPostulate.exonNums.first) &&
                    ( exonNums.first == otherPostulate.exonNums.second) );
        }

        public int hashCode() {
            return exonNums.first.hashCode() ^ exonNums.second.hashCode();
        }

        public String toString() {
            if ( clazz == PostulateClass.NONE ) {
                return "NONE";
            }

            return String.format("EX%d-EX%d",exonNums.first,exonNums.second);
        }

        public int getEventSize() {
            return exonLocs.first.minDistance(exonLocs.second);
        }

        public boolean isValid() {
            // todo -- this needs to somehow use the insert size distribution
            return getEventSize() > 350;
        }

        public byte[] getSequence(IndexedFastaSequenceFile ref) {
            // BIG TODO -- SPLICE DONOR/ACCEPTOR SITES AND STRAND OF TRANSCRIPT
            StringBuffer buf = new StringBuffer();
            buf.append(new String(ref.getSubsequenceAt(exonLocs.first.getContig(),exonLocs.first.getStart(),exonLocs.first.getStop()).getBases()));
            buf.append(new String(ref.getSubsequenceAt(exonLocs.second.getContig(),exonLocs.second.getStart(),exonLocs.second.getStop()).getBases()));
            return buf.toString().getBytes();
        }

        public boolean isRelevant(GenomeLoc loc1, GenomeLoc loc2) {
            if ( loc1 != null && isRelevant(loc1) ) {
                // loc 2 can't be way too far away
                return loc2 == null || loc2.overlapsP(exonLocs.first) || loc2.overlapsP(exonLocs.second);
            } else if ( loc2 != null && isRelevant(loc2) ) {
                return loc1 == null || loc1.overlapsP(exonLocs.first) || loc1.overlapsP(exonLocs.second);
            }

            return false;
        }

        @Requires("! loc.overlapsP(exonLocs.second)")
        public boolean isRelevant(GenomeLoc loc) {
            // basically if the read hangs partially in the intron
            boolean r1 =  loc.overlapsP(exonLocs.first) && loc.getStop() > exonLocs.first.getStop();
            boolean r2 = loc.overlapsP(exonLocs.second) && loc.getStart() < exonLocs.second.getStart();
            return r1 || r2;
        }
    }
}
