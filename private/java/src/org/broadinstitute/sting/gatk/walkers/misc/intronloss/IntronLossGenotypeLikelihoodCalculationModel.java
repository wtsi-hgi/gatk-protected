package org.broadinstitute.sting.gatk.walkers.misc.intronloss;

import com.google.java.contract.Requires;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordUtil;
import net.sf.samtools.SAMUtils;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.features.refseq.RefSeqFeature;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypePriors;
import org.broadinstitute.sting.gatk.walkers.genotyper.MultiallelicGenotypeLikelihoods;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: loaner
 * Date: 8/16/11
 * Time: 10:01 AM
 * To change this template use File | Settings | File Templates.
 */
public class IntronLossGenotypeLikelihoodCalculationModel extends GenotypeLikelihoodsCalculationModel {

    public static final int PLOIDY = 2;

    // todo -- deletion is with respect to a specific transcript, so intron loss can remove an exon
    //     ... so we could be looking 2 or 3 exons away, or even in the middle of an intron if there is
    //     ... an alternate transcript
    private boolean isInExon;
    private int exonNum;
    private GenomeLoc currentExon;
    private boolean hasNextExon;
    private boolean hasPreviousExon;
    private GenomeLoc nextExon;
    private GenomeLoc previousExon;
    private RodBinding<RefSeqFeature> refSeqBinding;
    private GenomeLocParser genomeLocParser;
    private Map<String,byte[]> insertSizeQualHistograms;

    public IntronLossGenotypeLikelihoodCalculationModel(RodBinding<RefSeqFeature> refSeqBinding,
                                                        GenomeLocParser parser, Map<String,byte[]> histograms) {
        super(null,null);
        this.refSeqBinding = refSeqBinding;
        this.genomeLocParser = parser;
        this.insertSizeQualHistograms = histograms;
    }

    public void resetRefSeqData(RefSeqFeature refSeqFeature, GenomeLoc loc) {
        List<GenomeLoc> exons = refSeqFeature.getExons();
        GenomeLoc closestPrev = null, closestNext = null;
        int prevDist = Integer.MIN_VALUE;
        int eNum = 0;
        for ( GenomeLoc exon : exons ) {
            if ( exon.overlapsP(loc) ) {
                isInExon = true;
                exonNum = eNum;
                currentExon = exon;
                hasNextExon = exonNum < exons.size()-1;
                hasPreviousExon = exonNum > 0;
                break;
            }
            int thisDist = loc.distance(exon) * (exon.isBefore(loc) ? -1 : 1);
            if ( thisDist > prevDist ) {
                if ( thisDist < 0 ) {
                    prevDist = thisDist;
                    closestPrev = exon;
                } else {
                    closestNext = exon;
                    prevDist = Integer.MAX_VALUE;
                }
            }
            eNum++;
        }

        if ( isInExon && hasNextExon ) {
            nextExon = exons.get(exonNum+1);
        } else {
            hasNextExon = closestNext != null;
            nextExon = closestNext;
        }

        if ( isInExon && hasPreviousExon ) {
            previousExon = exons.get(exonNum-1);
        } else {
            hasPreviousExon = closestPrev != null;
            previousExon = closestPrev;
        }
    }

    public Allele getLikelihoods(RefMetaDataTracker tracker,
                                          ReferenceContext ref,
                                          Map<String, AlignmentContext> contexts,
                                          AlignmentContextUtils.ReadOrientation contextType,
                                          GenotypePriors priors,
                                          Map<String, MultiallelicGenotypeLikelihoods> GLs,
                                          Allele alternateAlleleToUse, boolean useBAQedPileup) {
        if ( tracker.hasValues(refSeqBinding) ) {
            resetRefSeqData(tracker.getFirstValue(refSeqBinding),ref.getLocus());
        } else {
            throw new ReviewedStingException("Location "+ref.getLocus().toString()+" is not within a gene");
        }

        for ( Map.Entry<String,AlignmentContext> sample : contexts.entrySet() ) {
            MultiallelicGenotypeLikelihoods samGLs = GLs.get(sample.getKey());
            double[] samLikelihoods = samGLs.getLikelihoods();
            ReadBackedPileup pileup = AlignmentContextUtils.stratify(sample.getValue(), contextType).getBasePileup();
            updateLikelihoods(samLikelihoods, pileup);
            GLs.put(sample.getKey(), new MultiallelicGenotypeLikelihoods(sample.getKey(),
                    samGLs.getAlleles(),  samLikelihoods, pileup.size()));
        }

        return Allele.NO_CALL;
    }

    public void updateLikelihoods(double[] likArray, ReadBackedPileup stratifiedPileup) {
        for (PileupElement e : stratifiedPileup) {
            ILReadClass ilClass = getReadClass(e.getRead());
            double probIntDel = generativeIntronDeletionProbability(e.getRead(),ilClass);
            double probNoDel = generativeUnmodifiedGenomeProbability(e.getRead(),ilClass);
            // we may want to generalize ploidy in the future
            /* old
             * // hom ref:
             * likArray[0] += Math.log10(probNoDel);
             * // het:
             * likArray[1] += Math.log10(probIntDel/2 + probNoDel/2);
             * // hom var:
             * likArray[2] += Math.log10(probNoDel);
             */
            int ploidy = PLOIDY;
            for ( int nAlt = 0; nAlt <= ploidy; nAlt++ ) {
                likArray[nAlt] = Math.log10((nAlt*probIntDel)/ploidy + ((ploidy-nAlt)*probNoDel)/ploidy);
            }
        }
    }

    /**
     * The probability that an intronic deletion would have generated the given read [pair]. Utilizes
     * mate information in the SAMRecord (i.e. "MA:" the alignment position of the mate)
     * @param read - the read
     * @return probability that intron deletion would have generated the given read (pair)
     */
    private double generativeIntronDeletionProbability(SAMRecord read, ILReadClass readClass) {
        switch (readClass) {
            case IN_SAME_EXON: // just that both reads are properly aligned
                return readProperlyAligned(read)*mateProperlyAligned(read);
            case EXON_INTRON: // that one or both reads are IMPROPERLY aligned
                return readImproperlyAligned(read) + mateImproperlyAligned(read);
            case EXON_EXON: // both reads aligned, P[insert size] given the deletion event
                return readProperlyAligned(read)*mateProperlyAligned(read)*insertSizeProbability(read,getEventSize(read));
            case INTRON_EXON: // that one or both reads are improperly aligned
                return readImproperlyAligned(read) + mateImproperlyAligned(read);
            case INTRON_INTRON: // that both reads are improperly aligned
                return readImproperlyAligned(read)*mateImproperlyAligned(read);
        }

        throw new ReviewedStingException("This read has a read class which is not enumerated: "+readClass.name());
    }

    /**
     *
     * @param read
     * @return
     */
    private double generativeUnmodifiedGenomeProbability(SAMRecord read, ILReadClass readClass) {
        switch ( readClass ) {
            case IN_SAME_EXON: // both reads proper
                return readProperlyAligned(read)*mateProperlyAligned(read);
            case EXON_INTRON: // both reads proper
                return readProperlyAligned(read)*mateProperlyAligned(read);
            case EXON_EXON: // either read improperly aligned OR insert size
                double pImp = readImproperlyAligned(read)+mateImproperlyAligned(read);
                double pInsert = (1-pImp)*insertSizeProbability(read,0);
                return pImp+pInsert;
            case INTRON_EXON: // both reads proper
                return readProperlyAligned(read)*mateProperlyAligned(read);
            case INTRON_INTRON:
                return readProperlyAligned(read)*mateProperlyAligned(read);
        }

        throw new ReviewedStingException("This read has a read class which is not enumerated: "+readClass.name());
    }

    // todo -- there is a proper way to do these four calculations

    private double readProperlyAligned(SAMRecord read) {
        return QualityUtils.qualToProb(read.getMappingQuality());
    }

    private double mateProperlyAligned(SAMRecord read) {
        return QualityUtils.qualToProb(Byte.parseByte( (String) read.getAttribute("MQ")));
    }

    private double readImproperlyAligned(SAMRecord read) {
        return QualityUtils.qualToErrorProb((byte)read.getMappingQuality());
    }

    private double mateImproperlyAligned(SAMRecord read) {
        return QualityUtils.qualToErrorProb(Byte.parseByte( (String) read.getAttribute("MQ")));
    }

    // todo -- end todo

    private double insertSizeProbability(SAMRecord read, int adjustment) {
        byte[] quals = insertSizeQualHistograms.get(read.getReadGroup().getId());
        int size = Math.min(quals.length-1,Math.abs(read.getInferredInsertSize())-adjustment);
        return QualityUtils.qualToErrorProb(quals[size]);
    }

    @Requires("isInExon==true")
    private int getEventSize(SAMRecord read) {
        if ( genomeLocParser.createGenomeLoc(read).overlapsP(nextExon) )
            return nextExon.distance(currentExon);
        return currentExon.distance(previousExon);
    }

    private ILReadClass getReadClass(SAMRecord read) {
        GenomeLoc readLoc = genomeLocParser.createGenomeLoc(read);
        int iSize = read.getInferredInsertSize();
        GenomeLoc mate = genomeLocParser.createGenomeLoc(readLoc.getContig(),readLoc.getStart()+iSize,readLoc.getStop()+iSize);
        if ( isInExon && readLoc.overlapsP(currentExon) ) {
            // either IN_SAME_EXON, EXON_INTRON, or EXON_EXON
            if ( mate.overlapsP(currentExon) ) {
                return ILReadClass.IN_SAME_EXON;
            } else if ( mate.overlapsP(nextExon) || mate.overlapsP(previousExon) ) {
                return ILReadClass.EXON_EXON;
            }

            return ILReadClass.EXON_INTRON;
        } else {
            // either INTRON_INTRON or INTRON_EXON
            if ( mate.overlapsP(nextExon) || mate.overlapsP(previousExon) ) {
                return ILReadClass.INTRON_EXON;
            }
        }

        return ILReadClass.INTRON_INTRON;
    }

    enum ILReadClass {
        IN_SAME_EXON,EXON_INTRON,EXON_EXON,INTRON_EXON,INTRON_INTRON
    }
}
