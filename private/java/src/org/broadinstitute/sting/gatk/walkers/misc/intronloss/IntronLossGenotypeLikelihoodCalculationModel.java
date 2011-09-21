package org.broadinstitute.sting.gatk.walkers.misc.intronloss;

import com.google.java.contract.Requires;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import org.apache.commons.lang.ArrayUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.genotyper.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.clipreads.ReadClipper;
import org.broadinstitute.sting.utils.codecs.refseq.RefSeqFeature;
import org.broadinstitute.sting.utils.codecs.vcf.VCF3Codec;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.collections.PrimitivePair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import sun.rmi.rmic.iiop.ContextStack;

import java.lang.instrument.IllegalClassFormatException;
import java.util.*;

/**
 *
 */
public class IntronLossGenotypeLikelihoodCalculationModel {

    public static final int PLOIDY = 2;
    public static final double MAX_LIKELIHOOD = Math.log10(1.0-Math.pow(10.0,-12));
    public static final double MIN_LIKELIHOOD = Double.MIN_VALUE;
    public static final int ALIGNMENT_QUAL_CAP = 500;

    // todo -- deletion is with respect to a specific transcript, so intron loss can remove an exon
    //     ... so we could be looking 2 or 3 exons away, or even in the middle of an intron if there is
    //     ... an alternate transcript
    private RodBinding<RefSeqFeature> refSeqBinding;
    private GenomeLocParser genomeLocParser;
    private Map<String,byte[]> insertSizeQualHistograms;
    private Set<String> samples;
    private Logger logger;

    public IntronLossGenotypeLikelihoodCalculationModel(Logger logger) {
        this.logger = logger;
    }

    public void initialize(RodBinding<RefSeqFeature> refSeqBinding,
                           GenomeLocParser parser, Map<String,byte[]> histograms) {
        this.refSeqBinding = refSeqBinding;
        this.genomeLocParser = parser;
        this.insertSizeQualHistograms = histograms;
    }

    public void setSamples(Set<String> s) {
        samples =s;
    }

    public GenomeLoc getGenomeLoc(SAMRecord read) {
        //String debugStr = read == null ? "null" : String.format("%s %s:%d-%d",read.getReadName(),read.getReferenceName(),ReadUtils.getRefCoordSoftUnclippedStart(read),ReadUtils.getRefCoordSoftUnclippedEnd(read));
        //logger.debug(debugStr);
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
        for ( Pair<SAMRecord,SAMRecord> pair : readBin ) {
            IntronLossPostulate p = new IntronLossPostulate(geneFeature.getExons(),getGenomeLoc(pair.first),getGenomeLoc(pair.second));
            //logger.debug(String.format("%s %d %s %s", genomeLocParser.createGenomeLoc(r), r.getInferredInsertSize(), createMateLoc(r), Utils.join(",",geneFeature.getExons())));
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
                // update the likelihoods for the postulate itself
                //Map<String,double[]> oldLikelihoods = createPostulateLikelihoods(samples,postulate,postulates.get("NONE"),geneFeature,refFile);
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


        for ( Pair<SAMRecord,SAMRecord> supportingPair : myPostulate.supportingPairs ) {
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

        /*
        * Note: Early abort if we are more likely to be ref from SUPPORTING READS
        boolean probRef = true;
        for ( double[] lik : samLikelihoods.values() ) {
            probRef &= MathUtils.compareDoubles(MathUtils.arrayMax(lik),lik[0],1e-5) == 0;
        }
        if ( probRef ) {
            return samLikelihoods;
        }
         */

        if ( nonePostulate == null ) {
            return samLikelihoods;
        }

        for ( Pair<SAMRecord,SAMRecord> nonePair : nonePostulate.supportingPairs ) {
            double lossScore = 0.0;
            double noLossScore = 0.0;
            if ( myPostulate.isRelevant(getGenomeLoc(nonePair.first),getGenomeLoc(nonePair.second)) ) {
                double realignmentScore = scoreRealignment(nonePair,junctionSeq);
                double initialScore = scoreAlignment(nonePair);
                lossScore += realignmentScore;
                noLossScore += initialScore;
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

    private double scoreAlignment(Pair<SAMRecord,SAMRecord> initiallyAlignedPair) {
        // aligned reads are equipped already with necessary metrics
        double logScore = 0.0;
        for ( SAMRecord read : Arrays.asList(initiallyAlignedPair.first,initiallyAlignedPair.second) ) {
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

    private double scoreRealignment(Pair<SAMRecord,SAMRecord> inPair, byte[] sequence) {
        double logScore = 0.0;
        for ( SAMRecord rec : Arrays.asList(inPair.first,inPair.second) ) {
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
                phredBasesAreExonBases += 10*(2-Math.pow(2,1-elem.getLength()));
            } else if (elem.getOperator().equals(CigarOperator.I)) {
                // update the b
                b += elem.getLength();
                offset -= elem.getLength(); // compensate, since pos will move up by elem.getLength()
                phredBasesAreExonBases += 10+20*(2-Math.pow(2.0,-1.0*elem.getLength()));
            } else if ( elem.getOperator().equals(CigarOperator.S) ) {
                // just skip the bases and update b; soft clipping occurs at edges, so need to subtract from offset
                b += elem.getLength();
                offset -= elem.getLength();
                phredBasesAreExonBases += 5*elem.getLength();
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

    public Map<String,double[]> createPostulateLikelihoods(Set<String> samples, IntronLossPostulate myPostulate, IntronLossPostulate nonePostulate, RefSeqFeature geneFeature, IndexedFastaSequenceFile refFile) {
        Map<String,double[]> samLikelihoods = new HashMap<String,double[]>(samples.size());
        for ( String s : samples ) {
            samLikelihoods.put(s,new double[1+PLOIDY]);
        }

        // step one: update the likelihoods based on pairs supporting the postulate
        for ( Pair<SAMRecord,SAMRecord> supportingPair : myPostulate.supportingPairs ) {
            double deletionProbabilityLog = generativeIntronDeletionProbabilityLog(supportingPair, ILReadClass.EXON_EXON,myPostulate.getEventSize());
            double noDeletionProbabilityLog = generativeUnmodifiedGenomeProbabilityLog(supportingPair, ILReadClass.EXON_EXON,myPostulate.getEventSize());
            for ( SAMRecord read : Arrays.asList(supportingPair.first,supportingPair.second) ) {
                GenomeLoc readLoc = getGenomeLoc(read);
                if ( readProximalToEvent(read,myPostulate) && (readLoc.getStart() < myPostulate.exonLocs.second.getStart() || readLoc.getStop() > myPostulate.exonLocs.first.getStop())) {
                    // read overlaps either the intron or spans the event, so this is a supporting read
                    deletionProbabilityLog += overlapIntronDeletionProbabilityLog(read,myPostulate,geneFeature,refFile);
                    noDeletionProbabilityLog += overlapIntronNoDeletionProbabilityLog(read,myPostulate,geneFeature,refFile);
                }
            }
            for ( int ac = 0; ac <= PLOIDY; ac ++ ) {
                double logAc = Math.log10(ac);
                double logAnMAc = Math.log10(PLOIDY-ac);
                double logPloidy = Math.log10(PLOIDY); // todo -- should be global constant
                // todo -- numerically unstable: delta can be negative infinity
                double delta = MathUtils.log10sumLog10(new double[]{logAc + deletionProbabilityLog - logPloidy, logAnMAc + noDeletionProbabilityLog - logPloidy});
                String sam = supportingPair.first.getReadGroup().getSample();
                samLikelihoods.get(sam)[ac] += delta;
            }
        }

        // step two: update likelihoods based on read pairs possibly not supporting the postulate
        if ( nonePostulate == null )
            return samLikelihoods;

        for ( Pair<SAMRecord,SAMRecord> otherPair : nonePostulate.supportingPairs ) {
            // set the delta
            double deletionProbabilityLog = 0.0;
            double noDeletionProbabilityLog = 0.0;
            // first, check to see if the pair is an intron -> exon jump
            if ( otherPair.second != null && myPostulate.isIntronExon(getGenomeLoc(otherPair.first),getGenomeLoc(otherPair.second))) {
                deletionProbabilityLog += generativeIntronDeletionProbabilityLog(otherPair, ILReadClass.EXON_INTRON, myPostulate.getEventSize());
                noDeletionProbabilityLog += generativeIntronDeletionProbabilityLog(otherPair, ILReadClass.EXON_INTRON, myPostulate.getEventSize());
            }
            // next identify any reads spanning exon/intron or exon/exon boundary
            for ( SAMRecord read : Arrays.asList(otherPair.first,otherPair.second) ) {
                if ( read == null ) {
                    continue;
                }
                GenomeLoc readLoc = getGenomeLoc(read);
                // check if it extends beyind an exon of the postulate
                if ( readProximalToEvent(read,myPostulate) && (readLoc.getStart() < myPostulate.exonLocs.second.getStart() || readLoc.getStop() > myPostulate.exonLocs.first.getStop())) {
                    // read overlaps either the intron or spans the event, so this is a supporting read
                    deletionProbabilityLog += overlapIntronDeletionProbabilityLog(read, myPostulate, geneFeature, refFile);
                    noDeletionProbabilityLog += overlapIntronNoDeletionProbabilityLog(read, myPostulate, geneFeature, refFile);
                }
            }

            if ( deletionProbabilityLog < -1e-2 && noDeletionProbabilityLog < -1e-2 ) {
                // probabilities have been updated
                for ( int ac = 0; ac <= PLOIDY ; ac ++ ) {
                    double logAc = Math.log10(ac);
                    double logAnMAc = Math.log10(PLOIDY-ac);
                    double logPloidy = Math.log10(PLOIDY); // todo -- should be global constant
                    // todo -- numerically unstable: delta can be negative infinity
                    double delta = MathUtils.log10sumLog10(new double[]{logAc + deletionProbabilityLog - logPloidy, logAnMAc + noDeletionProbabilityLog - logPloidy});
                    String sam = otherPair.first.getReadGroup().getSample();
                    samLikelihoods.get(sam)[ac] += delta;
                }
            }
        }

        return samLikelihoods;
    }

    private double overlapIntronDeletionProbabilityLog(SAMRecord read, IntronLossPostulate postulate, RefSeqFeature gene, IndexedFastaSequenceFile ref) {
        // want the bases and qualities of the positions extending past the intron, and the bases and qualities of both the reference, and into the postulated exon
        GenomeLoc uClipReadLoc = getGenomeLoc(read);
        GenomeLoc overlExonLoc = uClipReadLoc.overlapsP(postulate.exonLocs.first) ? postulate.exonLocs.first : postulate.exonLocs.second;
        GenomeLoc nextExonLoc =   uClipReadLoc.overlapsP(postulate.exonLocs.first) ? postulate.exonLocs.second : postulate.exonLocs.first;
        GenomeLoc overlapLoc = uClipReadLoc.intersect(overlExonLoc);
        boolean lastBases = uClipReadLoc.getStop() > overlExonLoc.getStop();
        int numBases = lastBases ? uClipReadLoc.getStop() - overlExonLoc.getStop() : overlExonLoc.getStart()-uClipReadLoc.getStart();
        if ( numBases >= read.getReadBases().length) {
            return 1.0;
        }
        //logger.debug(uClipReadLoc.toString()+" "+overlExonLoc.toString() + " "+overlapLoc.toString()+" "+Boolean.toString(lastBases));
        //logger.debug(String.format("%d-%d=%d",uClipReadLoc.getStop(),overlExonLoc.getStop(),numBases));
        //byte[] intronBases = ref.getSubsequenceAt(overlapLoc.getContig(),overlapLoc.getStop()+1,overlapLoc.getStop()+1+numBases).getBases();
        StringBuffer exonSeq = new StringBuffer();
        if ( lastBases ) {
            exonSeq.append(new String(ref.getSubsequenceAt(overlExonLoc.getContig(),overlapLoc.getStart()-3,overlapLoc.getStop()).getBases()));
            exonSeq.append(new String(ref.getSubsequenceAt(nextExonLoc.getContig(),nextExonLoc.getStart(),nextExonLoc.getStart()+numBases+3).getBases()));
        } else {
            exonSeq.append(new String(ref.getSubsequenceAt(nextExonLoc.getContig(),nextExonLoc.getStop()-numBases-3,nextExonLoc.getStop()).getBases()));
            exonSeq.append(new String(ref.getSubsequenceAt(overlapLoc.getContig(),overlapLoc.getStart(),overlapLoc.getStop()+3).getBases()));
        }
        byte[] exonBases = exonSeq.toString().getBytes(); // lower-case bases are crap
        // see if the overlapping bases match the exon, and with what probability they do not
        // todo -- do we need to locally shift at all? -- yup!
        SWPairwiseAlignment exonAlignment = new SWPairwiseAlignment(exonBases,read.getReadBases());
        int offset = exonAlignment.getAlignmentStart2wrt1();
        double logProbBasesAreExonBases = 0.0;
        //logger.debug(new String(read.getReadBases()));
        //logger.debug(new String(exonBases));
        //logger.debug(new String(read.getReadBases()).substring(offset,offset+numBases));
        // todo -- CIGAR string!
        int b = 0;
        for ( CigarElement elem : exonAlignment.getCigar().getCigarElements() ) {
            if ( elem.getOperator().equals(CigarOperator.M) ) {
                for ( int ct = elem.getLength(); ct > 0; ct-- ) {
                    int pos = b + offset;
                    if (BaseUtils.basesAreEqual(read.getReadBases()[b],exonBases[pos])) {
                        // the probability that the base IS what it says it is
                        logProbBasesAreExonBases += Math.log10(QualityUtils.qualToProb(read.getBaseQualities()[b]));
                    } else {
                        // the probability that the base ISN'T what it says it is
                        logProbBasesAreExonBases += Math.log10(QualityUtils.qualToErrorProb(read.getBaseQualities()[b]))-3; // epsilon over 3
                    }
                    b++;
                }
            } else if ( elem.getOperator().equals(CigarOperator.D) ) {
                // just increment offset
                offset += elem.getLength();
            } else if (elem.getOperator().equals(CigarOperator.I)) {
                // update the b
                b += elem.getLength();
                offset -= elem.getLength(); // compensate, since pos will move up by elem.getLength()
            } else if ( elem.getOperator().equals(CigarOperator.S) || elem.getOperator().equals(CigarOperator.H) ) {
                // just skip the bases and update b; soft clipping occurs at edges, so need to subtract from offset
                b += elem.getLength();
                offset -= elem.getLength();
            } else {
                throw new StingException("Unsupported operator: "+elem.getOperator().toString());
            }
        }

        return logProbBasesAreExonBases; // todo -- see below ("make this better") ??
    }

    // todo -- code duplication. Merge with above.
    private double overlapIntronNoDeletionProbabilityLog(SAMRecord read, IntronLossPostulate postulate, RefSeqFeature gene, IndexedFastaSequenceFile ref) {
        // want the bases and qualities of the positions extending past the intron, and the bases and qualities of both the reference, and into the postulated exon
        GenomeLoc uClipReadLoc = getGenomeLoc(read);
        GenomeLoc overlExonLoc = uClipReadLoc.overlapsP(postulate.exonLocs.first) ? postulate.exonLocs.first : postulate.exonLocs.second;
        GenomeLoc overlapLoc = uClipReadLoc.intersect(overlExonLoc);
        boolean lastBases = uClipReadLoc.getStop() > overlExonLoc.getStop();
        int numBases = lastBases ? uClipReadLoc.getStop() - overlExonLoc.getStop() : overlExonLoc.getStart()-uClipReadLoc.getStart();
        if ( numBases >=  read.getReadBases().length ) {
            return 1.0;
        }
        byte[] intronBases = ref.getSubsequenceAt(overlapLoc.getContig(),overlapLoc.getStop()+1,overlapLoc.getStop()+numBases+1).getBases();
        //byte[] exonBases = ref.getSubsequenceAt(overlExonLoc.getContig(),overlExonLoc.getStart(),overlExonLoc.getStart()+numBases).toString().toUpperCase().getBytes(); // lower-case bases are crap
        // see if the overlapping bases match the exon, and with what probability they do not
        // todo -- do we need to locally shift at all?
        double logProbBasesAreExonBases = 0.0;
        int offset = lastBases ? read.getReadBases().length-numBases  : 0;
        for ( int b = 0; b < numBases; b++ ) {
            int pos = b+offset;
            if (BaseUtils.basesAreEqual(read.getReadBases()[pos],intronBases[b])) {
                // the probability that the base IS what it says it is
                logProbBasesAreExonBases += Math.log10(QualityUtils.qualToProb(read.getBaseQualities()[pos]));
            } else {
                // the probability that the base ISN'T what it says it is
                logProbBasesAreExonBases += Math.log10(QualityUtils.qualToErrorProb(read.getBaseQualities()[pos]))-3; // epsilon over 3
            }
        }

        return logProbBasesAreExonBases; // max of Q80 todo -- make this better ??
    }

    public boolean readProximalToEvent(SAMRecord read, IntronLossPostulate postulate) {
        GenomeLoc rLoc = getGenomeLoc(read);
        return  rLoc.overlapsP(postulate.exonLocs.first) && ! postulate.exonLocs.first.containsP(rLoc)|| rLoc.overlapsP(postulate.exonLocs.second) && ! postulate.exonLocs.second.containsP(rLoc);
    }

    private double generativeIntronDeletionProbabilityLog(Pair<SAMRecord,SAMRecord> readPair, ILReadClass readClass, int eventSize) {
        switch (readClass) {
            case EXON_INTRON: // that one or both reads are IMPROPERLY aligned
                return Math.log10( Math.pow(10.0,readImproperlyAligned(readPair.first)) + Math.pow(10.0, readImproperlyAligned(readPair.second)) );
            case EXON_EXON: // both reads aligned, P[insert size] given the deletion event
                double g =  readProperlyAligned(readPair.first) + readProperlyAligned(readPair.second) + insertSizeProbability(readPair.first,eventSize);
                return g;
            case INTRON_EXON: // that one or both reads are improperly aligned
                return Math.log10( Math.pow(10.0,readImproperlyAligned(readPair.first)) + Math.pow(10.0,readProperlyAligned(readPair.second)));
        }

        throw new ReviewedStingException("This read has a read class which is not enumerated: "+readClass.name());
    }

    private double generativeUnmodifiedGenomeProbabilityLog(Pair<SAMRecord,SAMRecord> pair, ILReadClass readClass, int eventSize) {
        switch ( readClass ) {
            case EXON_INTRON:
                return readProperlyAligned(pair.first)+readProperlyAligned(pair.second);
            case EXON_EXON: // either read improperly aligned OR insert size
                double pImp = Math.log10(Math.pow(10.0,readImproperlyAligned(pair.first))+Math.pow(10.0,readImproperlyAligned(pair.second)));
                double pInsert = Math.log10((1-pImp)) + insertSizeProbability(pair.first,0);
                return Math.log10(Math.pow(10.0,pImp)+Math.pow(10.0,pInsert));
            case INTRON_EXON: // both reads proper
                return readProperlyAligned(pair.first)+readProperlyAligned(pair.second);
        }

        throw new ReviewedStingException("This read has a read class which is not enumerated: "+readClass.name());

    }
    // todo -- there is a proper way to do these four calculations

    private double readProperlyAligned(SAMRecord read) {
        return Math.log10(QualityUtils.qualToProb(read.getMappingQuality()));
    }

    private double readImproperlyAligned(SAMRecord read) {
        return Math.log10(QualityUtils.qualToErrorProb(((Integer) read.getMappingQuality()).byteValue()));
    }

    // todo -- end todo

    private double insertSizeProbability(SAMRecord read, int adjustment) {
        byte[] quals = insertSizeQualHistograms.get(read.getReadGroup().getId());
        if ( quals == null ) {
            logger.warn("No insert size histogram for RGID: " + read.getReadGroup().getId()+" returning an uninformative probability");
            return QualityUtils.qualToErrorProb((byte) 4);
        }
        int size = Math.max(1,Math.min(quals.length-1,Math.abs(read.getInferredInsertSize())-adjustment));
        return Math.log10(QualityUtils.qualToErrorProb(quals[size]));
    }

    enum ILReadClass {
        IN_SAME_EXON,EXON_INTRON,EXON_EXON,INTRON_EXON,INTRON_INTRON
    }

    enum PostulateClass {
        LOSS, NONE
    }

    class IntronLossPostulate {

        Pair<GenomeLoc,GenomeLoc> exonLocs;
        Pair<Integer,Integer> exonNums;
        PostulateClass clazz;
        Set<Pair<SAMRecord,SAMRecord>> supportingPairs = new HashSet<Pair<SAMRecord,SAMRecord>>(1000);

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

        public boolean isIntronExon(GenomeLoc first, GenomeLoc second) {
            return ( first.overlapsP(exonLocs.first) && ! second.overlapsP(exonLocs.first) && ! second.overlapsP(exonLocs.second) ) ||
                     ( first.overlapsP(exonLocs.second) && ! second.overlapsP(exonLocs.second) && ! second.overlapsP(exonLocs.first) ) ||
                    ( second.overlapsP(exonLocs.first) && ! first.overlapsP(exonLocs.first) && ! first.overlapsP(exonLocs.second) ) ||
                    ( second.overlapsP(exonLocs.second) && ! first.overlapsP(exonLocs.first) && ! first.overlapsP(exonLocs.second));
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
            return loc.overlapsP(exonLocs.first) && loc.getStop() > exonLocs.first.getStop();
        }
    }
}
