/*
 * Copyright (c) 2011 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.walkers.indels.ConstrainedMateFixingManager;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SWPairwiseAlignment;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.InferredGeneticContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;

public class GenotypingEngine {

    // Smith-Waterman parameters copied from IndelRealigner
    private static final double SW_MATCH = 23.0;      // 1.0;
    private static final double SW_MISMATCH = -8.0;  //-1.0/3.0;
    private static final double SW_GAP = -14.0;       //-1.0-1.0/3.0;
    private static final double SW_GAP_EXTEND = -1.5; //-1.0/.0;

    public List<VariantContext> alignAndGenotype( final Pair<Haplotype, Haplotype> bestTwoHaplotypes, final byte[] ref, final GenomeLoc loc ) {
        final SWPairwiseAlignment swConsensus1 = new SWPairwiseAlignment( ref, bestTwoHaplotypes.first.bases, SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );
        final SWPairwiseAlignment swConsensus2 = new SWPairwiseAlignment( ref, bestTwoHaplotypes.second.bases, SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );

        System.out.println( bestTwoHaplotypes.first.toString() );
        System.out.println( "Cigar = " + swConsensus1.getCigar() );
        final List<VariantContext> vcs1 = generateVCsFromAlignment( swConsensus1, ref, bestTwoHaplotypes.first.bases, loc );

        System.out.println( bestTwoHaplotypes.second.toString() );
        System.out.println( "Cigar = " + swConsensus2.getCigar() );
        final List<VariantContext> vcs2 = generateVCsFromAlignment( swConsensus2, ref, bestTwoHaplotypes.second.bases, loc );

        return genotype( vcs1, vcs2 );
    }

    private static List<VariantContext> generateVCsFromAlignment( final SWPairwiseAlignment swConsensus, final byte[] ref, final byte[] read, final GenomeLoc loc ) {
        final ArrayList<VariantContext> vcs = new ArrayList<VariantContext>();

        int refPos = swConsensus.getAlignmentStart2wrt1();
        int readPos = 0;
        final int lookAhead = 5;

        for( final CigarElement ce : swConsensus.getCigar().getCigarElements() ) {
            final int elementLength = ce.getLength();
            switch( ce.getOperator() ) {
                case I:
                {
                    byte[] insertionBases = Arrays.copyOfRange( read, readPos, readPos + elementLength);
                    boolean allN = true;
                    for( byte b : insertionBases ) {
                        if( b != (byte) 'N' ) {
                            allN = false;
                            break;
                        }
                    }
                    if( !allN ) {
                        ArrayList<Allele> alleles = new ArrayList<Allele>();
                        alleles.add( Allele.create(Allele.NULL_ALLELE_STRING, true));
                        alleles.add( Allele.create(insertionBases, false));
                        System.out.println("> Insertion: " + alleles);
                        vcs.add(new VariantContext("HaplotypeCaller", loc.getContig(), loc.getStart() + refPos - 1, loc.getStart() + refPos - 1, alleles, VariantContext.NO_GENOTYPES, InferredGeneticContext.NO_NEG_LOG_10PERROR, null, null, ref[refPos-1]));
                    }
                    readPos += elementLength;
                    break;
                }
                case S:
                {
                    readPos += elementLength;
                    refPos += elementLength;
                    break;
                }
                case D:
                {
                    byte[] deletionBases = Arrays.copyOfRange( ref, refPos, refPos + elementLength);
                    ArrayList<Allele> alleles = new ArrayList<Allele>();
                    alleles.add( Allele.create(deletionBases, true) );
                    alleles.add( Allele.create(Allele.NULL_ALLELE_STRING, false) );
                    System.out.println( "> Deletion: " + alleles);
                    vcs.add( new VariantContext("HaplotypeCaller", loc.getContig(), loc.getStart() + refPos - 1, loc.getStart() + refPos + elementLength - 1, alleles, VariantContext.NO_GENOTYPES, InferredGeneticContext.NO_NEG_LOG_10PERROR, null, null, ref[refPos-1]) );
                    refPos += elementLength;
                    break;
                }
                case M:
                {
                    int numSinceMismatch = -1;
                    int stopOfMismatch = -1;
                    int startOfMismatch = -1;
                    int refPosStartOfMismatch = -1;
                    for( int iii = 0; iii < elementLength; iii++ ) {
                        if( ref[refPos] != read[readPos] && read[readPos] != ((byte) 'N') ) {
                            // SNP or start of possible MNP
                            if( stopOfMismatch == -1 ) {
                                startOfMismatch = readPos;
                                stopOfMismatch = readPos;
                                numSinceMismatch = 0;
                                refPosStartOfMismatch = refPos;
                            } else {
                                stopOfMismatch = readPos;
                            }
                        }

                        if( stopOfMismatch != -1) {
                            numSinceMismatch++;
                        }

                        if( numSinceMismatch > lookAhead || (iii == elementLength - 1 && stopOfMismatch != -1) ) {
                            byte[] refBases = Arrays.copyOfRange( ref, refPosStartOfMismatch, refPosStartOfMismatch + (stopOfMismatch - startOfMismatch) + 1 );
                            byte[] mismatchBases = Arrays.copyOfRange( read, startOfMismatch, stopOfMismatch + 1 );
                            ArrayList<Allele> alleles = new ArrayList<Allele>();
                            alleles.add( Allele.create( refBases, true ) );
                            alleles.add( Allele.create( mismatchBases, false ) );
                            System.out.println( "> SNP/MNP: " + alleles);
                            vcs.add( new VariantContext("HaplotypeCaller", loc.getContig(), loc.getStart() + refPosStartOfMismatch, loc.getStart() + refPosStartOfMismatch + (stopOfMismatch - startOfMismatch), alleles) );
                            numSinceMismatch = -1;
                            stopOfMismatch = -1;
                            startOfMismatch = -1;
                            refPosStartOfMismatch = -1;
                        }

                        refPos++;
                        readPos++;
                    }
                    break;
                }

                case N:
                case H:
                case P:
                default:
                    throw new ReviewedStingException( "Unsupported cigar operator: " + ce.getOperator() );
            }
        }

        if( vcs.size() == 0 ) {
            System.out.println("> Reference!");
        }

        return vcs;
    }

    private static List<VariantContext> genotype( final List<VariantContext> vcs1, final List<VariantContext> vcs2 ) {
        final ArrayList<VariantContext> vcs = new ArrayList<VariantContext>();

        final Iterator<VariantContext> vcs1Iter = vcs1.iterator();
        final Iterator<VariantContext> vcs2Iter = vcs2.iterator();

        VariantContext vc1Hold = null;
        VariantContext vc2Hold = null;

        do {
            final VariantContext vc1 = ( vc1Hold != null ? vc1Hold : (vcs1Iter.hasNext() ? vcs1Iter.next() : null) );
            final VariantContext vc2 = ( vc2Hold != null ? vc2Hold : (vcs2Iter.hasNext() ? vcs2Iter.next() : null) );

            vc1Hold = null;
            vc2Hold = null;

            if( vc1 == null && vc2 != null ) {
                ArrayList<Allele> alleles = new ArrayList<Allele>();
                alleles.addAll( vc2.getAlleles() );
                Genotype gt = new Genotype( "NA12878", alleles );
                HashMap<String,Genotype> genotypeMap = new HashMap<String,Genotype>();
                genotypeMap.put("NA12878", gt);
                vcs.add( VariantContext.modifyGenotypes( vc2, genotypeMap ) );
            } else if( vc1 != null && vc2 == null ) {
                ArrayList<Allele> alleles = new ArrayList<Allele>();
                alleles.addAll( vc1.getAlleles() );
                Genotype gt = new Genotype( "NA12878", alleles );
                HashMap<String,Genotype> genotypeMap = new HashMap<String,Genotype>();
                genotypeMap.put("NA12878", gt);
                vcs.add( VariantContext.modifyGenotypes( vc1, genotypeMap ) );
            } else if( vc1 != null ) { // && vc2 != null
                if( vc1.getStart() == vc2.getStart() ) {
                    ArrayList<Allele> alleles = new ArrayList<Allele>();
                    alleles.add( vc1.getAlternateAllele(0) );
                    alleles.add( vc2.getAlternateAllele(0) );
                    if( vc1.getAlleles().equals(vc2.getAlleles()) ) { // check if alleles match
                        Genotype gt = new Genotype( "NA12878", alleles );
                        HashMap<String,Genotype> genotypeMap = new HashMap<String,Genotype>();
                        genotypeMap.put("NA12878", gt);
                        vcs.add( VariantContext.modifyGenotypes( vc1, genotypeMap ) );
                    } else { // two alleles don't match, and don't call multialleleic records yet
                        vc2Hold = vc2;
                        ArrayList<Allele> theseAlleles = new ArrayList<Allele>();
                        theseAlleles.addAll( vc1.getAlleles() );
                        Genotype gt = new Genotype( "NA12878", theseAlleles );
                        HashMap<String,Genotype> genotypeMap = new HashMap<String,Genotype>();
                        genotypeMap.put("NA12878", gt);
                        vcs.add( VariantContext.modifyGenotypes( vc1, genotypeMap ) );
                    }
                } else if( vc1.getStart() < vc2.getStart()) {
                    vc2Hold = vc2;
                    ArrayList<Allele> alleles = new ArrayList<Allele>();
                    alleles.addAll( vc1.getAlleles() );
                    Genotype gt = new Genotype( "NA12878", alleles );
                    HashMap<String,Genotype> genotypeMap = new HashMap<String,Genotype>();
                    genotypeMap.put("NA12878", gt);
                    vcs.add( VariantContext.modifyGenotypes( vc1, genotypeMap ) );
                } else {
                    vc1Hold = vc1;
                    ArrayList<Allele> alleles = new ArrayList<Allele>();
                    alleles.addAll( vc2.getAlleles() );
                    Genotype gt = new Genotype( "NA12878", alleles );
                    HashMap<String,Genotype> genotypeMap = new HashMap<String,Genotype>();
                    genotypeMap.put("NA12878", gt);
                    vcs.add( VariantContext.modifyGenotypes( vc2, genotypeMap ) );
                }
            }


        } while ( vcs1Iter.hasNext() || vcs2Iter.hasNext() || vc1Hold != null || vc2Hold != null );

        return vcs;
    }







    //////////////////////////////////////////
    //
    // Code for debug output follows
    //
    //////////////////////////////////////////        

    public void alignAllHaplotypes( final List<Haplotype> haplotypes, final byte[] ref, final GenomeLoc loc, final StingSAMFileWriter writer, final SAMRecord exampleRead ) {

        int iii = 0;
        for( final Haplotype h : haplotypes ) {
            final SWPairwiseAlignment swConsensus = new SWPairwiseAlignment( ref, h.bases, SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );
            exampleRead.setReadName("Haplotype" + iii);
            exampleRead.setReadBases(h.bases);
            exampleRead.setAlignmentStart(loc.getStart() + swConsensus.getAlignmentStart2wrt1());
            exampleRead.setCigar(swConsensus.getCigar());
            final byte[] quals = new byte[h.bases.length];
            Arrays.fill(quals, (byte) 25);
            exampleRead.setBaseQualities(quals);
            writer.addAlignment(exampleRead);
            iii++;
        }
                
    }

    public void alignAllReads( final Pair<Haplotype,Haplotype> bestTwoHaplotypes, final byte[] ref, final GenomeLoc loc, final ConstrainedMateFixingManager manager, final List<SAMRecord> reads, final double[][] likelihoods ) {

        final SWPairwiseAlignment swConsensus0 = new SWPairwiseAlignment( ref, bestTwoHaplotypes.first.bases, SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );
        final SWPairwiseAlignment swConsensus1 = new SWPairwiseAlignment( ref, bestTwoHaplotypes.second.bases, SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );
        final Consensus consensus0 = new Consensus(bestTwoHaplotypes.first.bases, swConsensus0.getCigar(), swConsensus0.getAlignmentStart2wrt1());
        final Consensus consensus1 = new Consensus(bestTwoHaplotypes.second.bases, swConsensus1.getCigar(), swConsensus1.getAlignmentStart2wrt1());

        int iii = 0;
        for( final SAMRecord read : reads ) {
            final Consensus bestConsensus = ( likelihoods[iii][0] > likelihoods[iii][1] ? consensus0 : consensus1 );
            final AlignedRead aRead = new AlignedRead( read );
            bestConsensus.cigar = AlignmentUtils.leftAlignIndel(bestConsensus.cigar, ref, bestConsensus.str, bestConsensus.positionOnReference, bestConsensus.positionOnReference);
            Pair<Integer, Integer> altAlignment = findBestOffset(bestConsensus.str, aRead, loc.getStart());
            updateRead(bestConsensus.cigar, bestConsensus.positionOnReference, altAlignment.second, aRead, loc.getStart());
            aRead.finalizeUpdate();
            manager.addRead(aRead.getRead(), true);
            iii++;
        }
    }

    // private functions copied from IndelRealigner
    private Pair<Integer, Integer> findBestOffset(final byte[] ref, final AlignedRead read, final int leftmostIndex) {

        // optimization: try the most likely alignment first (to get a low score to beat)
        int originalAlignment = read.getOriginalAlignmentStart() - leftmostIndex;
        int bestScore = mismatchQualitySumIgnoreCigar(read, ref, originalAlignment, Integer.MAX_VALUE);
        int bestIndex = originalAlignment;

        // optimization: we can't get better than 0, so we can quit now
        if ( bestScore == 0 )
            return new Pair<Integer, Integer>(bestIndex, 0);

        // optimization: the correct alignment shouldn't be too far from the original one (or else the read wouldn't have aligned in the first place)
        for ( int i = 0; i < originalAlignment; i++ ) {
            int score = mismatchQualitySumIgnoreCigar(read, ref, i, bestScore);
            if ( score < bestScore ) {
                bestScore = score;
                bestIndex = i;
            }
            // optimization: we can't get better than 0, so we can quit now
            if ( bestScore == 0 )
                return new Pair<Integer, Integer>(bestIndex, 0);
        }

        final int maxPossibleStart = ref.length - read.getReadLength();
        for ( int i = originalAlignment + 1; i <= maxPossibleStart; i++ ) {
            int score = mismatchQualitySumIgnoreCigar(read, ref, i, bestScore);
            if ( score < bestScore ) {
                bestScore = score;
                bestIndex = i;
            }
            // optimization: we can't get better than 0, so we can quit now
            if ( bestScore == 0 )
                return new Pair<Integer, Integer>(bestIndex, 0);
        }

        return new Pair<Integer, Integer>(bestIndex, bestScore);
    }

    private static int mismatchQualitySumIgnoreCigar(final AlignedRead aRead, final byte[] refSeq, int refIndex, int quitAboveThisValue) {
        final byte[] readSeq = aRead.getReadBases();
        final byte[] quals = aRead.getBaseQualities();
        int sum = 0;
        for (int readIndex = 0 ; readIndex < readSeq.length ; refIndex++, readIndex++ ) {
            if ( refIndex >= refSeq.length ) {
                sum += 99;
                // optimization: once we pass the threshold, stop calculating
                if ( sum > quitAboveThisValue )
                    return sum;
            } else {
                byte refChr = refSeq[refIndex];
                byte readChr = readSeq[readIndex];
                if ( !BaseUtils.isRegularBase(readChr) || !BaseUtils.isRegularBase(refChr) )
                    continue; // do not count Ns/Xs/etc ?
                if ( readChr != refChr ) {
                    sum += (int)quals[readIndex];
                    // optimization: once we pass the threshold, stop calculating
                    if ( sum > quitAboveThisValue )
                        return sum;
                }
            }
        }
        return sum;
    }

    private boolean updateRead(final Cigar altCigar, final int altPosOnRef, final int myPosOnAlt, final AlignedRead aRead, final int leftmostIndex) {
        Cigar readCigar = new Cigar();

        // special case: there is no indel
        if ( altCigar.getCigarElements().size() == 1 ) {
            aRead.setAlignmentStart(leftmostIndex + myPosOnAlt);
            readCigar.add(new CigarElement(aRead.getReadLength(), CigarOperator.M));
            aRead.setCigar(readCigar);
            return true;
        }

        CigarElement altCE1 = altCigar.getCigarElement(0);
        CigarElement altCE2 = altCigar.getCigarElement(1);

        int leadingMatchingBlockLength = 0; // length of the leading M element or 0 if the leading element is I

        CigarElement indelCE;
        if ( altCE1.getOperator() == CigarOperator.I  ) {
            indelCE=altCE1;
            if ( altCE2.getOperator() != CigarOperator.M  ) {
                return false;
            }
        }
        else {
            if ( altCE1.getOperator() != CigarOperator.M  ) {
                return false;
            }
            if ( altCE2.getOperator() == CigarOperator.I  || altCE2.getOperator() == CigarOperator.D ) {
                indelCE=altCE2;
            } else {
                return false;
            }
            leadingMatchingBlockLength = altCE1.getLength();
        }

        // the easiest thing to do is to take each case separately
        int endOfFirstBlock = altPosOnRef + leadingMatchingBlockLength;
        boolean sawAlignmentStart = false;

        // for reads starting before the indel
        if ( myPosOnAlt < endOfFirstBlock) {
            aRead.setAlignmentStart(leftmostIndex + myPosOnAlt);
            sawAlignmentStart = true;

            // for reads ending before the indel
            if ( myPosOnAlt + aRead.getReadLength() <= endOfFirstBlock) {
                //readCigar.add(new CigarElement(aRead.getReadLength(), CigarOperator.M));
                //aRead.setCigar(readCigar);
                aRead.setCigar(null); // reset to original alignment
                return true;
            }
            readCigar.add(new CigarElement(endOfFirstBlock - myPosOnAlt, CigarOperator.M));
        }

        // forward along the indel
        //int indelOffsetOnRef = 0, indelOffsetOnRead = 0;
        if ( indelCE.getOperator() == CigarOperator.I ) {
            // for reads that end in an insertion
            if ( myPosOnAlt + aRead.getReadLength() < endOfFirstBlock + indelCE.getLength() ) {
                int partialInsertionLength = myPosOnAlt + aRead.getReadLength() - endOfFirstBlock;
                // if we also started inside the insertion, then we need to modify the length
                if ( !sawAlignmentStart )
                    partialInsertionLength = aRead.getReadLength();
                readCigar.add(new CigarElement(partialInsertionLength, CigarOperator.I));
                aRead.setCigar(readCigar);
                return true;
            }

            // for reads that start in an insertion
            if ( !sawAlignmentStart && myPosOnAlt < endOfFirstBlock + indelCE.getLength() ) {
                aRead.setAlignmentStart(leftmostIndex + endOfFirstBlock);
                readCigar.add(new CigarElement(indelCE.getLength() - (myPosOnAlt - endOfFirstBlock), CigarOperator.I));
                //indelOffsetOnRead = myPosOnAlt - endOfFirstBlock;
                sawAlignmentStart = true;
            } else if ( sawAlignmentStart ) {
                readCigar.add(indelCE);
                //indelOffsetOnRead = indelCE.getLength();
            }
        } else if ( indelCE.getOperator() == CigarOperator.D ) {
            if ( sawAlignmentStart )
                readCigar.add(indelCE);
            //indelOffsetOnRef = indelCE.getLength();
        }

        // for reads that start after the indel
        if ( !sawAlignmentStart ) {
            //aRead.setAlignmentStart(leftmostIndex + myPosOnAlt + indelOffsetOnRef - indelOffsetOnRead);
            //readCigar.add(new CigarElement(aRead.getReadLength(), CigarOperator.M));
            //aRead.setCigar(readCigar);
            aRead.setCigar(null); // reset to original alignment
            return true;
        }

        int readRemaining = aRead.getReadBases().length;
        for ( CigarElement ce : readCigar.getCigarElements() ) {
            if ( ce.getOperator() != CigarOperator.D )
                readRemaining -= ce.getLength();
        }
        if ( readRemaining > 0 )
            readCigar.add(new CigarElement(readRemaining, CigarOperator.M));
        aRead.setCigar(readCigar);

        return true;
    }


    // private classes copied from IndelRealigner
    private class AlignedRead {
        private final SAMRecord read;
        private byte[] readBases = null;
        private byte[] baseQuals = null;
        private Cigar newCigar = null;
        private int newStart = -1;
        private int mismatchScoreToReference = 0;
        private long alignerMismatchScore = 0;

        public AlignedRead(SAMRecord read) {
            this.read = read;
            mismatchScoreToReference = 0;
        }

        public SAMRecord getRead() {
               return read;
        }

        public int getReadLength() {
            return readBases != null ? readBases.length : read.getReadLength();
        }

        public byte[] getReadBases() {
            if ( readBases == null )
                getUnclippedBases();
            return readBases;
        }

        public byte[] getBaseQualities() {
            if ( baseQuals == null )
                getUnclippedBases();
            return baseQuals;
        }

        // pull out the bases that aren't clipped out
        private void getUnclippedBases() {
            readBases = new byte[getReadLength()];
            baseQuals = new byte[getReadLength()];
            byte[] actualReadBases = read.getReadBases();
            byte[] actualBaseQuals = read.getBaseQualities();
            int fromIndex = 0, toIndex = 0;

            for ( CigarElement ce : read.getCigar().getCigarElements() ) {
                int elementLength = ce.getLength();
                switch ( ce.getOperator() ) {
                    case S:
                        fromIndex += elementLength;
                        break;
                    case M:
                    case I:
                        System.arraycopy(actualReadBases, fromIndex, readBases, toIndex, elementLength);
                        System.arraycopy(actualBaseQuals, fromIndex, baseQuals, toIndex, elementLength);
                        fromIndex += elementLength;
                        toIndex += elementLength;
                    default:
                        break;
                }
            }

            // if we got clipped, trim the array
            if ( fromIndex != toIndex ) {
                byte[] trimmedRB = new byte[toIndex];
                byte[] trimmedBQ = new byte[toIndex];
                System.arraycopy(readBases, 0, trimmedRB, 0, toIndex);
                System.arraycopy(baseQuals, 0, trimmedBQ, 0, toIndex);
                readBases = trimmedRB;
                baseQuals = trimmedBQ;
            }
        }

        public Cigar getCigar() {
            return (newCigar != null ? newCigar : read.getCigar());
        }

        public void setCigar(Cigar cigar) {
            setCigar(cigar, true);
        }

        // tentatively sets the new Cigar, but it needs to be confirmed later
        public void setCigar(Cigar cigar, boolean fixClippedCigar) {
            if ( cigar == null ) {
                newCigar = null;
                return;
            }

            if ( fixClippedCigar && getReadBases().length < read.getReadLength() )
                cigar = reclipCigar(cigar);

            // no change?
            if ( read.getCigar().equals(cigar) ) {
                newCigar = null;
                return;
            }

            // no indel?
            String str = cigar.toString();
            if ( !str.contains("D") && !str.contains("I") ) {
                //    newCigar = null;
                //    return;
            }

            newCigar = cigar;
        }

        // pull out the bases that aren't clipped out
        private Cigar reclipCigar(Cigar cigar) {
            return reclipCigar(cigar, read);
        }

        private boolean isClipOperator(CigarOperator op) {
            return op == CigarOperator.S || op == CigarOperator.H || op == CigarOperator.P;
        }

        protected Cigar reclipCigar(Cigar cigar, SAMRecord read) {
            ArrayList<CigarElement> elements = new ArrayList<CigarElement>();

            int i = 0;
            int n = read.getCigar().numCigarElements();
            while ( i < n && isClipOperator(read.getCigar().getCigarElement(i).getOperator()) )
                elements.add(read.getCigar().getCigarElement(i++));

            elements.addAll(cigar.getCigarElements());

            i++;
            while ( i < n && !isClipOperator(read.getCigar().getCigarElement(i).getOperator()) )
                i++;

            while ( i < n && isClipOperator(read.getCigar().getCigarElement(i).getOperator()) )
                elements.add(read.getCigar().getCigarElement(i++));

            return new Cigar(elements);
        }

        // tentatively sets the new start, but it needs to be confirmed later
        public void setAlignmentStart(int start) {
            newStart = start;
        }

        public int getAlignmentStart() {
            return (newStart != -1 ? newStart : read.getAlignmentStart());
        }

        public int getOriginalAlignmentStart() {
            return read.getAlignmentStart();
        }

        // finalizes the changes made.
        // returns true if this record actually changes, false otherwise
        public boolean finalizeUpdate() {
            // if we haven't made any changes, don't do anything
            if ( newCigar == null )
                return false;
            if ( newStart == -1 )
                newStart = read.getAlignmentStart();

            read.setCigar(newCigar);
            read.setAlignmentStart(newStart);

            return true;
        }

        public void setMismatchScoreToReference(int score) {
            mismatchScoreToReference = score;
        }

        public int getMismatchScoreToReference() {
            return mismatchScoreToReference;
        }

        public void setAlignerMismatchScore(long score) {
            alignerMismatchScore = score;
        }

        public long getAlignerMismatchScore() {
            return alignerMismatchScore;
        }
    }

    private static class Consensus {
        public final byte[] str;
        public final ArrayList<Pair<Integer, Integer>> readIndexes;
        public final int positionOnReference;
        public int mismatchSum;
        public Cigar cigar;

        public Consensus(byte[] str, Cigar cigar, int positionOnReference) {
            this.str = str;
            this.cigar = cigar;
            this.positionOnReference = positionOnReference;
            mismatchSum = 0;
            readIndexes = new ArrayList<Pair<Integer, Integer>>();
        }

        @Override
        public boolean equals(Object o) {
            return ( this == o || (o instanceof Consensus && Arrays.equals(this.str,(((Consensus)o).str)) ) );
        }

        public boolean equals(Consensus c) {
            return ( this == c || Arrays.equals(this.str,c.str) ) ;
        }

        @Override
        public int hashCode() {
            return Arrays.hashCode(this.str);
        }
    }


}