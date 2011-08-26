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

import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SWPairwiseAlignment;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.InferredGeneticContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;

public class GenotypingEngine {

    // Smith-Waterman parameters copied from IndelRealigne
    private static final double SW_MATCH = 30.0;      // 1.0;
    private static final double SW_MISMATCH = -10.0;  //-1.0/3.0;
    private static final double SW_GAP = -10.0;       //-1.0-1.0/3.0;
    private static final double SW_GAP_EXTEND = -2.0; //-1.0/.0;

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
                        System.out.println("Insertion: " + alleles);
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
                    System.out.println( "Deletion: " + alleles);
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
                        if( ref[refPos] != read[readPos] ) {
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
                            System.out.println( "SNP/MNP: " + alleles);
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
                    if( alleles.get(0).equals(alleles.get(1)) ) { // check if alt allese match
                        Genotype gt = new Genotype( "NA12878", alleles );
                        HashMap<String,Genotype> genotypeMap = new HashMap<String,Genotype>();
                        genotypeMap.put("NA12878", gt);
                        vcs.add( VariantContext.modifyGenotypes( vc1, genotypeMap ) );
                    } else { // two alt alleles don't match, and don't call multialleleic records yet
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


        } while ( vcs1Iter.hasNext() || vcs2Iter.hasNext() );

        return vcs;
    }

    public void alignAllHaplotypes( final List<Haplotype> haplotypes, final byte[] ref, final GenomeLoc loc, final StingSAMFileWriter writer, final SAMRecord exampleRead ) {

        int iii = 0;
        for( Haplotype h : haplotypes ) {
            final SWPairwiseAlignment swConsensus = new SWPairwiseAlignment( ref, h.bases, SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );
            exampleRead.setReadName("Haplotype" + iii);
            exampleRead.setReadBases(h.bases);
            exampleRead.setAlignmentStart(loc.getStart() + swConsensus.getAlignmentStart2wrt1());
            exampleRead.setCigar(swConsensus.getCigar());
            byte[] quals = new byte[h.bases.length];
            Arrays.fill(quals, (byte) 25);
            exampleRead.setBaseQualities(quals);
            writer.addAlignment(exampleRead);
            iii++;
        }
                
    }

}