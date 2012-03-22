package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 3/15/12
 */

import net.sf.picard.reference.ReferenceSequenceFile;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * Unit tests for GenotypingEngine
 */
public class GenotypingEngineUnitTest extends BaseTest {

    private static ReferenceSequenceFile seq;
    private GenomeLocParser genomeLocParser;

    @BeforeClass
    public void init() throws FileNotFoundException {
        // sequence
        seq = new CachingIndexedFastaSequenceFile(new File(b37KGReference));
        genomeLocParser = new GenomeLocParser(seq);
    }

    @Test
    public void testFindHomVarEventAllelesInSample() {
        final List<Allele> eventAlleles = new ArrayList<Allele>();
        eventAlleles.add( Allele.create("A", true) );
        eventAlleles.add( Allele.create("C", false) );
        final List<Allele> haplotypeAlleles = new ArrayList<Allele>();
        haplotypeAlleles.add( Allele.create("AATA", true) );
        haplotypeAlleles.add( Allele.create("AACA", false) );
        haplotypeAlleles.add( Allele.create("CATA", false) );
        haplotypeAlleles.add( Allele.create("CACA", false) );
        final List<Allele> haplotypeAllelesForSample = new ArrayList<Allele>();
        haplotypeAllelesForSample.add( Allele.create("CATA", false) );
        haplotypeAllelesForSample.add( Allele.create("CACA", false) );
        final ArrayList<ArrayList<Integer>> alleleMapper = new ArrayList<ArrayList<Integer>>();
        ArrayList<Integer> Aallele = new ArrayList<Integer>();
        Aallele.add(0);
        Aallele.add(1);
        ArrayList<Integer> Callele = new ArrayList<Integer>();
        Callele.add(2);
        Callele.add(3);
        alleleMapper.add(Aallele);
        alleleMapper.add(Callele);
        final List<Allele> eventAllelesForSample = new ArrayList<Allele>();
        eventAllelesForSample.add( Allele.create("C", false) );
        eventAllelesForSample.add( Allele.create("C", false) );

        if(!compareAlleleLists(eventAllelesForSample, GenotypingEngine.findEventAllelesInSample(eventAlleles, haplotypeAlleles, haplotypeAllelesForSample, alleleMapper))) {
            logger.warn("calc alleles = " + GenotypingEngine.findEventAllelesInSample(eventAlleles, haplotypeAlleles, haplotypeAllelesForSample, alleleMapper));
            logger.warn("expected alleles = " + eventAllelesForSample);
        }
        Assert.assertTrue(compareAlleleLists(eventAllelesForSample, GenotypingEngine.findEventAllelesInSample(eventAlleles, haplotypeAlleles, haplotypeAllelesForSample, alleleMapper)));
    }

    @Test
    public void testFindHetEventAllelesInSample() {
        final List<Allele> eventAlleles = new ArrayList<Allele>();
        eventAlleles.add( Allele.create("A", true) );
        eventAlleles.add( Allele.create("C", false) );
        eventAlleles.add( Allele.create("T", false) );
        final List<Allele> haplotypeAlleles = new ArrayList<Allele>();
        haplotypeAlleles.add( Allele.create("AATA", true) );
        haplotypeAlleles.add( Allele.create("AACA", false) );
        haplotypeAlleles.add( Allele.create("CATA", false) );
        haplotypeAlleles.add( Allele.create("CACA", false) );
        haplotypeAlleles.add( Allele.create("TACA", false) );
        haplotypeAlleles.add( Allele.create("TTCA", false) );
        haplotypeAlleles.add( Allele.create("TTTA", false) );
        final List<Allele> haplotypeAllelesForSample = new ArrayList<Allele>();
        haplotypeAllelesForSample.add( Allele.create("TTTA", false) );
        haplotypeAllelesForSample.add( Allele.create("AATA", true) );
        final ArrayList<ArrayList<Integer>> alleleMapper = new ArrayList<ArrayList<Integer>>();
        ArrayList<Integer> Aallele = new ArrayList<Integer>();
        Aallele.add(0);
        Aallele.add(1);
        ArrayList<Integer> Callele = new ArrayList<Integer>();
        Callele.add(2);
        Callele.add(3);
        ArrayList<Integer> Tallele = new ArrayList<Integer>();
        Tallele.add(4);
        Tallele.add(5);
        Tallele.add(6);
        alleleMapper.add(Aallele);
        alleleMapper.add(Callele);
        alleleMapper.add(Tallele);
        final List<Allele> eventAllelesForSample = new ArrayList<Allele>();
        eventAllelesForSample.add( Allele.create("A", true) );
        eventAllelesForSample.add( Allele.create("T", false) );

        if(!compareAlleleLists(eventAllelesForSample, GenotypingEngine.findEventAllelesInSample(eventAlleles, haplotypeAlleles, haplotypeAllelesForSample, alleleMapper))) {
            logger.warn("calc alleles = " + GenotypingEngine.findEventAllelesInSample(eventAlleles, haplotypeAlleles, haplotypeAllelesForSample, alleleMapper));
            logger.warn("expected alleles = " + eventAllelesForSample);
        }
        Assert.assertTrue(compareAlleleLists(eventAllelesForSample, GenotypingEngine.findEventAllelesInSample(eventAlleles, haplotypeAlleles, haplotypeAllelesForSample, alleleMapper)));
    }

    private boolean compareAlleleLists(List<Allele> l1, List<Allele> l2) {
        if( l1.size() != l2.size() ) {
            return false; // sanity check
        }

        for( int i=0; i < l1.size(); i++ ){
            if ( !l2.contains(l1.get(i)) )
                return false;
        }
        return true;
    }

    
    private class BasicGenotypingTestProvider extends TestDataProvider {
        byte[] ref;
        byte[] hap;
        HashMap<Integer,Byte> expected;

        public BasicGenotypingTestProvider(String refString, String hapString, HashMap<Integer, Byte> expected) {
            super(BasicGenotypingTestProvider.class, String.format("Haplotype to VCF test: ref = %s, alignment = %s", refString,hapString));
            ref = refString.getBytes();
            hap = hapString.getBytes();
            this.expected = expected;
        }
        
        public HashMap<Integer,VariantContext> calcAlignment() {
            return GenotypingEngine.generateVCsFromAlignment( new SWPairwiseAlignment(ref, hap), ref, hap, genomeLocParser.createGenomeLoc("4",1,1+ref.length), "name");
        }
    }

    @DataProvider(name = "BasicGenotypingTestProvider")
    public Object[][] makeBasicGenotypingTests() {

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            HashMap<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(1 + contextSize, (byte)'M');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGACTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            HashMap<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(2 + contextSize, (byte)'M');
            map.put(21 + contextSize, (byte)'M');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG", "ATCTCGCATCGCGAGCATCGCCTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            HashMap<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(1 + contextSize, (byte)'M');
            map.put(20 + contextSize, (byte)'I');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGACACTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            HashMap<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(1 + contextSize, (byte)'M');
            map.put(20 + contextSize, (byte)'D');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGCTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{1,5,9,24,36} ) {
            HashMap<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(1, (byte)'M');
            map.put(20, (byte)'D');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider("AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGCTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            HashMap<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(2 + contextSize, (byte)'M');
            map.put(20 + contextSize, (byte)'I');
            map.put(30 + contextSize, (byte)'D');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "ACCTCGCATCGCGAGCATCGTTACTAGCCGATG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            HashMap<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(1 + contextSize, (byte)'M');
            map.put(20 + contextSize, (byte)'D');
            map.put(28 + contextSize, (byte)'M');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGCTAGCCCATAG", map);
        }

        return BasicGenotypingTestProvider.getTests(BasicGenotypingTestProvider.class);
    }
    
    @Test(dataProvider = "BasicGenotypingTestProvider", enabled = true)
    public void testHaplotypeToVCF(BasicGenotypingTestProvider cfg) {
        HashMap<Integer,VariantContext> calculatedMap = cfg.calcAlignment();
        HashMap<Integer,Byte> expectedMap = cfg.expected;
        logger.warn(String.format("Test: %s", cfg.toString()));
        if(!compareVCMaps(calculatedMap, expectedMap)) {
            logger.warn("calc map = " + calculatedMap);
            logger.warn("expected map = " + expectedMap);
        }
        Assert.assertTrue(compareVCMaps(calculatedMap, expectedMap));
    }
    
    /**
     * Private function to compare HashMap of VCs, it only checks the types and start locations of the VariantContext
     */
    private boolean compareVCMaps(HashMap<Integer, VariantContext> calc, HashMap<Integer, Byte> expected) {
        if( !calc.keySet().equals(expected.keySet()) ) { return false; } // sanity check
        for( Integer loc : expected.keySet() ) {
            Byte type = expected.get(loc);
            switch( type ) {
                case 'I':
                    if( !calc.get(loc).isSimpleInsertion() ) { return false; }
                    break;
                case 'D':
                    if( !calc.get(loc).isSimpleDeletion() ) { return false; }
                    break;
                case 'M':
                    if( !(calc.get(loc).isMNP() || calc.get(loc).isSNP()) ) { return false; }
                    break;
                default:
                    return false;
            }
        }
        return true;
    }
}