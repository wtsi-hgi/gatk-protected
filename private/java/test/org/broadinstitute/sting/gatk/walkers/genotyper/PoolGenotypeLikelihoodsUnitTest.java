package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.GenotypeLikelihoods;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: delangel
 * Date: 4/3/12
 * Time: 7:26 PM
 * To change this template use File | Settings | File Templates.
 */
public class PoolGenotypeLikelihoodsUnitTest {


    @Test
    public void testStoringLikelihoodElements() {


        // basic test storing a given PL vector in a PoolGenotypeLikelihoods object and then retrieving it back

        int ploidy = 20;
        int numAltAlleles = 4;
        int res = GenotypeLikelihoods.calculateNumLikelihoods(numAltAlleles, ploidy);
 //       System.out.format("Alt Alleles: %d, Ploidy: %d, #Likelihoods: %d\n", numAltAlleles, ploidy, res);

        List<Allele> alleles = new ArrayList<Allele>();
        alleles.add(Allele.create("T",true));
        alleles.add(Allele.create("C",false));
        alleles.add(Allele.create("A",false));
        alleles.add(Allele.create("G",false));

        double[] gls = new double[res];

        for (int k=0; k < gls.length; k++)
            gls[k]= (double)k;

        PoolGenotypeLikelihoods gl = new PoolSNPGenotypeLikelihoods(alleles, gls,ploidy, null, true);
        double[] glnew = gl.getLikelihoods();

        Assert.assertEquals(gls, glnew);
    }

    @Test
    public void testIndexIterator() {
        int[] seed = new int[]{1,2,3,4};
        PoolGenotypeLikelihoods.SumIterator iterator = runIterator(seed,-1);
        // Assert.assertTrue(compareIntArrays(iterator.getCurrentVector(), seed));
        Assert.assertEquals(iterator.getLinearIndex(),prod(seed)-1);

        seed = new int[]{1,0,1,1};
        iterator = runIterator(seed,-1);
        //  Assert.assertTrue(compareIntArrays(iterator.getCurrentVector(), seed));
        Assert.assertEquals(iterator.getLinearIndex(),prod(seed)-1);

        seed = new int[]{5};
        iterator = runIterator(seed,-1);
        //  Assert.assertTrue(compareIntArrays(iterator.getCurrentVector(), seed));
        Assert.assertEquals(iterator.getLinearIndex(),prod(seed)-1);

        // Diploid, # alleles = 4
        seed = new int[]{2,2,2,2};
        iterator = runIterator(seed,2);
        //  Assert.assertTrue(compareIntArrays(iterator.getCurrentVector(), seed));
        Assert.assertEquals(iterator.getLinearIndex(),9);

        // Diploid, # alleles = 2
        seed = new int[]{2,2};
        iterator = runIterator(seed,2);
        //  Assert.assertTrue(compareIntArrays(iterator.getCurrentVector(), seed));
        Assert.assertEquals(iterator.getLinearIndex(),2);

        // Diploid, # alleles = 3
        seed = new int[]{2,2,2};
        iterator = runIterator(seed,2);
        //  Assert.assertTrue(compareIntArrays(iterator.getCurrentVector(), seed));
        Assert.assertEquals(iterator.getLinearIndex(),5);

        // Triploid, # alleles = 2
        seed = new int[]{3,3};
        iterator = runIterator(seed,3);
        //  Assert.assertTrue(compareIntArrays(iterator.getCurrentVector(), seed));
        Assert.assertEquals(iterator.getLinearIndex(),3);
        // Triploid, # alleles = 3
        seed = new int[]{3,3,3};
        iterator = runIterator(seed,3);
        //  Assert.assertTrue(compareIntArrays(iterator.getCurrentVector(), seed));
        Assert.assertEquals(iterator.getLinearIndex(),9);

        // Triploid, # alleles = 4
        seed = new int[]{3,3,3,3};
        iterator = runIterator(seed,3);
        //  Assert.assertTrue(compareIntArrays(iterator.getCurrentVector(), seed));
        Assert.assertEquals(iterator.getLinearIndex(),9);



    }

    private PoolGenotypeLikelihoods.SumIterator runIterator(int[] seed, int restrictSumTo) {
        PoolGenotypeLikelihoods.SumIterator iterator = new PoolGenotypeLikelihoods.SumIterator(seed, restrictSumTo);

        while(iterator.hasNext()) {
            System.out.format("\n%d:",iterator.getLinearIndex());
            int[] a =  iterator.getCurrentVector();
            for (int i=0; i < seed.length; i++)
                System.out.format("%d ",a[i]);

            iterator.next();
        }

        return iterator;

    }

    private static int prod(int[] x) {
        int prod = 1;
        for (int k=0; k < x.length; k++) {
            prod *= (1+x[k]);
        }
        return prod;
    }

}
