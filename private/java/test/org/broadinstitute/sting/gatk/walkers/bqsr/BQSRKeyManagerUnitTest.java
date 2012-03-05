package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.BitSet;

/**
 * @author Mauricio Carneiro
 * @since 3/7/12
 */
public class BQSRKeyManagerUnitTest {
    RecalibrationArgumentCollection RAC;

    @BeforeClass
    public void init() {
        RAC = new RecalibrationArgumentCollection();
    }

    @Test(enabled = true)
    public void testCombineBitSets() {
        final ArrayList<Covariate> covariateList = new ArrayList<Covariate>();
        covariateList.add(new ReadGroupCovariate());
        covariateList.add(new QualityScoreCovariate());
        covariateList.add(new CycleCovariate());
        covariateList.add(new ContextCovariate());

        int readLength = 1000;
        GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(ReadUtils.createRandomReadBases(readLength, true), ReadUtils.createRandomReadQuals(readLength), readLength + "M");
        read.setReadGroup(new GATKSAMReadGroupRecord("MY.ID"));
        read.getReadGroup().setPlatform("illumina");

        final BitSet[][][] covariateKeys = new BitSet[covariateList.size()][RecalDataManager.BaseRecalibrationType.values().length][];
        int i = 0;
        for (Covariate cov : covariateList) {
            cov.initialize(RAC);
            covariateKeys[i][RecalDataManager.BaseRecalibrationType.BASE_SUBSTITUTION.index] = cov.getValues(read).getMismatches();
            covariateKeys[i][RecalDataManager.BaseRecalibrationType.BASE_INSERTION.index] = cov.getValues(read).getInsertions();
            covariateKeys[i][RecalDataManager.BaseRecalibrationType.BASE_DELETION.index] = cov.getValues(read).getDeletions();
            i++;
        }

        BQSRKeyManager keyManager = new BQSRKeyManager(covariateList);

        for (int l = 0; l < readLength; l++) {
            for (int eventType = 0; eventType < RecalDataManager.BaseRecalibrationType.values().length; eventType++) {
                BitSet[] keySet = new BitSet[covariateList.size()];
                Object[] expected = new Object[covariateList.size()];
                for (int j = 0; j < covariateList.size(); j++) {
                    keySet[j] = covariateKeys[j][eventType][l];
                    expected[j] = covariateList.get(j).keyFromBitSet(keySet[j]);
                }

                BitSet hashKey = keyManager.bitSetFrom(keySet, RecalDataManager.BaseRecalibrationType.eventFrom(eventType));
                Object[] actual = keyManager.keySetFrom(hashKey);

                for (Object o : actual)
                    System.out.print(o + ", ");
                System.out.println();
                for (Object o : expected)
                    System.out.print(o + ", ");
                System.out.println();
                System.out.println();

                for (int k = 0; k < expected.length; k++)
                    Assert.assertEquals(actual[k], expected[k]);
            }
        }
    }

}
