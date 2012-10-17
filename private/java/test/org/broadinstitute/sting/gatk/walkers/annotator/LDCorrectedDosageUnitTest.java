package org.broadinstitute.sting.gatk.walkers.annotator;

import org.testng.Assert;
import org.testng.annotations.Test;
import org.apache.commons.math.linear.*;

import java.util.Arrays;

/**
 * Test out pieces of LDCorrectedUnitTest (in particular the OLS computation and rectify function)
 */
public class LDCorrectedDosageUnitTest {

    private final RealMatrix predictor1 =
          new Array2DRowRealMatrix(new double[][]{ {1.0,1.0,1.0,1.0,1.0},
                                                   {2.0,0.0,2.0,2.0,1.0},
                                                   {0.0,1.0,1.0,1.0,0.0},
                                                   {1.0,1.0,2.0,0.0,2.0}});

    private final RealVector response1 = new ArrayRealVector(new double[]{ 1.0,2.0,2.0,1.0,0.0 });

    private final RealMatrix predictor2 =
          new Array2DRowRealMatrix(( new double[][] { { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 } ,
                                                      { 2.0, 1.0, 2.0, 0.0, 0.0, 1.0, 2.0, 2.0, 0.0 } ,
                                                      { 1.0, 2.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 } ,
                                                      { 2.0, 2.0, 2.0, 0.0, 0.0, 1.0, 2.0, 2.0, 1.0 } ,
                                                      { 0.0, 1.0, 0.0, 1.0, 1.0, 2.0, 0.0, 0.0, 0.0 } ,
                                                      { 0.0, 1.0, 0.0, 1.0, 1.0, 2.0, 1.0, 1.0, 0.0 } ,
                                                      { 2.0, 2.0, 1.0, 0.0, 0.0, 1.0, 0.0, 2.0, 0.0 } }));

    private final RealVector response2 = new ArrayRealVector(new double[] {2.0,2.0,0.0,1.0,0.0,0.0,1.0,2.0,0.0} );

    @Test
    public void testResiduals1() {
        double[] test1Residuals = LDCorrectedDosage.runOLS(predictor1,response1,true).getData()[0];
        double[] test1Expected = {0.6086957,0.3043478,0.1521739,-0.4565217,-0.6086957};
        for ( int idx = 0; idx < test1Residuals.length; idx++ ) {
            Assert.assertEquals(test1Residuals[idx],test1Expected[idx],1e-7);
        }
    }

    @Test
    public void testResiduals2() {
        double[] test2Residuals = LDCorrectedDosage.runOLS(predictor2,response2,false).getData()[0];
        double[] test2Expected = {-0.11764706, 0.23529412,  0.11764706,  0.05882353,  0.17647059, -0.23529412,
                                   0.05882353, -0.05882353, -0.23529412};
        double[] test2ResidualsRectified = LDCorrectedDosage.runOLS(predictor2,response2,true).getData()[0];
        double[] test2ExpectedRectified = { 0.0, 0.23529412, 0.0,  0.05882353, 0.0, -0.23529412,
                                   0.05882353, 0.0, -0.23529412};
        for ( int idx = 0; idx < test2Residuals.length; idx++ ) {
            Assert.assertEquals(test2Residuals[idx],test2Expected[idx],1e-7);
            Assert.assertEquals(test2ResidualsRectified[idx],test2ExpectedRectified[idx],1e-7);
        }
    }

}
