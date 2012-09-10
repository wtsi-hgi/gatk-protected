package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: rpoplin
 * Date: 8/28/12
 */

public class DelocalizedBaseRecalibratorUnitTest {

    @Test
    public void basicDBQSRFractionalErrorTestEnd() {
        byte[] baq = "@@@@@@@FGH".getBytes();
        boolean[] skip = new boolean[baq.length];
        Arrays.fill(skip, false);
        int[] errors = new int[baq.length];
        Arrays.fill(errors, 0);
        errors[7] = 1;
        double[] answer = new double[baq.length];
        Arrays.fill(answer,0.0);
        answer[6] = answer[7] = answer[8] = answer[9] = 1.0 / 4.0;
        double[] result = DelocalizedBaseRecalibrator.calculateFractionalErrorArray(skip, errors, baq);
        for( int iii = 0; iii < answer.length; iii++) {
            Assert.assertEquals(result[iii], answer[iii], 1E-6);
        }
    }

    @Test
    public void basicDBQSRFractionalErrorTestStart() {
        byte[] baq = "FFF@@@@@@@".getBytes();
        boolean[] skip = new boolean[baq.length];
        Arrays.fill(skip, false);
        int[] errors = new int[baq.length];
        Arrays.fill(errors, 0);
        errors[2] = 1;
        double[] answer = new double[baq.length];
        Arrays.fill(answer,0.0);
        answer[0] = answer[1] = answer[2] = answer[3] = 1.0 / 4.0;
        double[] result = DelocalizedBaseRecalibrator.calculateFractionalErrorArray(skip, errors, baq);
        for( int iii = 0; iii < answer.length; iii++) {
            Assert.assertEquals(result[iii], answer[iii], 1E-6);
        }
    }

    @Test
    public void basicDBQSRFractionalErrorTestNoBAQ() {
        byte[] baq = "@@@@@@@@@@".getBytes();
        boolean[] skip = new boolean[baq.length];
        Arrays.fill(skip, false);
        int[] errors = new int[baq.length];
        Arrays.fill(errors, 0);
        errors[7] = 1;
        double[] answer = new double[baq.length];
        Arrays.fill(answer,0.0);
        answer[7] = 1.0;
        double[] result = DelocalizedBaseRecalibrator.calculateFractionalErrorArray(skip, errors, baq);
        for( int iii = 0; iii < answer.length; iii++) {
            Assert.assertEquals(result[iii], answer[iii], 1E-6);
        }
    }

    @Test
    public void basicDBQSRFractionalErrorTestBAQOffset() {
        byte[] baq = "@FGH@@@@@@".getBytes();
        boolean[] skip = new boolean[baq.length];
        Arrays.fill(skip, false);
        int[] errors = new int[baq.length];
        Arrays.fill(errors, 0);
        errors[7] = 1;
        double[] answer = new double[baq.length];
        Arrays.fill(answer,0.0);
        answer[7] = 1.0;
        double[] result = DelocalizedBaseRecalibrator.calculateFractionalErrorArray(skip, errors, baq);
        for( int iii = 0; iii < answer.length; iii++) {
            Assert.assertEquals(result[iii], answer[iii], 1E-6);
        }
    }

    @Test
    public void basicDBQSRFractionalErrorTestMiddle() {
        byte[] baq = "@@@FGH@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@".getBytes();
        boolean[] skip = new boolean[baq.length];
        Arrays.fill(skip, false);
        int[] errors = new int[baq.length];
        Arrays.fill(errors, 0);
        errors[4] = 1;
        double[] answer = new double[baq.length];
        Arrays.fill(answer,0.0);
        answer[2] = answer[3] = answer[4] = answer[5] = answer[6] = 1.0 / 5.0;
        double[] result = DelocalizedBaseRecalibrator.calculateFractionalErrorArray(skip, errors, baq);
        for( int iii = 0; iii < answer.length; iii++) {
            Assert.assertEquals(result[iii], answer[iii], 1E-6);
        }
    }

    @Test
    public void basicDBQSRFractionalErrorTestMiddleSkip() {
        byte[] baq = "@@@FGH@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@".getBytes();
        boolean[] skip = new boolean[baq.length];
        Arrays.fill(skip, false);
        skip[4] = true;
        int[] errors = new int[baq.length];
        Arrays.fill(errors, 0);
        errors[4] = 1;
        double[] answer = new double[baq.length];
        Arrays.fill(answer,0.0);
        double[] result = DelocalizedBaseRecalibrator.calculateFractionalErrorArray(skip, errors, baq);
        for( int iii = 0; iii < answer.length; iii++) {
            Assert.assertEquals(result[iii], answer[iii], 1E-6);
        }
    }

    @Test
    public void basicDBQSRFractionalErrorTestMiddleUninformativeSkip() {
        byte[] baq = "@@@FGH@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@".getBytes();
        boolean[] skip = new boolean[baq.length];
        Arrays.fill(skip, false);
        skip[3] = true;
        skip[13] = true;
        int[] errors = new int[baq.length];
        Arrays.fill(errors, 0);
        errors[4] = 1;
        double[] answer = new double[baq.length];
        Arrays.fill(answer,0.0);
        answer[2] = answer[3] = answer[4] = answer[5] = answer[6] = 1.0 / 5.0;
        double[] result = DelocalizedBaseRecalibrator.calculateFractionalErrorArray(skip, errors, baq);
        for( int iii = 0; iii < answer.length; iii++) {
            Assert.assertEquals(result[iii], answer[iii], 1E-6);
        }
    }

}
