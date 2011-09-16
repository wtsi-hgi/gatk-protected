package org.broadinstitute.sting.gatk.walkers.poolcaller;

import org.broadinstitute.sting.utils.MathUtils;

/**
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 7/22/11
 * Time: 3:48 PM
 *
 * This class implements the basic probability model used by Allele Count and Error Models.
 */
public class ProbabilityModel {
    protected double [] model;

    public String toString() {
        String result = "( ";
        boolean skipComma = true;
        for (double v : model) {
            if (skipComma) {
                skipComma = false;
            }
            else {
                result += ", ";
            }
            result += String.format("%.4f", v);
        }
        return result + " )";
    }

    public int size() {
        return model.length;
    }

    public double getCumulativeSum (int upTo) {
        return MathUtils.log10CumulativeSumLog10(model, upTo);
    }

    public double getMaximumLikelihood() {
        return MathUtils.arrayMax(model);
    }

    public int getMaximumLikelihoodIndex() {
        return MathUtils.maxElementIndex(model);

    }

}
