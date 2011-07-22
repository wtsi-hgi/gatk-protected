package org.broadinstitute.sting.gatk.walkers.replication_validation;

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

}
