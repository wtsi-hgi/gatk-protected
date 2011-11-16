package org.broadinstitute.sting.gatk.walkers.IndelCountCovariates;

import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 9/26/11
 */

public class ContextCovariate implements Covariate {

    final int CONTEXT_SIZE = 9;
    String allN = "";

    // Initialize any member variables using the command-line arguments passed to the walkers
    public void initialize( final RecalibrationArgumentCollection RAC ) {
        for( int iii = 0; iii < CONTEXT_SIZE; iii++ ) {
            allN += "N";
        }
    }

    public void getValues(GATKSAMRecord read, Comparable[] comparable) {
        byte[] bases = read.getReadBases();
        for(int i = 0; i < read.getReadLength(); i++) {
            comparable[i] = ( i-CONTEXT_SIZE < 0 ? allN : new String(Arrays.copyOfRange(bases,i-CONTEXT_SIZE,i)) );
        }
    }

    // Used to get the covariate's value from input csv file in TableRecalibrationWalker
    public final Comparable getValue( final String str ) {
        return str;
    }

}
