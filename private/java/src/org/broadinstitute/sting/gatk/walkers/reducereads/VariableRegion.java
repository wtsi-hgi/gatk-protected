package org.broadinstitute.sting.gatk.walkers.reducereads;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

/**
 * Created by IntelliJ IDEA.
 * User: roger
 * Date: 8/8/11
 * Time: 12:02 PM
 * To change this template use File | Settings | File Templates.
 */
public class VariableRegion {
    public int start;
    public int end;

    public VariableRegion( int Start, int End) {
            if ( Start <= End ) {
                this.start = Start;
                this.end = End;
            }
            else
                throw new ReviewedStingException("Variable Region must have start before end");
        }

    public VariableRegion merge( VariableRegion that ) {
        return new VariableRegion( Math.min(this.start,that.start), Math.max(this.end,that.end));
    }
}
