package org.broadinstitute.sting.utils.recalibration;

/**
 * Created with IntelliJ IDEA.
 * User: depristo
 * Date: 8/3/12
 * Time: 9:33 AM
 * To change this template use File | Settings | File Templates.
 */
public final class ContextDatum extends RecalDatum {
    public final static String ROOT_CONTEXT = "";
    public final String context;

    public ContextDatum(final String context, final long observations, final long errors) {
        super(observations, errors, (byte)30); // TODO -- should use default value?
        this.context = context;
    }

    @Override
    public String toString() {
        return context == ROOT_CONTEXT ? "x" : context;
    }

    public String getParentContext() {
        return size() == 1 ? ROOT_CONTEXT : context.substring(0, size() - 1);
    }

    public int size() { return context.length(); }
}
