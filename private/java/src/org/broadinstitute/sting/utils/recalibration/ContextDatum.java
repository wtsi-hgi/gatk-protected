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
        return isRootContext() ? "x" : context;
    }

    /**
     * What's the parent context of this ContextDatum?  The parent is the context excluding
     * the last base of this parent's context.  Because contexts run from 0 [more context] to N [least]
     * the parent context is the substring from 1 to N.
     *
     * getParent(ACGT) => CGT
     * getParent(CGT) => GT
     * getParent(GT) => T
     * getParent(T) => ROOT_CONTEXT
     * getParent(ROOT_CONTEXT) => IllegalArgumentException
     *
     * @return
     */
    public String getParentContext() {
        if ( isRootContext() )
            throw new IllegalArgumentException("cannot get parent of root context");
        else
            return size() == 1 ? ROOT_CONTEXT : context.substring(1, size());
    }

    public boolean isRootContext() {
        return context.equals(ROOT_CONTEXT);
    }

    public int size() { return context.length(); }
}
