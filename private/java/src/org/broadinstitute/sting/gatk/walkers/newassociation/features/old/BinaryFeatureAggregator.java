package org.broadinstitute.sting.gatk.walkers.newassociation.features.old;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.walkers.newassociation.RFAArgumentCollection;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 5/4/11
 * Time: 12:52 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class BinaryFeatureAggregator {

    protected int nAberrant;
    protected int nNotAberrant;

    public BinaryFeatureAggregator(RFAArgumentCollection collection) {
        init(collection);
        nAberrant = 0;
        nNotAberrant = 0;
    }

    public void aggregate(SAMRecord record) {
        if ( featureDefined(record) ) {
            aggregate(extractFeature(record));
        }
    }

    protected void init(RFAArgumentCollection col) { }

    protected void aggregate(boolean feature) {
        if ( feature ) {
            nAberrant ++;
        } else {
            nNotAberrant ++;
        }
    }

    protected abstract boolean featureDefined(SAMRecord record);

    protected abstract boolean extractFeature(SAMRecord record);

    public int getnAberrant() { return nAberrant; }

    public int getnNotAberrant() { return nNotAberrant; }

    public int getnReads() { return nAberrant + nNotAberrant; }

    public double getMean() { return (nAberrant)/(1.0+getnReads()); }

    public double getVar() { return getnReads()*getMean()*(1-getMean()); }

    public double getUnbiasedVar() { return getVar(); }

    public Boolean parse(SAMRecord read) {
        if ( featureDefined(read) ) {
            return extractFeature(read);
        } else {
            return null;
        }
    }

    public String parseStr(SAMRecord read) {
        if ( featureDefined(read) ) {
            return Boolean.toString(extractFeature(read));
        } else {
            return "undefined";
        }
    }

}
