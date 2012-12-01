package org.broadinstitute.sting.gatk.walkers.na12878kb.core;

import com.mongodb.ReflectionDBObject;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;

import java.util.Date;
import java.util.LinkedList;
import java.util.List;

/**
 *
 * User: depristo
 * Date: 11/4/12
 * Time: 7:07 AM
 * To change this template use File | Settings | File Templates.
 */
public class CallSet extends ReflectionDBObject {
    private static final String CALLSET_VCF_HEADERLINE_KEY = "CallSet";

    private String name;
    private Date date;
    private boolean isReviewer;

    public CallSet() {}

    public CallSet(String name, Date date, boolean isReviewer) {
        this.name = name;
        this.date = date;
        this.isReviewer = isReviewer;
    }

    public static boolean isVCFHeaderLine(final VCFHeaderLine line) {
        return line.getKey().equals(CALLSET_VCF_HEADERLINE_KEY);
    }

    public CallSet( final VCFHeaderLine headerLine) {
        if ( ! isVCFHeaderLine(headerLine) )
            throw new IllegalArgumentException("VCFHeaderLine key must be " + CALLSET_VCF_HEADERLINE_KEY + " but saw " + headerLine.getKey());

        for ( final String keyValue : Utils.split(headerLine.getValue(), ";")) {
            final String[] parts = keyValue.split("=");
            if ( parts[0].equals("Name") ) setName(parts[1]);
            else if ( parts[0].equals("Date") ) setDate(new Date(Long.valueOf(parts[1])));
            else if ( parts[0].equals("IsReviewer") ) setReviewer(Boolean.valueOf(parts[1]));
            else
                throw new IllegalArgumentException("Unexpected key/value " + keyValue + " in " + headerLine);
        }
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public Date getDate() {
        return date;
    }

    public void setDate(Date date) {
        this.date = date;
    }

    public boolean isReviewer() {
        return isReviewer;
    }

    public void setReviewer(boolean reviewer) {
        isReviewer = reviewer;
    }

    public VCFHeaderLine asVCFHeaderLine() {
        final List<String> values = new LinkedList<String>();

        values.add("Name" + "=" + getName());
        values.add("Date" + "=" + getDate().getTime());
        values.add("IsReviewer" + "=" + isReviewer());

        return new VCFHeaderLine(CALLSET_VCF_HEADERLINE_KEY, Utils.join(";", values));
    }

    @Override
    public String toString() {
        return "CallSet{" +
                "name='" + name + '\'' +
                ", date=" + date +
                ", isReviewer=" + isReviewer +
                '}';
    }
}
