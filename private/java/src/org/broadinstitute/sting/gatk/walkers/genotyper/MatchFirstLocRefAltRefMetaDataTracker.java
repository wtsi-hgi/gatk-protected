package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

// Use this class to store a single vc to be retrieved and genotyped:
public class MatchFirstLocRefAltRefMetaDataTracker extends RefMetaDataTracker {
    public MatchFirstLocRefAltRefMetaDataTracker(RefMetaDataTracker tracker, VariantContext vc) {
        super();
        this.tracker = tracker;
        this.vc = vc;
    }

    public <T extends Feature> List<T> getValues(final RodBinding<T> type, final GenomeLoc onlyAtThisLoc) {
        List<T> l = new LinkedList<T>();

        for (T t : tracker.getValues(type, onlyAtThisLoc)) {
            VariantContext tVC = (VariantContext) t;
            if (tVC.getChr().equals(vc.getChr()) && tVC.getStart() == vc.getStart() && tVC.getEnd() == vc.getEnd() && tVC.getReference().equals(vc.getReference()) && tVC.getAlternateAlleles().equals(vc.getAlternateAlleles())) {
                l.add((T) tVC);
                break;
            }
        }

        return l;
    }

    public <T extends Feature> List<T> getValues(final Collection<RodBinding<T>> rodBindings, final GenomeLoc onlyAtThisLoc) {
        List<T> l = new LinkedList<T>();

        for (RodBinding<T> type : rodBindings) {
            l.addAll(getValues(type, onlyAtThisLoc));
            if (!l.isEmpty())
                break;
        }

        return l;
    }

    private RefMetaDataTracker tracker;
    private VariantContext vc;
}
