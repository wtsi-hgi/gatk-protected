package org.broadinstitute.sting.gatk.walkers.poolcaller;

import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.Collection;

/**
 * A support class to facilitate future addition/removal of parameters to the Reference Sample class
 *
 * @author Mauricio Carneiro
 * @since 7/27/11
 */
public class ReferenceSampleParameters {
    public String name;
    public ReadBackedPileup pileup;
    public Collection<Byte> trueBases;

    public ReferenceSampleParameters(String name, ReadBackedPileup pileup, Collection<Byte> trueBases) {
        this.name = name;
        this.pileup = pileup;
        this.trueBases = trueBases;
    }
}
