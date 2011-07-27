package org.broadinstitute.sting.gatk.walkers.replication_validation;

import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.Collection;

/**
 * Short one line description of the walker.
 * <p/>
 * <p>
 * [Long description of the walker]
 * </p>
 * <p/>
 * <p/>
 * <h2>Input</h2>
 * <p>
 * [Description of the Input]
 * </p>
 * <p/>
 * <h2>Output</h2>
 * <p>
 * [Description of the Output]
 * </p>
 * <p/>
 * <h2>Examples</h2>
 * <pre>
 *    java
 *      -jar GenomeAnalysisTK.jar
 *      -T [walker name]
 *  </pre>
 *
 * @author Mauricio Carneiro
 * @since 7/27/11
 */
public class LaneParameters {
    public String name;
    public ReadBackedPileup lanePileup;
    public String referenceSampleName;
    public Collection<Byte> trueReferenceBases;
    public byte referenceSequenceBase;
    public byte minQualityScore;
    public byte maxQualityScore;
    public byte phredScaledPrior;
    public int maxAlleleCount;

    public LaneParameters(String name, ReadBackedPileup lanePileup, String referenceSampleName, Collection<Byte> trueReferenceBases, byte referenceSequenceBase, byte minQualityScore, byte maxQualityScore, byte phredScaledPrior, int maxAlleleCount) {
        this.name = name;
        this.lanePileup = lanePileup;
        this.referenceSampleName = referenceSampleName;
        this.trueReferenceBases = trueReferenceBases;
        this.referenceSequenceBase = referenceSequenceBase;
        this.minQualityScore = minQualityScore;
        this.maxQualityScore = maxQualityScore;
        this.phredScaledPrior = phredScaledPrior;
        this.maxAlleleCount = maxAlleleCount;
    }
}
