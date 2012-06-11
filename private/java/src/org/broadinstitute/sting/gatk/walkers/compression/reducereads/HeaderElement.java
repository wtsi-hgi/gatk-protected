package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.LinkedList;

/**
 * The element that describes the header of the sliding window.
 *
 * Each site has a header element containing the counts of each base, it's reference based location and whether or
 * not the site has insertions (to it's right). It also contains information about the bases that have been filtered
 * out due to mapping or base quality.
 */
public class HeaderElement {
    private BaseCounts consensusBaseCounts;      // How many A,C,G,T (and D's) are in this site.
    private BaseCounts filteredBaseCounts;       // How many A,C,G,T (and D's) were filtered out in this site.
    private int insertionsToTheRight;            // How many reads in this site had insertions to the immediate right
    private int location;                        // Genome location of this site (the sliding window knows which contig we're at
    private LinkedList<Integer> mappingQuality;  // keeps the mapping quality of each read that contributed to this element (site)

    public int getLocation() {
        return location;
    }

    public BaseCounts getFilteredBaseCounts() {
        return filteredBaseCounts;
    }

    public BaseCounts getConsensusBaseCounts() {
        return consensusBaseCounts;
    }

    /**
     * Creates a new HeaderElement with the following default values: - empty consensusBaseCounts - empty
     * filteredBaseCounts - 0 insertions to the right - empty mappingQuality list
     *
     * @param location the reference location for the new element
     */
    public HeaderElement(int location) {
        this(new BaseCounts(), new BaseCounts(), 0, location, new LinkedList<Integer>());
    }

    /**
     * Creates a new HeaderElement with all given parameters
     *
     * @param consensusBaseCounts  the BaseCounts object for the running consensus synthetic read
     * @param filteredBaseCounts   the BaseCounts object for the filtered data synthetic read
     * @param insertionsToTheRight number of insertions to the right of this HeaderElement
     * @param location             the reference location of this reference element
     * @param mappingQuality       the list of mapping quality values of all reads that contributed to this
     *                             HeaderElement
     */
    public HeaderElement(BaseCounts consensusBaseCounts, BaseCounts filteredBaseCounts, int insertionsToTheRight, int location, LinkedList<Integer> mappingQuality) {
        this.consensusBaseCounts = consensusBaseCounts;
        this.filteredBaseCounts = filteredBaseCounts;
        this.insertionsToTheRight = insertionsToTheRight;
        this.location = location;
        this.mappingQuality = mappingQuality;
    }

    /**
     * Whether or not the site represented by this HeaderElement is variant according to the definitions of variant
     * by insertion, deletion and mismatches.
     *
     * @return true if site is variant by any definition. False otherwise.
     */
    public boolean isVariant(double minVariantProportion, double minIndelProportion) {
        return hasConsensusData() && (isVariantFromInsertions(minIndelProportion) || isVariantFromMismatches(minVariantProportion) || isVariantFromDeletions(minIndelProportion));
    }

    /**
     * Adds a new base to the HeaderElement updating all counts accordingly
     *
     * @param base           the base to add
     * @param baseQual           the base quality
     * @param baseMappingQuality the mapping quality of the read this base belongs to
     */
    public void addBase(byte base, byte baseQual, int baseMappingQuality, int minBaseQual, int minMappingQual) {
        if (basePassesFilters(baseQual, minBaseQual, baseMappingQuality, minMappingQual))
            consensusBaseCounts.incr(base, baseQual);                                                                   // If the base passes filters, it is included in the consensus base counts
        else
            filteredBaseCounts.incr(base, baseQual);                                                                    // If the base fails filters, it is included with the filtered data base counts

        this.mappingQuality.add(baseMappingQuality);                                                                    // Filtered or not, the RMS mapping quality includes all bases in this site
    }

    public void removeBase(byte base, byte baseQual, int baseMappingQuality, int minBaseQual, int minMappingQual) {
        if (basePassesFilters(baseQual, minBaseQual, baseMappingQuality, minMappingQual))
            consensusBaseCounts.decr(base, baseQual);                                                                   // If the base passes filters, it is included in the consensus base counts
        else
            filteredBaseCounts.decr(base, baseQual);                                                                    // If the base fails filters, it is included with the filtered data base counts

        this.mappingQuality.remove((Integer) baseMappingQuality);                                                       // Filtered or not, the RMS mapping quality includes all bases in this site

    }
    /**
     * Adds an insertions to the right of the HeaderElement and updates all counts accordingly. All insertions
     * should be added to the right of the element.
     */
    public void addInsertionToTheRight() {
        insertionsToTheRight++;
    }

    /**
     * Does this HeaderElement contain consensus data?
     *
     * @return whether or not this HeaderElement contains consensus data
     */
    public boolean hasConsensusData() {
        return consensusBaseCounts.totalCount() > 0;
    }

    /**
     * Does this HeaderElement contain filtered data?
     *
     * @return whether or not this HeaderElement contains filtered data
     */
    public boolean hasFilteredData() {
        return filteredBaseCounts.totalCount() > 0;
    }

    /**
     * A HeaderElement is empty if it has no consensus or filtered data
     *
     * @return whether or not this HeaderElement has no data
     */
    public boolean isEmpty() {
        return (!hasFilteredData() && !hasConsensusData());
    }

    /**
     * The RMS of the mapping qualities of all reads that contributed to this HeaderElement
     *
     * @return the RMS of the mapping qualities of all reads that contributed to this HeaderElement
     */
    public double getRMS() {
        return MathUtils.rms(mappingQuality);
    }

    /**
     * Whether or not the HeaderElement is variant due to excess insertions
     *
     * @return whether or not the HeaderElement is variant due to excess insertions
     */
    public boolean isVariantFromInsertions(double minIndelProportion) {
        int numberOfBases = consensusBaseCounts.totalCount();
        if (numberOfBases == 0 && insertionsToTheRight > 0)
            return true;  // we only have insertions
        else if (numberOfBases == 0)
            return false; // we don't have anything

        // if we have bases and insertions, check the ratio
        return ((double) insertionsToTheRight / numberOfBases) > minIndelProportion;
    }

    /**
     * Whether or not the HeaderElement is variant due to excess deletions
     *
     * @return whether or not the HeaderElement is variant due to excess insertions
     */
    public boolean isVariantFromDeletions(double minIndelProportion) {
        return consensusBaseCounts.baseIndexWithMostCounts() == BaseIndex.D || consensusBaseCounts.baseCountProportion(BaseIndex.D) > minIndelProportion;
    }

    /**
     * Whether or not the HeaderElement is variant due to excess mismatches
     *
     * @return whether or not the HeaderElement is variant due to excess insertions
     */
    public boolean isVariantFromMismatches(double minVariantProportion) {
        BaseIndex mostCommon = consensusBaseCounts.baseIndexWithMostCountsWithoutIndels();
        double mostCommonProportion = consensusBaseCounts.baseCountProportionWithoutIndels(mostCommon);
        return mostCommonProportion != 0.0 && mostCommonProportion < (1 - minVariantProportion);
    }

    private boolean basePassesFilters(byte baseQual, int minBaseQual, int baseMappingQuality, int minMappingQual) {
        return baseQual >= minBaseQual && baseMappingQuality >= minMappingQual;
    }

    /**
     * removes an insertion from this element (if you removed a read that had an insertion)
     */
    public void removeInsertionToTheRight() {
        this.insertionsToTheRight--;
        if (insertionsToTheRight < 0)
            throw new ReviewedStingException("Removed too many insertions, header is now negative!");
    }
}