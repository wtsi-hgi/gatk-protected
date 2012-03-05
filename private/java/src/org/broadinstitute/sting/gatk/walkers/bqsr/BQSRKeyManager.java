package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.BitSetUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * This class provides all the functionality for the BitSet representation of the keys to the hash table of BQSR
 *
 * It also handles the event type "covariate" which is not exactly a covariate, but is added as a key to the hashmap. The Key Manager will
 * add the event type as a bitset to the end of the covariate bitset key. This way, it won't get int the way of masking the information
 * out of the key for the actual covariates, and having the covariates handle it. The key manager handles the event type.
 *
 * @author Mauricio Carneiro
 * @since 3/6/12
 */
public class BQSRKeyManager {
    private HashMap<Covariate, CovariateKeyInfo> covariateInfoMap;
    private CovariateKeyInfo[] covariatesInOrder;

    private int totalNumberOfBits;
    private int totalNumberOfCovariates;

    private int[] bitsBefore;

    /**
     * Initializes the KeyManager with the total number of covariates to use
     *
     * @param covariateList the ordered list of covariates
     */
    public BQSRKeyManager(List<Covariate> covariateList) {
        totalNumberOfCovariates = covariateList.size();
        bitsBefore = new int[totalNumberOfCovariates];
        totalNumberOfBits = bitsInEventType();                                     // The event type will use some number of bits and won't be added later as a covariate, so we have to add it now

        covariateInfoMap = new HashMap<Covariate, CovariateKeyInfo>(totalNumberOfCovariates * 2);
        covariatesInOrder = new CovariateKeyInfo[totalNumberOfCovariates];

        int order = 0;
        for (Covariate covariate : covariateList) {
            int nBits = covariate.numberOfBits();
            bitsBefore[order] = (order == 0) ? 0 : bitsBefore[order - 1] + covariatesInOrder[order - 1].nBits();
            CovariateKeyInfo info = new CovariateKeyInfo(order, maskFor(order, nBits), covariate);
            covariateInfoMap.put(covariate, info);
            covariatesInOrder[order] = info;
            totalNumberOfBits += nBits;
            order++;
        }
    }

    /**
     * Combines several keys in bitset representation into one bitset key and adds the event type
     *
     * @param keys      The keys in bitset representation for each covariate
     * @param eventType The type of event described by this keyset (e.g. mismatches, insertions, deletions)
     * @return one key in bitset representation that aggregates all keys
     */
    public BitSet bitSetFrom(BitSet[] keys, RecalDataManager.BaseRecalibrationType eventType) {
        BitSet fullKey = new BitSet(totalNumberOfBits);
        for (int order = 0; order < covariatesInOrder.length; order++) {
            for (int j = keys[order].nextSetBit(0); j >= 0; j = keys[order].nextSetBit(j + 1)) {
                if (j + bitsBefore[order] > totalNumberOfBits)
                    throw new ReviewedStingException("This is going beyond the total number of possible bits, something is wrong here!");
                fullKey.set(j + bitsBefore[order]);                                     // translate the bits set in the key to their corresponding position in the full key
            }
        }

        int lastBitIndexBeforeEventType = totalNumberOfBits - bitsInEventType();
        BitSet eventBitSet = BitSetUtils.bitSetFrom(eventType.index);                   // create a bitset with the event type
        for (int j = eventBitSet.nextSetBit(0); j >= 0; j = eventBitSet.nextSetBit(j + 1))
            fullKey.set(j + lastBitIndexBeforeEventType);                               // add the event type bits to the end of the full keyset
        return fullKey;
    }

    /**
     * Generates a key set of objects from a combined bitset key.
     *
     * Masks out each covariate independently and decodes their values (Object) into a keyset
     *
     * @param key the bitset representation of the keys
     * @return an object array with the values for each key
     */
    public Object[] keySetFrom(BitSet key) {
        Object[] keySet = new Object[totalNumberOfCovariates + 1];                      // We are creating one object per covariate plus the event type
        for (Map.Entry<Covariate, CovariateKeyInfo> covInfo : covariateInfoMap.entrySet()) {
            Covariate covariate = covInfo.getKey();
            CovariateKeyInfo covariateKeyInfo = covInfo.getValue();
            BitSet covariateValue = (BitSet) key.clone();                               // Make a copy of the key to extract the value for this covariate
            covariateValue.and(covariateKeyInfo.mask);                                  // Extract the bits for this particular covariate only
            covariateValue = chopKeyFor(covariateValue, covariate);                     // Chop it for just the covariate (to make it a single bitset the covariate can decode)
            keySet[covariateKeyInfo.order] = covariate.keyFromBitSet(covariateValue);   // Use the covariate's machinery to turn the BitSet into an object
        }
        keySet[totalNumberOfCovariates] = eventFrom(key);                               // Add the event to the key set
        return keySet;
    }

    /**
     * Translates a masked bitset into a bitset starting at 0
     *
     * @param key       the masked out bitset
     * @param covariate the covariate
     * @return a translated bitset starting at 0 for the covariate machinery to decode
     */
    private BitSet chopKeyFor(BitSet key, Covariate covariate) {
        CovariateKeyInfo covariateKeyInfo = covariateInfoMap.get(covariate);
        int order = covariateKeyInfo.order;
        int nBits = covariateKeyInfo.nBits();
        BitSet choppedKey = new BitSet(nBits);
        if (order > 0) {                                                                // No need to chop keys if this is the first covariate
            for (int i = key.nextSetBit(0); i >= 0; i = key.nextSetBit(i + 1))
                choppedKey.set(i - bitsBefore[order]);                                    // Set every bit translocated to the beginning of the BitSet
        }
        return choppedKey;                                                              // we don't need to clip the end of the bitset because it will all be zeroes anyway.
    }

    /**
     * Creates a mask for the requested covariate to extract the relevant bitset from a combined bitset key
     *
     * @param order the index of the covariate in the ordered covariate list
     * @param nBits the number of bits needed by the Covariate to represent its values in BitSet form
     * @return the bitset relevant to the covariate
     */
    private BitSet maskFor(int order, int nBits) {
        BitSet mask = new BitSet(totalNumberOfBits);
        int from = bitsBefore[order];
        int to = from + nBits;
        mask.set(from, to);
        return mask;
    }

    /**
     * Decodes the event type (enum) from the full bitset key
     *
     * @param fullKey the full key of all covariates + event type
     * @return the decoded event type.
     */
    private RecalDataManager.BaseRecalibrationType eventFrom(BitSet fullKey) {
        int nBits = bitsInEventType();
        BitSet eventKey = new BitSet(nBits);
        int firstBitIndex = totalNumberOfBits - nBits;
        for (int i = fullKey.nextSetBit(firstBitIndex); i >= 0; i = fullKey.nextSetBit(i + 1))
            eventKey.set(i - firstBitIndex);
        return RecalDataManager.BaseRecalibrationType.eventFrom((int) BitSetUtils.longFrom(eventKey));
    }

    private int bitsInEventType() {
        return BitSetUtils.numberOfBitsToRepresent(RecalDataManager.BaseRecalibrationType.values().length);
    }

    /**
     * Aggregate information for each Covariate
     */
    private class CovariateKeyInfo {
        public int order;
        public BitSet mask;
        public Covariate covariate;     // This allows reverse lookup of the Covariates in order

        private CovariateKeyInfo(int order, BitSet mask, Covariate covariate) {
            this.order = order;
            this.mask = mask;
            this.covariate = covariate;
        }

        /**
         * Just to make life a bit easier for the other members of the BQSRKeyManager class
         *
         * @return the number of bits used to represent this covariate's values
         */
        public int nBits() {
            return covariate.numberOfBits();
        }
    }
}
