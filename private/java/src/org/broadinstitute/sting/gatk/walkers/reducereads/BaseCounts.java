package org.broadinstitute.sting.gatk.walkers.reducereads;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;

import java.util.EnumMap;
import java.util.Map;

/**
* Created by IntelliJ IDEA.
* User: depristo
* Date: 4/8/11
* Time: 2:55 PM
*/

final public class BaseCounts {
    public final static BaseIndex MAX_BASE_INDEX_WITH_NO_COUNTS = BaseIndex.A;
    public final static byte MAX_BASE_WITH_NO_COUNTS = MAX_BASE_INDEX_WITH_NO_COUNTS.getByte();

    private final Map<BaseIndex, Integer> counts;
    private final Map<BaseIndex, Long> sumQuals;

    public BaseCounts() {
        counts = new EnumMap<BaseIndex, Integer>(BaseIndex.class);
        sumQuals = new EnumMap<BaseIndex, Long>(BaseIndex.class);
        for ( BaseIndex i : BaseIndex.values() ) {
            counts.put(i,0);
            sumQuals.put(i, 0L);
        }
    }

    public static BaseCounts createWithCounts(int [] countsACGT) {
        BaseCounts baseCounts = new BaseCounts();
        baseCounts.counts.put(BaseIndex.A, countsACGT[0]);
        baseCounts.counts.put(BaseIndex.C, countsACGT[1]);
        baseCounts.counts.put(BaseIndex.G, countsACGT[2]);
        baseCounts.counts.put(BaseIndex.T, countsACGT[3]);
        return baseCounts;
    }

    @Requires("other != null")
    public void add(BaseCounts other) {
        for (BaseIndex i : BaseIndex.values() )
            counts.put(i, counts.get(i) + other.counts.get(i));
    }

    @Requires("other != null")
    public void sub(BaseCounts other) {
        for (BaseIndex i : BaseIndex.values() )
            counts.put(i, counts.get(i) - other.counts.get(i));
    }

    @Ensures("totalCount() == old(totalCount()) || totalCount() == old(totalCount()) + 1")
    public void incr(byte base) {
        BaseIndex i = BaseIndex.byteToBase(base);
        if ( i != null ) // no Ns
            counts.put(i, counts.get(i) + 1);
    }

    @Ensures("totalCount() == old(totalCount()) || totalCount() == old(totalCount()) + 1")
    public void incr(byte base, byte qual) {
        BaseIndex i = BaseIndex.byteToBase(base);
        if ( i != null ) { // no Ns
            counts.put(i, counts.get(i) + 1);
            sumQuals.put(i, sumQuals.get(i) + qual);
        }
    }

    @Ensures("result >= 0")
    public int getCount(byte base) {
        return counts.get(BaseIndex.byteToBase(base));
    }

    public byte baseWithMostCounts() {
        return maxBaseIndex().getByte();
    }

    @Ensures("result >= 0")
    public int countOfMostCommonBase() {
        return counts.get(maxBaseIndex());
    }

    @Ensures("result >= 0")
    public long sumQualsOfMostCommonBase() {
        return sumQuals.get(maxBaseIndex());
    }

    @Ensures("result >= 0")
    public byte averageQualsOfMostCommonBase() {
        return (byte) (sumQualsOfMostCommonBase() / countOfMostCommonBase());
    }


    @Ensures("result >= 0")
    public int totalCount() {
        int sum = 0;

        for ( int c : counts.values() ) {
            sum += c;
        }

        return sum;
    }

    /**
     * Given a base , it returns the proportional count of this base compared to all other bases
     * @param base
     * @return the proportion of this base over all other bases
     */
    @Ensures("result >=0")
    public double baseCountProportion(byte base) {
        return (double) counts.get(BaseIndex.byteToBase(base)) / totalCount();
    }

    /**
     * Given a base , it returns the proportional count of this base compared to all other bases
     * @param baseIndex
     * @return the proportion of this base over all other bases
     */
    @Ensures("result >=0")
    public double baseCountProportion(BaseIndex baseIndex) {
        return (double) counts.get(baseIndex) / totalCount();
    }


    @Ensures("result != null")
    public String toString() {
        StringBuilder b = new StringBuilder();
        for ( Map.Entry<BaseIndex,Integer> elt : counts.entrySet() ) {
            b.append(elt.toString()).append("=").append(elt.getValue()).append(",");
        }
        return b.toString();
    }

    @Ensures({"result != null", "totalCount() != 0 || result == MAX_BASE_INDEX_WITH_NO_COUNTS"})
    public BaseIndex maxBaseIndex() {
        BaseIndex maxI = MAX_BASE_INDEX_WITH_NO_COUNTS;
        for ( BaseIndex i : counts.keySet() )
            if ( counts.get(i) > counts.get(maxI) )
                maxI = i;
        return maxI;
    }
}
