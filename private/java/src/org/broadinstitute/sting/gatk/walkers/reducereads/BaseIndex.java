package org.broadinstitute.sting.gatk.walkers.reducereads;

/**
 * Simple byte / base index conversions
 *
 * @author carneiro
 * @since 8/26/11
 */
public enum BaseIndex {
    A  ( 'A', 0 ),
    C  ( 'C', 1 ),
    G  ( 'G', 2 ),
    T  ( 'T', 3 ),
    D  ( 'D', 4 ),
    I  ( 'I', 5 ), // insertion to the right of the base
    EQ ( '=', 6 );

    final byte b;
    final int index;

    private BaseIndex(char base, int index) {
        this.b = (byte)base;
        this.index = index;
    }

    public byte getByte() { return b; }

    public static BaseIndex byteToBase(final byte base) {
        switch (base) {
            case 'A':
            case 'a':
                return A;
            case 'C':
            case 'c':
                return C;
            case 'G':
            case 'g':
                return G;
            case 'T':
            case 't':
                return T;
            case 'D':
            case 'd':
            case '-':
                return D;
            case 'I':
            case 'i':
                return I;
            case '=':
                return EQ;
            default: return null;
        }
    }
}
