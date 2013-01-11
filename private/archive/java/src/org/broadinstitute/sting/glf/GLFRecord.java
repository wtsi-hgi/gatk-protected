/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.utils.genotype.glf;

import net.sf.samtools.util.BinaryCodec;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;


/*
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @author aaron
 *         <p/>
 *         Class RecordType
 *         <p/>
 *         The base record type for all GLF entries. Each record has a number of fields
 *         common to the record set.  This is also the source of the REF_BASE enumeration,
 *         which represents the accepted FASTA nucleotide symbols and their assocated GLF
 *         field values.
 */
// TODO -- DELETE ME GLF
public abstract class GLFRecord {
    public final static double LIKELIHOOD_SCALE_FACTOR = 10;


    // fields common to all records
    protected String contig;
    protected REF_BASE refBase;
    protected long position = 1;
    protected int readDepth = 0;
    protected short rmsMapQ = 0;

    /** the reference base enumeration, with their short (type) values for GLF */
    public enum REF_BASE {
        X((short) 0x00),
        A((short) 0x01),
        C((short) 0x02),
        M((short) 0x03),
        G((short) 0x04),
        R((short) 0x05),
        S((short) 0x06),
        V((short) 0x07),
        T((short) 0x08),
        W((short) 0x09),
        Y((short) 0x0A),
        H((short) 0x0B),
        K((short) 0x0C),
        D((short) 0x0D),
        B((short) 0x0E),
        N((short) 0x0F);

        private final short fieldValue;

        /**
         * private constructor, used by the enum class to makes each enum value
         *
         * @param value the short values specified in the enum listing
         */
        REF_BASE(short value) {
            fieldValue = value;
        }

        /**
         * return the character representation
         *
         * @return the char for the reference base
         */
        public char toChar() {
            return this.toString().charAt(0);
        }

        /**
         * static method from returning a REF_BASE given the character representation
         *
         * @param value the character representation of a REF_BASE
         *
         * @return the corresponding REF_BASE
         * @throws IllegalArgumentException if the value passed can't be converted
         */
        public static REF_BASE toBase(char value) {
            // for the case where they're passing in the enumated value
            if (value <= 0x0F && value >= 0) {
                return REF_BASE.values()[value];
            }
            String str = String.valueOf(value).toUpperCase();
            for (int x = 0; x < REF_BASE.values().length; x++) {
                if (REF_BASE.values()[x].toString().equals(str)) {
                    return REF_BASE.values()[x];
                }
            }
            throw new IllegalArgumentException("Counldn't find matching reference base for " + str);
        }

        /** @return the hex value of the given REF_BASE */
        public short getBaseHexValue() {
            return fieldValue;
        }
    }

    /** the record type enum, which enumerates the different records we can have in a GLF */
    public enum RECORD_TYPE {
        SINGLE((short) 1),
        VARIABLE((short) 2);

        private final short fieldValue;   // a place to store the type

        RECORD_TYPE(short value) {
            fieldValue = value;
        }

        public short getReadTypeValue() {
            return fieldValue;
        }
    }


    /**
     * Constructor, given the base a character reference base
     *
     * @param contig            the contig string
     * @param base              the reference base in the reference
     * @param position          the distance from the beginning of the reference seq
     * @param readDepth         the read depth at this position
     * @param rmsMapQ           the root mean square of the mapping quality
     */
    public GLFRecord(String contig, char base, long position, int readDepth, short rmsMapQ) {
        REF_BASE newBase = REF_BASE.toBase(base);
        validateInput(contig, newBase, position, readDepth, rmsMapQ);
    }

    /**
     * Constructor, given the base a REF_BASE
     *
     * @param contig            the contig string
     * @param base              the reference base in the reference
     * @param position          the distance from the beginning of the reference seq
     * @param readDepth         the read depth at this position
     * @param rmsMapQ           the root mean square of the mapping quality
     */
    GLFRecord(String contig, REF_BASE base, long position, int readDepth, short rmsMapQ) {
        validateInput(contig, base, position, readDepth, rmsMapQ);
    }

    /**
     * validate the input during construction, and store valid values
     *
     * @param chromosome        the reference contig, as a String
     * @param base              the reference base in the reference, as a REF_BASE
     * @param position          the distance from the beginning of the reference seq
     * @param readDepth         the read depth at this position
     * @param rmsMapQ           the root mean square of the mapping quality
     */
    private void validateInput(String chromosome, REF_BASE base, long position, int readDepth, short rmsMapQ) {
        // add any validation to the contig string here
        this.contig = chromosome;

        this.refBase = base;

        if (position > 4294967295L || position < 0) {
            throw new IllegalArgumentException("Position is out of bounds (0 to 0xffffffff) value passed = " + position);
        }
        this.position = position;

//        if (minimumLikelihood > 255 || minimumLikelihood < 0) {
//            throw new IllegalArgumentException("minimumLikelihood is out of bounds (0 to 0xffffffff) value passed = " + minimumLikelihood);
//        }
//        this.minimumLikelihood = GLFRecord.toCappedShort(minimumLikelihood);

        if (readDepth > 16777215 || readDepth < 0) {
            throw new IllegalArgumentException("readDepth is out of bounds (0 to 0xffffff) value passed = " + readDepth);
        }
        this.readDepth = readDepth;

        if (rmsMapQ > 255 || rmsMapQ < 0) {
            throw new IllegalArgumentException("rmsMapQ is out of bounds (0 to 0xff) value passed = " + rmsMapQ);
        }
        this.rmsMapQ = rmsMapQ;
    }

    /**
     * write the this record to a binary codec output.
     *
     * @param out the binary codec to write to
     * @param lastRecord the record to write
     */
    void write(BinaryCodec out, GLFRecord lastRecord) {
        long offset;
        if (lastRecord != null && lastRecord.getContig().equals(this.getContig()))
            offset = this.position - lastRecord.getPosition();
        else
            offset = this.position - 1; // we start at one, we need to subtract that off
        out.writeUByte((short) (this.getRecordType().getReadTypeValue() << 4 | (refBase.getBaseHexValue() & 0x0f)));
        out.writeUInt(((Long) (offset)).intValue()); // we have to subtract one, we're an offset
        long write = ((long) (readDepth & 0xffffff) | (long) (this.getMinimumLikelihood() & 0xff) << 24);
        out.writeUInt(write);
        out.writeUByte(rmsMapQ);
    }

    /**
     * get the record type
     *
     * @return the record type enumeration
     */
    public abstract RECORD_TYPE getRecordType();

    /**
     * Return the size of this record in bytes.
     *
     * @return the size of this record type, in bytes
     */
    public int getByteSize() {
        return 10; // the record type field (1), offset (4), the min depth field (4), and the rms mapping (1)
    }

    /**
     * convert a double to a byte, capping it at the maximum value of 255
     *
     * @param d a double value
     *
     * @return a byte, capped at
     */
    protected static short toCappedShort(double d) {
        return (d > 255.0) ? (short) 255 : (short) Math.round(d);
    }

    /**
     * find the minimum value in a set of doubles
     *
     * @param vals the array of values
     *
     * @return the minimum value
     */
    protected static double findMin(double vals[]) {
        if (vals.length < 1) throw new ReviewedStingException("findMin: an array of size < 1 was passed in");

        double min = vals[0];
        for (double d : vals)
            if (d < min) min = d;

        return min;
    }

    public REF_BASE getRefBase() {
        return refBase;
    }

    public long getPosition() {
        return position;
    }

    public short getMinimumLikelihood() {
        return calculateMinLikelihood();
    }

    public int getReadDepth() {
        return readDepth;
    }

    public short getRmsMapQ() {
        return rmsMapQ;
    }

    public String getContig() {
        return this.contig;
    }

    /**
     * this method had to be abstracted so that the underlying records could set the minimum likelihood (ML) in the event
     * that the ML is above 255.  In this case the records need to scale the value appropriately, and warn the users.
     * @return a short of the minimum likelihood.
     */
    protected abstract short calculateMinLikelihood();

    public boolean equals(GLFRecord rec) {
        return (rec != null) &&
                contig.equals(rec.getContig()) &&
                (refBase == rec.getRefBase()) &&
                (position == rec.getPosition()) &&
                (readDepth == rec.getReadDepth()) &&
                (rmsMapQ == rec.getRmsMapQ());
    }


}

