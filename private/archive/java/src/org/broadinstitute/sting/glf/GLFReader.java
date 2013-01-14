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
import net.sf.samtools.util.BlockCompressedInputStream;
import net.sf.samtools.util.RuntimeEOFException;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.genotype.LikelihoodObject;

import java.io.DataInputStream;
import java.io.File;
import java.io.IOException;
import java.util.Iterator;

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

/** an object for reading in GLF files */
// TODO -- DELETE ME GLF
public class GLFReader implements Iterator<GLFRecord> {

    // our next record
    private GLFRecord nextRecord = null;

    // the glf magic number, which identifies a properly formatted GLF file
    public static final short[] glfMagic = {'G', 'L', 'F', '\3'};

    // our input codec
    private final BinaryCodec inputBinaryCodec;

    // our header string
    private String headerStr;

    // our reference name
    private String referenceName;

    // reference length
    private int referenceLength;

    // the current location, keeping track of the offsets
    private long currentLocation = 1;

    // we have this variable becuase there is no eof for glf's
    private int lastRecordType = -1;

    private File myFile;

    /**
     * create a glf reader
     *
     * @param readFrom the file to read from
     */
    public GLFReader(File readFrom) {
        myFile = readFrom;

        try {
            inputBinaryCodec = new BinaryCodec(new DataInputStream(new BlockCompressedInputStream(readFrom)));
        } catch (IOException e) {
            throw new UserException.CouldNotReadInputFile(myFile, e);
        }

        inputBinaryCodec.setInputFileName(readFrom.getName());

        // first verify that it's a valid GLF
        for (short s : glfMagic) {
            if (inputBinaryCodec.readUByte() != s)
                throw new UserException.MalformedFile(myFile, "Verification of GLF format failed: magic string doesn't match)");
        }

        // get the header string
        headerStr = inputBinaryCodec.readLengthAndString(false);

        if (advanceContig()) {
            // setup the next record
            next();
        }

    }

    /**
     * read in a single point call
     *
     * @param refBase          the reference base
     * @param inputBinaryCodec the binary codec
     *
     * @return a single point call object
     */
    private GLFSingleCall generateSPC(char refBase, BinaryCodec inputBinaryCodec) {
        int offset = (int) inputBinaryCodec.readUInt();
        long depth = inputBinaryCodec.readUInt();
        short min_lk = (short) ((depth & 0x00000000ff000000) >> 24);
        int readDepth = (int) (depth & 0x0000000000ffffff);
        short rmsMapping = inputBinaryCodec.readUByte();
        double[] lkValues = new double[LikelihoodObject.GENOTYPE.values().length];
        for (int x = 0; x < LikelihoodObject.GENOTYPE.values().length; x++) {
            lkValues[x] = ((double)inputBinaryCodec.readUByte() / GLFRecord.LIKELIHOOD_SCALE_FACTOR + (double)min_lk);
        }
        return new GLFSingleCall(referenceName, refBase, (int)(offset+currentLocation), readDepth, rmsMapping, lkValues);
    }

    /**
     * read in a variable length call, and generate a VLC object from the data
     *
     * @param refBase          the reference base
     * @param inputBinaryCodec the input codex
     *
     * @return a GLFVariableLengthCall object
     */
    private GLFVariableLengthCall generateVLC(char refBase, BinaryCodec inputBinaryCodec) {
        int offset = (int) inputBinaryCodec.readUInt();
        int depth = (int) inputBinaryCodec.readUInt();
        short min_lk = (short) ((depth & 0x00000000ff000000) >> 24);
        int readDepth = (depth & 0x0000000000ffffff);
        short rmsMapping = inputBinaryCodec.readUByte();
        short lkHom1 = inputBinaryCodec.readUByte();
        short lkHom2 = inputBinaryCodec.readUByte();
        short lkHet = inputBinaryCodec.readUByte();
        int indelLen1 = (int) inputBinaryCodec.readShort();
        int indelLen2 = (int) inputBinaryCodec.readShort();

        int readCnt = Math.abs(indelLen1);
        short indelSeq1[] = new short[readCnt];
        for (int x = 0; x < readCnt; x++) {
            indelSeq1[x] = inputBinaryCodec.readUByte();
        }
        readCnt = Math.abs(indelLen2);
        short indelSeq2[] = new short[readCnt];
        for (int x = 0; x < readCnt; x++) {
            indelSeq2[x] = inputBinaryCodec.readUByte();
        }
        return new GLFVariableLengthCall(referenceName, refBase, offset+currentLocation, readDepth, rmsMapping, lkHom1, lkHom2, lkHet, indelLen1, indelSeq1, indelLen2, indelSeq2);
    }

    public boolean hasNext() {
        return (nextRecord != null);
    }

    public GLFRecord next() {
        GLFRecord ret = nextRecord;
        short firstBase = protectedByteReadForFile();
        if (firstBase == -1) return ret;

        // parse out the record type and reference base
        byte recordType = (byte) ((firstBase & 0x0f0) >> 4);
        char refBase = (char) (firstBase & 0x000f);
        lastRecordType = recordType;

        if (recordType == 1) {
            nextRecord = generateSPC(refBase, inputBinaryCodec);
        } else if (recordType == 2) {
            nextRecord = generateVLC(refBase, inputBinaryCodec);
        } else if (recordType == 0) {
            if (advanceContig()) {
                return next();
            }
            //nextRecord = null;
        } else {
            throw new UserException.MalformedFile(myFile, "Unknown GLF record type (type = " + recordType + ")");
        }
        if (nextRecord != null) currentLocation = nextRecord.getPosition();
        return ret;
    }

    /**
     * read a short, and if we see an exception only throw it if it's unexpected (not after a zero)
     * @return a short
     */
    private short protectedByteReadForFile() {
        short st = -1;
        try {
            st = inputBinaryCodec.readUByte();
        } catch (RuntimeEOFException exp) {
            nextRecord = null;
            if (lastRecordType != 0) {
                throw exp; // if the last record was a zero, this is an ok condition.  Otherwise throw an exception
            }
        }
        return st;
    }

    /**
     * advance to the next contig
     *
     * @return true if we could advance
     */
    private boolean advanceContig() {
        // try to read the next sequence record
        try {
            // get the reference name
            referenceName = inputBinaryCodec.readLengthAndString(true);

            // get the reference length - this may be a problem storing an unsigned int into a signed int.  but screw it.
            referenceLength = (int) inputBinaryCodec.readUInt();
            //System.err.println(referenceName.length());
            currentLocation = 1;
            return true;
        } catch (RuntimeException e) {
            if (lastRecordType != 0) {
                throw e; // if the last record was a zero, this is an ok condition.  Otherwise throw an exception
            }
            nextRecord = null;
        }
        return false;
    }

    public void remove() {
        throw new ReviewedStingException("GLFReader doesn't support remove()");
    }

    public void close() {
        inputBinaryCodec.close();
    }
    
    public String getHeaderStr() {
        return headerStr;
    }

}
