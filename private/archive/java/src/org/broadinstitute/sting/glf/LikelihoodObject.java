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

package org.broadinstitute.sting.utils.genotype;

import edu.mit.broad.picard.genotype.DiploidGenotype;
import edu.mit.broad.picard.genotype.geli.GenotypeLikelihoods;
import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.HashMap;


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
 *         Class LikelyhoodObject
 *         <p/>
 *         An object used to store likelyhood information for genotypes.  Genotype
 *         likelihoods are assumed to be infinite (negitive log likelihood), unless set.
 *         This allows the consumer to make an empty LikelihoodObject, and just set
 *         those values which have associated likelihood values.
 */

// TODO -- DELETE ME GLF
public class LikelihoodObject {


    // our possible genotypes, in order according to GLFv3
    public enum GENOTYPE {
        AA, AC, AG, AT, CC, CG, CT, GG, GT, TT
    }

    // our pileup of bases
    //final private String basePileup;

    // possible types of likihoods to store

    public enum LIKELIHOOD_TYPE {
        NEGATIVE_LOG, LOG, RAW;
    }

    // our liklihood storage type
    protected LIKELIHOOD_TYPE mLikelihoodType = LIKELIHOOD_TYPE.NEGATIVE_LOG;

    // default the bestGenotype likelihoods to the allele AA
    protected GENOTYPE bestGenotype = GENOTYPE.AA;

    // how many genotypes we're storing
    public static final int genoTypeCount = GENOTYPE.values().length;

    // the associated negitive log likelihoods for each genotype
    protected final HashMap<GENOTYPE, Double> likelihoods = new HashMap<GENOTYPE, Double>();

    /** create a blank likelihood object */
    public LikelihoodObject() {
        for (GENOTYPE type : GENOTYPE.values()) {
            likelihoods.put(type, Double.MAX_VALUE);
        }
    }

    /**
     * create a likelihood object, given a picard style GenotypeLikelihoods object.  The
     * GenotypeLikelihoods stores likelihoods in log likelihood format, and we want them in
     * negitive log likelihood
     *
     * @param lk the likelihood object
     */
    public LikelihoodObject(GenotypeLikelihoods lk) {
        mLikelihoodType = LIKELIHOOD_TYPE.LOG;
        Double minValue = Double.MAX_VALUE;
        for (GENOTYPE type : GENOTYPE.values()) {
            byte[] bases = new byte[2];
            bases[0] = (byte) type.toString().charAt(0);
            bases[1] = (byte) type.toString().charAt(1);
            double val = -1.0d * lk.getLikelihood(DiploidGenotype.fromBases(bases));
            likelihoods.put(type, val);
            if (val < minValue) {
                bestGenotype = type;
            }
        }
    }

    /**
     * create a likelyhood object, given an array of genotype scores in GLFv3 ordering
     *
     * @param values an array of int's from 0 to 255, representing the negitive log likelihoods.
     * @param type   the likelihood storage type
     */
    public LikelihoodObject(double[] values, LIKELIHOOD_TYPE type) {
        mLikelihoodType = type;
        if (values.length != GENOTYPE.values().length) {
            throw new IllegalArgumentException("invalid array passed to LikelihoodObject, should be size " + GENOTYPE.values().length);
        }
        findBestLikelihood(values);
    }

    /**
     * find the best likelihood
     * @param values
     */
    private void findBestLikelihood(double[] values) {
        int index = 0;
        double lowestScore = Double.MAX_VALUE;
        for (GENOTYPE t : GENOTYPE.values()) {
            likelihoods.put(t, values[index]);
            if (values[index] < lowestScore) {
                lowestScore = values[index];
                bestGenotype = t;
            }
            ++index;
        }
    }

    /**
     * set the likelihood, given it's probability and the genotype
     *
     * @param type the genotype
     * @param lh   the likelihood as a double
     */
    public void setLikelihood(GENOTYPE type, double lh) {
        likelihoods.put(type, lh);
        if (lh < likelihoods.get(this.bestGenotype)) {
            this.bestGenotype = type;
        }
    }

    /**
     * find the minimum likelihood value stored in the set.  This represents the most likely genotype,
     * since genotypes are represented as negitive log likeihoods
     *
     * @return the min value
     */
    public double getBestLikelihood() {
        return likelihoods.get(this.bestGenotype);
    }

    /**
     * return a byte array representation of the likelihood object, in GLFv3 specified order.
     * The return type is short[] instead of byte[], since signed bytes only store -127 to 127,
     * not the 255 range we need.
     *
     * @return a byte array of the genotype values
     */
    public short[] toByteArray() {
        short ret[] = new short[GENOTYPE.values().length];
        int index = 0;
        for (GENOTYPE type : GENOTYPE.values()) {
            ret[index] = (likelihoods.get(type).intValue() > 254) ? 255 : (short) likelihoods.get(type).intValue();
            ++index;
        }
        return ret;
    }

    /**
     * create a float array of our genotype values, in order specified in the GENOTYPE enum (currently the GLF and
     * geli ordering).
     *
     * @return a float array containing our genotype likelihoods, as negitive log likelihoods
     */
    public double[] toDoubleArray() {
        // make an array of floats
        double[] ft = new double[10];
        int index = 0;
        for (GENOTYPE T : GENOTYPE.values()) {
            ft[index] = this.likelihoods.get(T).doubleValue();
            index++;
        }
        return ft;
    }

    /**
     * convert this object, with aditional information, to a GenotypeLikelihoods object.  This involves determining
     * what our underlying storage type is, and coverting our values to the appropriate (log likelihood) format.
     *
     * @return a GenotypeLikelihoods object representing our data
     */
    public GenotypeLikelihoods convertToGenotypeLikelihoods(SAMFileHeader samHeader, int seqIndex, int seqPosition, byte refBase) {
        double[] ft = toDoubleArray();
        float[] db = new float[ft.length];
        int index = 0;
        if (this.mLikelihoodType == LIKELIHOOD_TYPE.NEGATIVE_LOG) {
            for (; index < ft.length; index++) {
                db[index] = ((float) ft[index] * -1.0f);
            }
        } else if (this.mLikelihoodType == LIKELIHOOD_TYPE.RAW) {
            for (; index < ft.length; index++) {
                db[index] = (float) Math.log(ft[index]);
            }
        } else {
            for (int x = 0; x < ft.length; x++)
                db[x] = (float)ft[x];
        }
        return new GenotypeLikelihoods(samHeader, seqIndex, seqPosition, refBase, db);
    }

    /**
     * getter for the likelihood type
     *
     * @return our likelihood storage type
     */
    public LIKELIHOOD_TYPE getLikelihoodType() {
        return mLikelihoodType;
    }


    /**
     * validate a genotype score
     *
     * @param score the score to validate
     */
    public void validateScore(double score) {
        int x = 0;
        switch (mLikelihoodType) {
            case NEGATIVE_LOG:
                if (score < 0)
                    throw new ReviewedStingException("Likelikhood score of " + score + " is invalid, for NEGATIVE_LOG it must be greater than or equal to 0");
                break;
            case LOG:
                if (score > 0)
                    throw new ReviewedStingException("Likelikhood score of " + score + " is invalid, for LOG it must be less than or equal to 0");
                break;
            case RAW:
                if (score < 0 || score > 1)
                    throw new ReviewedStingException("Likelikhood score of " + score + " is invalid, for RAW it must be [0,1]");
                break;
        }
    }


    /**
     * set our likelihood storage type, and adjust our current likelihood values to reflect
     * the new setting.
     *
     * @param likelihood the type to set the values to.
     */
    public void setLikelihoodType(LIKELIHOOD_TYPE likelihood) {
        if (likelihood == mLikelihoodType)
            return;
        if (mLikelihoodType == LIKELIHOOD_TYPE.RAW) {
            double mult = 1.0;
            if (likelihood == LIKELIHOOD_TYPE.NEGATIVE_LOG) {
                mult = -1.0;
            }
            // one of us in log, the other negitive log, it doesn't matter which
            for (GENOTYPE g : likelihoods.keySet()) {
                likelihoods.put(g, -1.0 * Math.log(likelihoods.get(g)));
            }
        } else if (likelihood == LIKELIHOOD_TYPE.RAW) {
            double mult = 1.0;
            if (mLikelihoodType == LIKELIHOOD_TYPE.NEGATIVE_LOG) {
                mult = -1.0;
            }
            // one of us in log, the other negitive log, it doesn't matter which
            for (GENOTYPE g : likelihoods.keySet()) {
                likelihoods.put(g, Math.pow(likelihoods.get(g) * mult, 10));
            }
        } else {
            // one of us in log, the other negitive log, it doesn't matter which
            for (GENOTYPE g : likelihoods.keySet()) {
                likelihoods.put(g, -1.0 * likelihoods.get(g));
            }
        }
        this.mLikelihoodType = likelihood;
    }
}


