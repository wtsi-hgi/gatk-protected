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

package org.broadinstitute.sting.utils.recalibration;

import org.broadinstitute.sting.gatk.walkers.bqsr.*;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.collections.NestedIntegerArray;
import org.broadinstitute.sting.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Unit tests for on-the-fly recalibration.
 *
 * @author Mauricio Carneiro
 * @since 3/16/12
 */
public class BaseRecalibrationUnitTest {

    private RecalDataManager dataManager;

    private ReadGroupCovariate rgCovariate;
    private QualityScoreCovariate qsCovariate;
    private ContextCovariate cxCovariate;
    private CycleCovariate cyCovariate;

    private GATKSAMRecord read = ReadUtils.createRandomRead(10000);
    private BaseRecalibration baseRecalibration;
    private ReadCovariates readCovariates;


    @BeforeClass
    public void init() {
        GATKSAMReadGroupRecord rg = new GATKSAMReadGroupRecord("rg");
        rg.setPlatform("illumina");
        read.setReadGroup(rg);

        byte[] quals = new byte[read.getReadLength()];
        for (int i = 0; i < read.getReadLength(); i++)
            quals[i] = 20;
        read.setBaseQualities(quals);

        RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();
        List<Covariate> requiredCovariates = new ArrayList<Covariate>();
        List<Covariate> optionalCovariates = new ArrayList<Covariate>();

        dataManager = new RecalDataManager(true, 4);

        rgCovariate = new ReadGroupCovariate();
        rgCovariate.initialize(RAC);
        requiredCovariates.add(rgCovariate);

        qsCovariate = new QualityScoreCovariate();
        qsCovariate.initialize(RAC);
        requiredCovariates.add(qsCovariate);

        cxCovariate = new ContextCovariate();
        cxCovariate.initialize(RAC);
        optionalCovariates.add(cxCovariate);
        cyCovariate = new CycleCovariate();
        cyCovariate.initialize(RAC);
        optionalCovariates.add(cyCovariate);

        final Covariate[] requestedCovariates = new Covariate[requiredCovariates.size() + optionalCovariates.size()];
        int covariateIndex = 0;
        for (final Covariate cov : requiredCovariates)
            requestedCovariates[covariateIndex++] = cov;
        for (final Covariate cov : optionalCovariates)
            requestedCovariates[covariateIndex++] = cov;

        readCovariates = RecalDataManager.computeCovariates(read, requestedCovariates);

        RecalibrationTables recalibrationTables = new RecalibrationTables(requestedCovariates);
        final NestedIntegerArray<RecalDatum> rgTable = recalibrationTables.getTable(RecalibrationTables.TableType.READ_GROUP_TABLE);
        final NestedIntegerArray<RecalDatum> qualTable = recalibrationTables.getTable(RecalibrationTables.TableType.QUALITY_SCORE_TABLE);

        for (int i=0; i<read.getReadLength(); i++) {
            final int[] bitKeys = readCovariates.getMismatchesKeySet(i);
            final Object[] objKey = buildObjectKey(bitKeys);

            Random random = new Random();
            final int nObservations = random.nextInt(10000);
            final int nErrors = random.nextInt(10);
            final byte estimatedQReported = 30;
            double empiricalQuality = calcEmpiricalQual(nObservations, nErrors);

            RecalDatum oldDatum = new RecalDatum(nObservations, nErrors, estimatedQReported, empiricalQuality);
            dataManager.addToAllTables(objKey, oldDatum, QualityUtils.MIN_USABLE_Q_SCORE);

            RecalDatum newDatum = new RecalDatum(nObservations, nErrors, estimatedQReported);

            rgTable.put(newDatum, bitKeys[0], EventType.BASE_SUBSTITUTION.index);
            qualTable.put(newDatum, bitKeys[0], bitKeys[1], EventType.BASE_SUBSTITUTION.index);
            for (int j = 0; j < optionalCovariates.size(); j++) {
                final NestedIntegerArray<RecalDatum> covTable = recalibrationTables.getTable(RecalibrationTables.TableType.OPTIONAL_COVARIATE_TABLES_START.index + j);
                covTable.put(newDatum, bitKeys[0], bitKeys[1], j, bitKeys[RecalibrationTables.TableType.OPTIONAL_COVARIATE_TABLES_START.index + j], EventType.BASE_SUBSTITUTION.index);
            }
        }

        dataManager.generateEmpiricalQualities(1, QualityUtils.MAX_RECALIBRATED_Q_SCORE);

        List<Byte> quantizedQuals = new ArrayList<Byte>();
        List<Long> qualCounts = new ArrayList<Long>();
        for (byte i = 0; i <= QualityUtils.MAX_QUAL_SCORE; i++) {
            quantizedQuals.add(i);
            qualCounts.add(1L);
        }
        QuantizationInfo quantizationInfo = new QuantizationInfo(quantizedQuals, qualCounts);
        quantizationInfo.noQuantization();
        baseRecalibration = new BaseRecalibration(quantizationInfo, recalibrationTables, requestedCovariates);

    }


    @Test(enabled=false)
    public void testGoldStandardComparison() {
        for (int i = 0; i < read.getReadLength(); i++) {
            int [] bitKey = readCovariates.getKeySet(i, EventType.BASE_SUBSTITUTION);
            Object [] objKey = buildObjectKey(bitKey);
            byte v2 = baseRecalibration.performSequentialQualityCalculation(bitKey, EventType.BASE_SUBSTITUTION);
            byte v1 = goldStandardSequentialCalculation(objKey);
            Assert.assertEquals(v2, v1);
        }
    }

    private Object[] buildObjectKey(final int[] bitKey) {
        Object[] key = new Object[bitKey.length];
        key[0] = rgCovariate.formatKey(bitKey[0]);
        key[1] = qsCovariate.formatKey(bitKey[1]);
        key[2] = cxCovariate.formatKey(bitKey[2]);
        key[3] = cyCovariate.formatKey(bitKey[3]);
        return key;
    }

    /**
     * Implements a serial recalibration of the reads using the combinational table.
     * First, we perform a positional recalibration, and then a subsequent dinuc correction.
     *
     * Given the full recalibration table, we perform the following preprocessing steps:
     *
     * - calculate the global quality score shift across all data [DeltaQ]
     * - calculate for each of cycle and dinuc the shift of the quality scores relative to the global shift
     * -- i.e., DeltaQ(dinuc) = Sum(pos) Sum(Qual) Qempirical(pos, qual, dinuc) - Qreported(pos, qual, dinuc) / Npos * Nqual
     * - The final shift equation is:
     *
     * Qrecal = Qreported + DeltaQ + DeltaQ(pos) + DeltaQ(dinuc) + DeltaQ( ... any other covariate ... )
     *
     * @param key The list of Comparables that were calculated from the covariates
     * @return A recalibrated quality score as a byte
     */
    private byte goldStandardSequentialCalculation(final Object... key) {

        final byte qualFromRead = (byte) Integer.parseInt(key[1].toString());
        final Object[] readGroupCollapsedKey = new Object[1];
        final Object[] qualityScoreCollapsedKey = new Object[2];
        final Object[] covariateCollapsedKey = new Object[3];

        // The global quality shift (over the read group only)
        readGroupCollapsedKey[0] = key[0];
        final org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatum globalRecalDatum = ((org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatum) dataManager.getCollapsedTable(0).get(readGroupCollapsedKey));
        double globalDeltaQ = 0.0;
        if (globalRecalDatum != null) {
            final double globalDeltaQEmpirical = globalRecalDatum.getEmpiricalQuality();
            final double aggregrateQReported = globalRecalDatum.getEstimatedQReported();
            globalDeltaQ = globalDeltaQEmpirical - aggregrateQReported;
        }

        // The shift in quality between reported and empirical
        qualityScoreCollapsedKey[0] = key[0];
        qualityScoreCollapsedKey[1] = key[1];
        final org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatum qReportedRecalDatum = ((org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatum) dataManager.getCollapsedTable(1).get(qualityScoreCollapsedKey));
        double deltaQReported = 0.0;
        if (qReportedRecalDatum != null) {
            final double deltaQReportedEmpirical = qReportedRecalDatum.getEmpiricalQuality();
            deltaQReported = deltaQReportedEmpirical - qualFromRead - globalDeltaQ;
        }

        // The shift in quality due to each covariate by itself in turn
        double deltaQCovariates = 0.0;
        double deltaQCovariateEmpirical;
        covariateCollapsedKey[0] = key[0];
        covariateCollapsedKey[1] = key[1];
        for (int iii = 2; iii < key.length; iii++) {
            covariateCollapsedKey[2] = key[iii]; // The given covariate
            final org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatum covariateRecalDatum = ((org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatum) dataManager.getCollapsedTable(iii).get(covariateCollapsedKey));
            if (covariateRecalDatum != null) {
                deltaQCovariateEmpirical = covariateRecalDatum.getEmpiricalQuality();
                deltaQCovariates += (deltaQCovariateEmpirical - qualFromRead - (globalDeltaQ + deltaQReported));
            }
        }

        final double newQuality = qualFromRead + globalDeltaQ + deltaQReported + deltaQCovariates;
        return QualityUtils.boundQual((int) Math.round(newQuality), QualityUtils.MAX_RECALIBRATED_Q_SCORE);

        // Verbose printouts used to validate with old recalibrator
        //if(key.contains(null)) {
        //    System.out.println( key  + String.format(" => %d + %.2f + %.2f + %.2f + %.2f = %d",
        //                 qualFromRead, globalDeltaQ, deltaQReported, deltaQPos, deltaQDinuc, newQualityByte));
        //}
        //else {
        //    System.out.println( String.format("%s %s %s %s => %d + %.2f + %.2f + %.2f + %.2f = %d",
        //                 key.get(0).toString(), key.get(3).toString(), key.get(2).toString(), key.get(1).toString(), qualFromRead, globalDeltaQ, deltaQReported, deltaQPos, deltaQDinuc, newQualityByte) );
        //}

        //return newQualityByte;
    }

    public static double calcEmpiricalQual(final int observations, final int errors) {
        final int smoothing = 1;
        final double doubleMismatches = (double) (errors + smoothing);
        final double doubleObservations = (double) ( observations + smoothing );
        double empiricalQual = -10 * Math.log10(doubleMismatches / doubleObservations);
        return Math.min(QualityUtils.MAX_RECALIBRATED_Q_SCORE, empiricalQual);
    }
}
