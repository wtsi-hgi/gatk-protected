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

package org.broadinstitute.sting.gatk.walkers.annotator;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.apache.commons.math.linear.*;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ExperimentalAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.GenotypeAnnotation;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.variant.vcf.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.util.*;

/**
 * This annotation attempts to generate a genotype dosage that is controlled for LD for analyses
 * such as variance analysis (e.g. animal model using relationship matrices) or association tests.
 *
 * Ideally one would want to include only the component of the genotypes not explained by the rest of the genome, e.g.
 * use in analysis the residuals of
 *
 * SNP_target ~ SNP1 + SNP2 + ... + SNP_K
 *
 * Under the assumption that LD is local (within say 150KB), we can correct for the LD structure
 * using only the SNPs in a window of 150KB.
 *
 * Furthermore, this can be done in a sliding fashion, with the already-calculated residuals used to model the remaining
 * that is:
 *
 * SNP_{k+1} = SNP_1 + Resid_2 + Resid_3 + ... + Resid_k
 *
 * It follows by symmetry that the sum of squared residuals should be the same.
 */
public class LDCorrectedDosage extends GenotypeAnnotation implements ExperimentalAnnotation {

    private static final int MAX_DISTANCE_IN_BP = 7500;
    private static final String VCF_KEY_NAME = "LDCD";
    private static final String VCF_HEADER_DESCRIPTION = "LD-Controlled Genotype Dosages";
    private TreeSet<GenomeLoc> matrixPositions = new TreeSet<GenomeLoc>();
    private HashMap<String,Integer> sampleOffsets = null;
    private RealMatrix predictorMatrix = null ;
    private Logger logger = Logger.getLogger(this.getClass());

    public List<String> getKeyNames() { return Arrays.asList(VCF_KEY_NAME); }

    public List<VCFFormatHeaderLine> getDescriptions() {
        VCFFormatHeaderLine formatLine = new VCFFormatHeaderLine(VCF_KEY_NAME, VCFHeaderLineCount.A, VCFHeaderLineType.Float,VCF_HEADER_DESCRIPTION);
        return Arrays.asList(formatLine);
    }

     public void annotate(final RefMetaDataTracker tracker,
                         final AnnotatorCompatible walker,
                         final ReferenceContext ref,
                         final AlignmentContext stratifiedContext,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final PerReadAlleleLikelihoodMap alleleLikelihoodMap) {

         // initialize the matrix if we have to
         if ( predictorMatrix == null ) {
             predictorMatrix = initializeMatrix(vc);
             sampleOffsets = initializeSampleMap(vc);
         }

         // check if we need to calculate the residuals
         if ( ! matrixPositions.contains(ref.getLocus()) ) {

             // check the VC for validity
             String invalidVCReason = checkVCIsValid(vc);
             if (invalidVCReason != null ) {
                 throw new UserException.MalformedVCF(String.format("The input VCF is currently not usable with LDCorrectedDosage: %s",invalidVCReason));
             }

             // at this point, VC must be biallelic and have an alternate allele, however there must be some variance
             if ( vc.getCalledChrCount(vc.getAlternateAllele(0)) < 1 || vc.getCalledChrCount(vc.getAlternateAllele(0)) == vc.getCalledChrCount() )
                 return;

             // Extract the genotypes into a matrix
             RealVector newGenotypes = extractGenotypeVector(vc);

             /*double[][] matrix = predictorMatrix.getData();
             for ( int i = 0; i < matrix.length; i++ ) {
                 StringBuffer buf = new StringBuffer();
                 for ( int j = 0; j < matrix[0].length; j++ ) {
                    buf.append(matrix[i][j]);
                    buf.append(" ");
                 }
                 System.out.println(buf);
             }*/
             // Calculate the residuals
             RealMatrix residuals = calculateLeastSquaresResiduals(newGenotypes);

             // remove any genotypes beyond the window
             Set<GenomeLoc> toRemove = new HashSet<GenomeLoc>();
             for ( GenomeLoc distalLoc : matrixPositions ) {
                 if ( distalLoc.distance(ref.getLocus()) > MAX_DISTANCE_IN_BP )
                     toRemove.add(distalLoc);
                 else
                    break;
             }
             predictorMatrix = removeRows(predictorMatrix,toRemove.size());
             matrixPositions.removeAll(toRemove);

             // insert residuals into the predictor matrix
             predictorMatrix = extendMatrix(predictorMatrix,residuals);

             // add the locus to the list
             matrixPositions.add(ref.getLocus());
         }

         // grab the residual from the matrix. The row is the last row, and the column is the sample offset.
         double resid = predictorMatrix.getEntry(predictorMatrix.getRowDimension()-1,sampleOffsets.get(g.getSampleName())-1);

         // stick the residual into the builder
         gb.attribute(VCF_KEY_NAME,resid);

     }

    /**
     * Wrapper method to safely call into OLS library and deal with potential errors
     * coming out of the OLS solver which may necessitate modifications to the
     * predictor matrix. Passes @genotypeDosages and @predictorMatrix through to
     * runOLS().
     * @param genotypeDosages - a (N_SAMPLE X 1) VECTOR of genotype dosages
     * @modified predictorMatrix - a (N_VARIANT X N_SAMPLE) MATRIX of genotype dosages. May be modified in response
     *           to errors caught from the OLS method.
     * @return - residuals from an OLS of genotype dosages on the predictor matrix
     */
    @Requires({"genotypeDosages != null","genotypeDosages.getDimension() == sampleOffsets.size()"})
    @Ensures({"result != null"})
    private RealMatrix calculateLeastSquaresResiduals(final RealVector genotypeDosages) {
        try {
            return runOLS(predictorMatrix,genotypeDosages,true);
        } catch ( InvalidMatrixException e ) {
            logger.warn("Error in solving OLS: Decomposition takes too many iterations. Removing most recently added variant residuals.");
            predictorMatrix = removeRows(predictorMatrix,1);
            if ( predictorMatrix.getRowDimension() <= 1 ) {
                throw new ReviewedStingException("All rows of predictor matrix removed, yet SVD not converging.");
            }
            return calculateLeastSquaresResiduals(genotypeDosages);
            //logger.warn(Arrays.toString(predictorMatrix.getRow(predictorMatrix.getRowDimension() - 1)));
            //throw new ReviewedStingException(e.getMessage(),e);
        }
    }

    /**
     * Runs ordinary least squares regression of @response on @predict, and returns the residuals
     * (response - predict*beta)
     * @param predict - a (N_VARIANT X N_SAMPLE) MATRIX of (possibly corrected) genotype dosages
     * @param response - a (N_SAMPLE x 1) VECTOR of genotype dosages
     * @param rectify - rectify the predicted values (cap at 0, 2) prior to calculating residuals
     * @return a (1 X N_SAMPLE) MATRIX holding the OLS residuals of @response on @predict, hinge rectified if @rectify
     */
    @Requires({"predict != null", "response != null","predict.getColumnDimension()==response.getDimension()"})
    @Ensures({"result != null","result.getRowDimension() == 1","result.getColumnDimension() == response.getDimension()"})
    protected static RealMatrix runOLS(final RealMatrix predict, final RealVector response, boolean rectify) {
        final RealMatrix tPredict = predict.transpose();
        DecompositionSolver solver = new SingularValueDecompositionImpl(tPredict).getSolver();
        RealVector coefficients = solver.solve(response);
        RealVector predictedResponse = tPredict.operate(coefficients);
        if ( rectify )
            predictedResponse = hingeRectify(predictedResponse);
        RealVector residuals = response.subtract(predictedResponse);
        return new Array2DRowRealMatrix(residuals.getData()).transpose();
    }

    /**
     * Rectifies the input vector via the hinge function
     *             |- 0   x < 0
     *     f(x) =  |- x   0 < x < 2
     *             |- 2   x > 2
     * @param vector - an input vector to be rectified
     * @return the vector formed by taking f(@vector) elementwise
     */
    @Requires({"vector != null"})
    @Ensures({"result != null","result.getDimension() == vector.getDimension()"})
    protected static RealVector hingeRectify(RealVector vector) {
        double[] rectifiedData = new double[vector.getDimension()];
        int idx = 0;
        for ( double datum : vector.getData() ) {
            if ( datum < 0.0 ) {
                rectifiedData[idx++] = 0.;
            } else if ( datum > 2.0 ) {
                rectifiedData[idx++] = 2.0;
            } else {
                rectifiedData[idx++] = datum;
            }
        }

        return new ArrayRealVector(rectifiedData);
    }

    /**
     * Removes a @numRowsToRemove from the matrix @baseMatrix, starting with the second row (the first row
     * is always reserved as the row of ones), and returns the result, which is not necessarily a clone.
     * @param baseMatrix - contractually a (N_VARIANT X N_SAMPLE) MATRIX of genotype dosages
     * @param numRowsToRemove - The number of rows to remove
     * @return - the submatrix formed by repeatedly removing the second row of @baseMatrix @numRowsToRemove number of times
     */
    @Requires({"baseMatrix != null","numRowsToRemove >= 0","numRowsToRemove < baseMatrix.getRowDimension()",
               "baseMatrix.getColumnDimension() == sampleOffsets.size()"})
    @Ensures({"result != null","result.getRowDimension() == baseMatrix.getRowDimension()-numRowsToRemove",
              "result.getColumnDimension() == baseMatrix.getColumnDimension()"})
    private RealMatrix removeRows(RealMatrix baseMatrix, final int numRowsToRemove) {
        int[] rowsToKeep = new int[baseMatrix.getRowDimension()-numRowsToRemove];
        rowsToKeep[0]=0;
        int idx = 1;
        for ( int i = 1+numRowsToRemove; i < baseMatrix.getRowDimension(); i++ ) {
            rowsToKeep[idx++]=i;
        }
        int[] colsToKeep = new int[baseMatrix.getColumnDimension()];
        for ( int i = 0; i < colsToKeep.length; i++ ) {
            colsToKeep[i] = i;
        }
        return baseMatrix.getSubMatrix(rowsToKeep,colsToKeep);
    }

    /**
     * Appends the matrix @extension row-wise to the bottom of @baseMatrix, and returns the result
     * @param baseMatrix - contractually a (N_VARIANT X N_SAMPLE) MATRIX of genotype dosages
     * @param extension - contractually a (N_VARIANT X 1) MATRIX of genotype dosages
     * @return - the matrix formed by concatenating @extension to the bottom of @baseMatrix
     */
    @Requires({"baseMatrix != null","extension != null","baseMatrix.getColumnDimension() == sampleOffsets.size()",
                "baseMatrix.getRowDimension() > 0","extension.getColumnDimension() == sampleOffsets.size()",
                "extension.getRowDimension() == 1"})
    @Ensures({"result != null","result.getRowDimension() == baseMatrix.getRowDimension() + extension.getRowDimension()"})
    private RealMatrix extendMatrix(final RealMatrix baseMatrix, final RealMatrix extension) {
        int nRow = baseMatrix.getRowDimension() + extension.getRowDimension();
        RealMatrix extendedMatrix = baseMatrix.createMatrix(nRow,baseMatrix.getColumnDimension());
        extendedMatrix.setSubMatrix(baseMatrix.getData(),0,0);
        extendedMatrix.setSubMatrix(extension.getData(),baseMatrix.getRowDimension(),0);
        return extendedMatrix;
    }

    /**
     * Extracts the genotype dosages from the given variant context, @vc, and returns the result as a
     * RealVector object. Missing genotypes are assigned the average dosage. Requires bi-allelic sites.
     * @param vc - a VariantContext containing genotypes
     * @return - a (N_SAMPLES X 1) VECTOR of genotype dosages
     */
    @Requires({"vc != null","vc.hasGenotypes()","vc.getNSamples() > vc.getNoCallCount()","vc.getNSamples() == sampleOffsets.size()"})
    @Ensures({"result != null","result.getDimension() == vc.getNSamples()"})
    private RealVector extractGenotypeVector(final VariantContext vc) {
        double averageDosage = 0.0;
        final double denom = 1.0/(vc.getNSamples()-vc.getNoCallCount());
        for ( Genotype genotype : vc.getGenotypes() ) {
            averageDosage += genotype.isNoCall() ? 0 : genotype.countAllele(vc.getAlternateAllele(0))/denom;
        }

        double[] genotypeVector = new double[vc.getNSamples()];
        int idx = 0;
        for ( Genotype genotype : vc.getGenotypes() ) {
            genotypeVector[idx++] = genotype.isNoCall() ? averageDosage : (double) genotype.countAllele(vc.getAlternateAllele(0));
        }

        return new ArrayRealVector(genotypeVector);
    }

    /**
     * Initialzes a (1 X N_SAMPLE) MATRIX from a given variant context @context
     * @param context - a variant context containing genotypes
     * @return - the (1 X N_SAMPLE) MATRIX of ones
     */
    @Requires({"context != null","context.hasGenotypes()"})
    @Ensures({"result != null","result.getRowDimension() == 1","result.getColumnDimension() == context.getNSamples()",
              "result.getEntry(0,0)==1.0"})
    private RealMatrix initializeMatrix(final VariantContext context) {
        // need to know the number of samples in order to initialize the matrix properly
        double[][] ones = new double[1][context.getNSamples()];
        for ( int i = 0; i < ones[0].length; i++ ) {
            ones[0][i] = 1.0;
        }
        // the zeroth row of the matrix will always be a row of ones (the intercept)
        return new Array2DRowRealMatrix(ones);
    }

    @Requires({"context != null","context.hasGenotypes()"})
    @Ensures({"result.size() == context.getNSamples()"})
    private HashMap<String,Integer> initializeSampleMap(final VariantContext context) {
        HashMap<String,Integer> sampleMap = new HashMap<String,Integer>(context.getNSamples());
        int offset = 1; // remember, offset 0 is the ones vector
        for ( Genotype genotype : context.getGenotypes()) {
            sampleMap.put(genotype.getSampleName(),offset);
            ++offset;
        }

        return sampleMap;
    }

    private String checkVCIsValid(final VariantContext context) {
        if ( ! context.hasGenotypes() )
            return String.format("Context at position %s:%d does not have genotypes, which is unsupported by this annotation.",context.getChr(),context.getStart());
        if ( context.getNoCallCount() == context.getNSamples() )
            return String.format("Context at position %s:%d consists of entirely no calls, which is unsupported by this annotation.",context.getChr(),context.getStart());
        if ( ! context.isBiallelic() )
            return String.format("Context at position %s:%d is not bi-allelic, which is unsupported by this annotation.",context.getChr(),context.getStart());
        return null;
    }

}
