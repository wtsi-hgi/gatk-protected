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
import org.broadinstitute.sting.gatk.walkers.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypeBuilder;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

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
