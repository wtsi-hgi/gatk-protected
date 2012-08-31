package org.broadinstitute.sting.gatk.walkers.annotator;

import Jama.Matrix;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math.linear.MatrixUtils;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ExperimentalAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.GenotypeAnnotation;
import org.broadinstitute.sting.gatk.walkers.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.codecs.vcf.*;
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

    private static final int MAX_DISTANCE_IN_BP = 75000;
    private static final String VCF_KEY_NAME = "LDCD";
    private static final String VCF_HEADER_DESCRIPTION = "LD-Controlled Genotype Dosages";
    private TreeSet<GenomeLoc> matrixPositions = new TreeSet<GenomeLoc>();
    private HashMap<String,Integer> sampleOffsets = null;
    private Matrix predictorMatrix = null ;

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

         System.out.println("Foo1");

         // initialize the matrix if we have to
         if ( predictorMatrix == null ) {
             predictorMatrix = initializeMatrix(vc);
             sampleOffsets = initializeSampleMap(vc);
         }

         System.out.println("Foo2");

         // check if we need to calculate the residuals
         if ( ! matrixPositions.contains(ref.getLocus()) ) {

             // check the VC for validity
             String invalidVCReason = checkVCIsValid(vc);
             if (invalidVCReason != null ) {
                 throw new UserException.MalformedVCF(String.format("The input VCF is currently not usable with LDCorrectedDosage: %s",invalidVCReason));
             }

             // Extract the genotypes into a matrix
             Matrix newGenotypes = extractGenotypeVector(vc);

             // Calculate the residuals
             Matrix residuals = runOLS(predictorMatrix,newGenotypes);

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

         System.out.println("Foo3");

         // grab the residual from the matrix. The row is the last row, and the column is the sample offset.
         double resid = predictorMatrix.get(predictorMatrix.getRowDimension(),sampleOffsets.get(g.getSampleName()));

         // stick the residual into the builder
         gb.attribute(VCF_KEY_NAME,resid);

     }

    private Matrix runOLS(Matrix predict, Matrix response) {
        Matrix beta = predict.solve(response);
        return response.minus(beta.times(predict));
    }

    // remove rows 1 to numRowsToRemove
    private Matrix removeRows(Matrix baseMatrix, int numRowsToRemove) {
        int[] rowsToKeep = new int[baseMatrix.getRowDimension()-numRowsToRemove];
        rowsToKeep[0]=0;
        int idx = 1;
        for ( int i = 1+numRowsToRemove; i < baseMatrix.getRowDimension(); i++ ) {
            rowsToKeep[idx]=i;
        }
        return baseMatrix.getMatrix(rowsToKeep, 0, baseMatrix.getColumnDimension() - 1);
    }

    private Matrix extendMatrix(Matrix baseMatrix, Matrix extension) {
        int nRow = baseMatrix.getRowDimension() + extension.getRowDimension();
        return new Matrix(ArrayUtils.addAll(baseMatrix.getRowPackedCopy(),extension.getRowPackedCopy()),nRow);
    }

    private Matrix extractGenotypeVector(VariantContext vc) {
        double[] genotypeVector = new double[vc.getNSamples()];
        int idx = 0;
        for ( Genotype genotype : vc.getGenotypes() ) {
            genotypeVector[idx++] = (double) genotype.countAllele(vc.getAlternateAllele(0));
        }

        return new Matrix(genotypeVector,1);
    }

    private Matrix initializeMatrix(VariantContext context) {
        // need to know the number of samples in order to initialize the matrix properly
        double[] ones = new double[context.getNSamples()];
        for ( int i = 0; i < ones.length; i++ ) {
            ones[i] = 1.0;
        }
        // the zeroth row of the matrix will always be a row of ones (the intercept)
        return new Matrix(ones,1);
    }

    private HashMap<String,Integer> initializeSampleMap(VariantContext context) {
        HashMap<String,Integer> sampleMap = new HashMap<String,Integer>(context.getNSamples());
        int offset = 1; // remember, offset 0 is the ones vector
        for ( Genotype genotype : context.getGenotypes()) {
            sampleMap.put(genotype.getSampleName(),offset);
            ++offset;
        }

        return sampleMap;
    }

    private String checkVCIsValid(VariantContext context) {
        if ( context.getNoCallCount() > 0 )
            return "Genotypes may not be no-call.";
        if ( ! context.isBiallelic() )
            return "Sites must be bi-allelic.";
        return null;
    }

}
