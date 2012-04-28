/*
 * Copyright (c) 2011 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.variantrecalibrator3.variantrecalibration;

import ca.mcgill.mcb.pcingola.collections.ArrayUtil;
import org.apache.commons.lang.ArrayUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.collections.ExpandingArrayList;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Mar 4, 2011
 */

public class VariantDataManager {
    private ExpandingArrayList<VariantDatum> data;
    private final double[] meanVector;
    private final double[] varianceVector;
    public final List<String> initialKeys;
    private final VariantRecalibratorArgumentCollection VRAC;
    protected final static Logger logger = Logger.getLogger(VariantDataManager.class);
    protected final List<String> finalKeys;
    protected final List<String> allKeys;
    protected final List<TrainingSet> trainingSets;
    protected final boolean[] finalTraining;
    protected final boolean[] initialTraining;
    protected int[] obsCounts;

    public VariantDataManager( final List<String> initialKeys, final List<String> finalKeys, final VariantRecalibratorArgumentCollection VRAC ) {
        this.data = null;
        this.initialKeys = new ArrayList<String>( initialKeys );
        this.finalKeys = new ArrayList<String>(finalKeys);
        List<String> allKeys = new ArrayList<String>(initialKeys);
        for ( String key : finalKeys ) {
            if ( ! allKeys.contains(key) ) {
                allKeys.add(key);
            }
        }
        this.allKeys = allKeys;
        boolean[] ft = new boolean[allKeys.size()];
        boolean[] it = new boolean[allKeys.size()];
        int idx = 0;
        for ( String k : allKeys ) {
            it[idx] = initialKeys.contains(k);
            ft[idx] = finalKeys.contains(k);
            ++idx;
        }
        this.finalTraining = ft;
        this.initialTraining = it;
        this.VRAC = VRAC;
        obsCounts = new int[this.allKeys.size()];
        meanVector = new double[this.allKeys.size()];
        varianceVector = new double[this.allKeys.size()];
        trainingSets = new ArrayList<TrainingSet>();
    }

    public void setData( final ExpandingArrayList<VariantDatum> data ) {
        this.data = data;
    }

    public ExpandingArrayList<VariantDatum> getData() {
        return data;
    }

    public void normalizeData() {
        boolean foundZeroVarianceAnnotation = false;
        for( int iii = 0; iii < meanVector.length; iii++ ) {
            final double theMean = meanVector[iii];
            final double theSTD = Math.sqrt(varianceVector[iii]/(obsCounts[iii]-1));
            logger.info( allKeys.get(iii) + String.format(": \t mean = %.8f\t standard deviation = %.8f", theMean, theSTD) );
            if( obsCounts[iii] == 0 ) {
                throw new UserException.BadInput("Values for " + allKeys.get(iii) + " annotation not detected for ANY training variant in the input callset. VariantAnnotator may be used to add these annotations. See http://www.broadinstitute.org/gsa/wiki/index.php/VariantAnnotator");
            }

            foundZeroVarianceAnnotation |= (theSTD < 1E-6);
            int initIdx = -1;
            int finIdx = -1;
            if ( initialTraining[iii] )
                initIdx = initialKeys.lastIndexOf(allKeys.get(iii));
            if ( finalTraining[iii] )
                finIdx = finalKeys.lastIndexOf(allKeys.get(iii));
            for( final VariantDatum datum : data ) {
                // Transform each data point via: (x - mean) / standard deviation
                normalizeDatum(datum,theMean,theSTD,initIdx,finIdx);
            }
        }
        if( foundZeroVarianceAnnotation ) {
            throw new UserException.BadInput( "Found annotations with zero variance. They must be excluded before proceeding." );
        }

        // trim data by standard deviation threshold and mark failing data for exclusion later
        for( final VariantDatum datum : data ) {
            boolean remove = false;
            for( final double val : datum.initialAnnotations ) {
                remove = remove || (Math.abs(val) > VRAC.STD_THRESHOLD);
            }
            /**
             * todo -- i'm pretty sure we don't want to threshold on the final annotations
             *
            for ( final double val : datum.finalAnnotations ) {
                remove = remove || (Math.abs(val) > VRAC.STD_THRESHOLD);
            }
             */
            datum.failingSTDThreshold = remove;
        }
    }

    private void normalizeDatum(VariantDatum datum, double mean, double std, int initIdx, int finalIdx) {
        if ( initIdx > -1 ) {
            datum.initialAnnotations[initIdx] = datum.isNullInitial[initIdx] ? GenomeAnalysisEngine.getRandomGenerator().nextGaussian() : (datum.initialAnnotations[initIdx] - mean)/std;
        }

        /*if ( finalIdx > -1 ) {
            datum.finalAnnotations[finalIdx] = datum.isNullFinal[finalIdx] ? GenomeAnalysisEngine.getRandomGenerator().nextGaussian() : (datum.finalAnnotations[finalIdx] - mean)/std;
        }*/ // don't actually want to normalize the final data
    }

     public void addTrainingSet( final TrainingSet trainingSet ) {
         trainingSets.add( trainingSet );
     }

     public boolean checkHasTrainingSet() {
         for( final TrainingSet trainingSet : trainingSets ) {
             if( trainingSet.isTraining ) { return true; }
         }
         return false;
     }

     public boolean checkHasTruthSet() {
         for( final TrainingSet trainingSet : trainingSets ) {
             if( trainingSet.isTruth ) { return true; }
         }
         return false;
     }

     public boolean checkHasKnownSet() {
         for( final TrainingSet trainingSet : trainingSets ) {
             if( trainingSet.isKnown ) { return true; }
         }
         return false;
     }

    public ExpandingArrayList<VariantDatum> getTrainingData() {
        final ExpandingArrayList<VariantDatum> trainingData = new ExpandingArrayList<VariantDatum>();
        for( final VariantDatum datum : data ) {
            if( (datum.atPositiveTrainingSite || datum.atNegativeTrainingSite) && !datum.failingSTDThreshold && datum.originalQual > VRAC.QUAL_THRESHOLD ) {
                trainingData.add( datum );
            }
        }
        logger.info( "Training with " + trainingData.size() + " variants after standard deviation thresholding." );
        if( trainingData.size() < VRAC.MIN_NUM_BAD_VARIANTS ) {
            logger.warn( "WARNING: Training with very few variant sites! Please check the model reporting PDF to ensure the quality of the model is reliable." );
        }
        return trainingData;
    }

    public ExpandingArrayList<VariantDatum> selectWorstVariants( double bottomPercentage, final int minimumNumber ) {
        // The return value is the list of training variants
        final ExpandingArrayList<VariantDatum> trainingData = new ExpandingArrayList<VariantDatum>();

        // First add to the training list all sites overlapping any bad sites training tracks
        for( final VariantDatum datum : data ) {
            if( datum.atNegativeTrainingSite && !datum.failingSTDThreshold && !Double.isInfinite(datum.getLod()) ) {
                trainingData.add( datum );
            }
        }
        final int numBadSitesAdded = trainingData.size();
        logger.info( "Found " + numBadSitesAdded + " variants overlapping bad sites training tracks." );

        // Next sort the variants by the LOD coming from the positive model and add to the list the bottom X percent of variants
        Collections.sort( data );
        final int numToAdd = Math.max( minimumNumber - trainingData.size(), Math.round((float)bottomPercentage * data.size()) );
        if( numToAdd > data.size() ) {
            throw new UserException.BadInput( "Error during negative model training. Minimum number of variants to use in training is larger than the whole call set. One can attempt to lower the --minNumBadVariants arugment but this is unsafe." );
        } else if( numToAdd == minimumNumber - trainingData.size() ) {
            logger.warn( "WARNING: Training with very few variant sites! Please check the model reporting PDF to ensure the quality of the model is reliable." );
            bottomPercentage = ((float) numToAdd) / ((float) data.size());
        }
        int index = 0, numAdded = 0;
        while( numAdded < numToAdd ) {
            final VariantDatum datum = data.get(index++);
            if( !datum.atNegativeTrainingSite && !datum.failingSTDThreshold && !Double.isInfinite(datum.getLod()) ) {
                datum.atNegativeTrainingSite = true;
                trainingData.add( datum );
                numAdded++;
            }
        }
        logger.info( "Additionally training with worst " + String.format("%.3f", (float) bottomPercentage * 100.0f) + "% of passing data --> " + (trainingData.size() - numBadSitesAdded) + " variants with LOD <= " + String.format("%.4f", data.get(index).getLod()) + "." );
        return trainingData;
    }

    public ExpandingArrayList<VariantDatum> getRandomDataForPlotting( int numToAdd ) {
        numToAdd = Math.min(numToAdd, data.size());
        final ExpandingArrayList<VariantDatum> returnData = new ExpandingArrayList<VariantDatum>();
        for( int iii = 0; iii < numToAdd; iii++) {
            final VariantDatum datum = data.get(GenomeAnalysisEngine.getRandomGenerator().nextInt(data.size()));
            if( !datum.failingSTDThreshold ) {
                returnData.add(datum);
            }
        }

        // Add an extra 5% of points from bad training set, since that set is small but interesting
        for( int iii = 0; iii < Math.floor(0.05*numToAdd); iii++) {
            final VariantDatum datum = data.get(GenomeAnalysisEngine.getRandomGenerator().nextInt(data.size()));
            if( datum.atNegativeTrainingSite && !datum.failingSTDThreshold ) { returnData.add(datum); }
            else { iii--; }
        }

        return returnData;
    }

    public void decodeAnnotations( final VariantDatum datum, final VariantContext vc, final boolean jitter ) {
        final double[] annotations = new double[allKeys.size()];
        final double[] rawAnnotations = new double[allKeys.size()];
        final boolean[] isNull = new boolean[allKeys.size()];
        int iii = 0;
        for( final String key : allKeys ) {
            isNull[iii] = false;
            annotations[iii] = decodeAnnotation( key, vc, jitter,true ); // jitters and does special hard-coded stuff
            rawAnnotations[iii] = decodeAnnotation(key,vc,false,false); // no jitter, no hard-coded stuff
            if( Double.isNaN(annotations[iii]) ) { isNull[iii] = true; }
            iii++;
        }
        if ( datum.atPositiveTrainingSite ) {
            updateStatistics(annotations,isNull);
        }
        segregateInitialFinalAnnotations(annotations,rawAnnotations,isNull,datum);
    }

    private void updateStatistics(double[] vector, boolean[] isNull) {
        for ( int idx = 0; idx < vector.length; idx++ ) {
            if ( isNull[idx] ) {
                continue;
            }
            obsCounts[idx]++;
            double oldMean = meanVector[idx];
            meanVector[idx] += (vector[idx] - meanVector[idx])/obsCounts[idx];
            varianceVector[idx] += (vector[idx] - oldMean)*(vector[idx]-meanVector[idx]);
        }
    }

    public void removeWorstPositiveTrainingSites(double ignorePercentage) {
        ExpandingArrayList<VariantDatum> trainingData = new ExpandingArrayList<VariantDatum>();
        for ( VariantDatum datum : data ) {
            if ( datum.atPositiveTrainingSite ) {
                trainingData.add(datum);
            }
        }

        Collections.sort(trainingData);
        int stop = (int) Math.round(ignorePercentage*trainingData.size());
        logger.info(String.format("Un-classifying the worst %.2f percent of positive training data with LOD<=%.2f",100*ignorePercentage,trainingData.get(stop).getLod()));
        for ( int idx = 0; idx < stop; idx++ ) {
            trainingData.get(idx).atPositiveTrainingSite=false;
        }
    }

    public void reclassifyBestPositiveSites(double topPercentage) {
        int stopIdx = data.size()-1-((int) (topPercentage*data.size()));
        Collections.sort(data);
        logger.info(String.format("Classifying the best %.2f percent of classified sites with LOD >= %.2f",100*topPercentage,data.get(stopIdx).getLod()));
        for ( int ctr = 0; ctr < (int) (topPercentage*data.size()); ctr++ ) {
            int idx = data.size()-1-ctr;
            data.get(idx).atPositiveTrainingSite = true;
        }
    }

    public ExpandingArrayList<VariantDatum> selectWorstVariantsNoThresholding( double bottomPercentage, final int minimumNumber ) {
        // The return value is the list of training variants
        final ExpandingArrayList<VariantDatum> trainingData = new ExpandingArrayList<VariantDatum>();

        // First add to the training list all sites overlapping any bad sites training tracks
        for( final VariantDatum datum : data ) {
            if( datum.atNegativeTrainingSite && !datum.failingSTDThreshold && !Double.isInfinite(datum.getLod()) ) {
                trainingData.add( datum );
            }
        }
        final int numBadSitesAdded = trainingData.size();
        logger.info( "Found " + numBadSitesAdded + " variants overlapping bad sites training tracks." );

        // Next sort the variants by the LOD coming from the positive model and add to the list the bottom X percent of variants
        Collections.sort( data );
        final int numToAdd = Math.max( minimumNumber - trainingData.size(), Math.round((float)bottomPercentage * data.size()) );
        if( numToAdd > data.size() ) {
            throw new UserException.BadInput( "Error during negative model training. Minimum number of variants to use in training is larger than the whole call set. One can attempt to lower the --minNumBadVariants arugment but this is unsafe." );
        } else if( numToAdd == minimumNumber - trainingData.size() ) {
            logger.warn( "WARNING: Training with very few variant sites! Please check the model reporting PDF to ensure the quality of the model is reliable." );
            bottomPercentage = ((float) numToAdd) / ((float) data.size());
        }
        int index = 0, numAdded = 0;
        while( numAdded < numToAdd ) {
            final VariantDatum datum = data.get(index++);
            if( !datum.atNegativeTrainingSite ) {
                datum.atNegativeTrainingSite = true;
                trainingData.add( datum );
                numAdded++;
            }
        }
        logger.info( "Additionally training with worst " + String.format("%.3f", (float) bottomPercentage * 100.0f) + "% of passing data --> " + (trainingData.size() - numBadSitesAdded) + " variants with LOD <= " + String.format("%.4f", data.get(index).getLod()) + "." );
        return trainingData;
    }

    private void segregateInitialFinalAnnotations(double[] annotations, double[] rawAnnotations, boolean[] isNull, VariantDatum datum) {
        double[] init = new double[initialKeys.size()];
        double[] finT = new double[finalKeys.size()];
        boolean[] nullInit = new boolean[initialKeys.size()];
        boolean[] nullFinal = new boolean[finalKeys.size()];
        int initIdx = 0;
        int finIdx = 0;
        for ( int idx = 0; idx < annotations.length; idx++) {
            if ( initialTraining[idx] ) {
                init[initIdx] = annotations[idx];
                nullInit[initIdx] = isNull[idx];
                initIdx++;
            }
            if ( finalTraining[idx] ) {
                finT[finIdx] = rawAnnotations[idx];
                nullFinal[finIdx] = isNull[idx];
                finIdx++;
            }
        }
        datum.initialAnnotations = init;
        datum.isNullInitial = nullInit;
        datum.finalAnnotations = finT;
        datum.isNullFinal = nullFinal;
    }

    private static double decodeAnnotation( final String annotationKey, final VariantContext vc, final boolean jitter, final boolean specialSauce ) {
        double value;

        try {
            value = Double.parseDouble( (String)vc.getAttribute( annotationKey ) );
            if( Double.isInfinite(value) ) { value = Double.NaN; }
            if( jitter && annotationKey.equalsIgnoreCase("HRUN") ) { // Integer valued annotations must be jittered a bit to work in this GMM
                if ( specialSauce ) {
                    value += -0.25 + 0.5 * GenomeAnalysisEngine.getRandomGenerator().nextDouble();
                }
            }

            if (vc.isIndel() && annotationKey.equalsIgnoreCase("QD")) {
            // normalize QD by event length for indel case
                int eventLength = Math.abs(vc.getAlternateAllele(0).getBaseString().length() - vc.getReference().getBaseString().length()); // ignore multi-allelic complication here for now
                if (eventLength > 0 && specialSauce ) { // sanity check
                    value /= (double)eventLength;
                }
            }

            if( specialSauce && jitter && annotationKey.equalsIgnoreCase("HaplotypeScore") && MathUtils.compareDoubles(value, 0.0, 0.0001) == 0 ) { value = -0.2 + 0.4*GenomeAnalysisEngine.getRandomGenerator().nextDouble(); }
            if( specialSauce && jitter && annotationKey.equalsIgnoreCase("FS") && MathUtils.compareDoubles(value, 0.0, 0.001) == 0 ) { value = -0.2 + 0.4*GenomeAnalysisEngine.getRandomGenerator().nextDouble(); }
        } catch( Exception e ) {
            value = Double.NaN; // The VQSR works with missing data by marginalizing over the missing dimension when evaluating the Gaussian mixture model
        }

        return value;
    }

    public void parseTrainingSets( final RefMetaDataTracker tracker, final GenomeLoc genomeLoc, final VariantContext evalVC, final VariantDatum datum, final boolean TRUST_ALL_POLYMORPHIC ) {
        datum.isKnown = false;
        datum.atTruthSite = false;
        datum.atPositiveTrainingSite = false;
        datum.atNegativeTrainingSite = false;
        datum.prior = 2.0;

        for( final TrainingSet trainingSet : trainingSets ) {
            for( final VariantContext trainVC : tracker.getValues(trainingSet.rodBinding, genomeLoc) ) {
                if( isValidVariant( evalVC, trainVC, TRUST_ALL_POLYMORPHIC ) ) {
                    datum.isKnown = datum.isKnown || trainingSet.isKnown;
                    datum.atTruthSite = datum.atTruthSite || trainingSet.isTruth;
                    datum.atPositiveTrainingSite = datum.atPositiveTrainingSite || trainingSet.isTraining;
                    datum.prior = Math.max( datum.prior, trainingSet.prior );
                    datum.consensusCount += ( trainingSet.isConsensus ? 1 : 0 );
                    datum.atNegativeTrainingSite = datum.atNegativeTrainingSite || trainingSet.isAntiTraining;
                    datum.atMonomorphicSite = datum.atMonomorphicSite || trainingSet.isMonomorphic;
                }
            }
        }
    }

    private boolean isValidVariant( final VariantContext evalVC, final VariantContext trainVC, final boolean TRUST_ALL_POLYMORPHIC) {
        return trainVC != null && trainVC.isNotFiltered() && trainVC.isVariant() &&
                        ((evalVC.isSNP() && trainVC.isSNP()) || ((evalVC.isIndel()||evalVC.isMixed()) && (trainVC.isIndel()||trainVC.isMixed()))) &&
                        (TRUST_ALL_POLYMORPHIC || !trainVC.hasGenotypes() || trainVC.isPolymorphicInSamples());
    }

    public void writeOutRecalibrationTable( final PrintStream RECAL_FILE ) {
        for( final VariantDatum datum : data ) {
            RECAL_FILE.println(String.format("%s,%d,%d,%.4f,%s",
                    datum.contig, datum.start, datum.stop, datum.getLod(),
                    (datum.worstAnnotation != -1 ? allKeys.get(datum.worstAnnotation) : "NULL")));
        }
    }
}
