/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.performance;

import org.broad.tribble.Tribble;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.broad.tribble.source.BasicFeatureSource;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.SimpleTimer;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Iterator;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;


/**
 * Tests low-level multi-threading performance of tribble
 *
 * <p>
 * Creates a thread pool that reads the input VCF file in parallel with
 * N threads from 1 to maxThreads and emits the wall time needed to process
 * the entire file.  Assumes the VCF file has a chromosome
 * named 1 that has at least 250 Mb.
 * </p>
 *
 * <h2>Input</h2>
 * <p>
 * An already indexed VCF file for evaluation
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * A simple table of runtimes and detailed information
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 *    java
 *      -jar GenomeAnalysisTK.jar
 *      -T $WalkerName
 *      -vcf my.b37.vcf
 *  </pre>
 *
 * @author Your Name
 * @since Date created
 */
public class RodThreadingPerformanceTest extends RodWalker<Integer, Integer> {
    private final static String CHROMOSOME = "1";
    private final static int MB_TO_BP = 1000000;

    @Output(doc="File to which results should be written",required=true)
    protected PrintStream out;

    /**
     * Chr1 of this file will be used to evaluate performance of multi-threading in the rod system
     */
    @Argument(fullName="vcf", shortName="vcf", doc="VCF input to evalute reading", required=true)
    public File VCFFile;

    /**
     * Evaluate performance of the system using 1, 2, ... maxThreads.
     */
    @Argument(fullName="verbose", shortName="verbose", doc="Should we output per thread statistics", required=false)
    public boolean VERBOSE = false;

    /**
     * Evaluate performance of the system using 1, 2, ... maxThreads.
     */
    @Argument(fullName="maxThreads", shortName="maxThreads", doc="Up to this many threads will be run in parallel", required=false)
    public int maxThreads = 4;

    /**
     * We read --vcf from 1 to maxMB basepairs for evaluation.  For example, -maxMB
     * 10 means that we read from 1-10,000,000 from VCF.
     */
    @Argument(fullName="maxMB", shortName="maxMB", doc="File to test", required=false)
    public int maxMB = 250;

    private Index index;

    public void initialize() {
        index = IndexFactory.loadIndex(Tribble.indexFile(VCFFile).getAbsolutePath());

        for ( int nThreads = 1; nThreads <= maxThreads; nThreads++ ) {
            run(nThreads);
        }
    }

    @Override
    public boolean isDone() {
        return true;
    }

    private void run(int nThreads) {
        ExecutorService pool = Executors.newFixedThreadPool(nThreads);
        int chunkSize = maxMB / nThreads * MB_TO_BP;
        SimpleTimer timer = new SimpleTimer();
        timer.start();
        for ( int i = 0; i < nThreads; i++ ) {
            int start = i * chunkSize + 1;
            int end = start + chunkSize;
            ReadRegion region = new ReadRegion(start, end);
            pool.execute(region);
        }

        pool.shutdown();

        try {
            pool.awaitTermination(1, TimeUnit.MINUTES);
        } catch ( InterruptedException e ) {
            ;
        }
        System.out.printf("TIME: %d thread(s) runtime %.2f seconds%n", nThreads, timer.getElapsedTime());
    }

    public class ReadRegion implements Runnable {
        final int start, end;

        public ReadRegion(final int start, final int end) {
            this.start = start;
            this.end = end;
        }

        public void run() {
            try {
                FeatureCodec codec = new VCFCodec();
                BasicFeatureSource source = new BasicFeatureSource(VCFFile.getAbsolutePath(), index, codec);

                // now read iterate over the file
                long recordCount = 0l;

                // this call could be replaced with a query
                Iterator<Feature> iter = source.query(CHROMOSOME, start, end);

                // cycle through the iterators
                while (iter.hasNext()) {
                    Feature feat = iter.next();
                    ++recordCount;
                }

                if ( VERBOSE )
                    System.out.printf("  THREAD: %d-%d read %d objects%n", start, end, recordCount);
            } catch (IOException e) {
                throw new RuntimeException("Something went wrong while reading feature file " + VCFFile, e);
            }
        }
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) { return 0; }
    public Integer reduceInit() { return 0; }
    public Integer reduce(Integer counter, Integer sum) { return counter + sum; }
}