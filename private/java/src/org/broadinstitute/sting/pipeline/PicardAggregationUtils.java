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

package org.broadinstitute.sting.pipeline;

import edu.mit.broad.picard.util.PicardAggregationFsUtil;
import net.sf.picard.io.IoUtil;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;

public class PicardAggregationUtils {
    private static final PicardAggregationFsUtil aggregationFsUtil = new PicardAggregationFsUtil();
    public static final String PICARD_AGGREGATION_DIR = aggregationFsUtil.AGGREGATION_DIRECTORY.getAbsolutePath() + "/";

    /**
     * Returns the path to the sample BAM.
     * @param project Project
     * @param sample Sample
     * @param version Version
     * @return The path to the sample BAM.
     */
    public static String getSampleBam(String project, String sample, int version) {
        return getSampleDir(project, sample, version) + IoUtil.makeFileNameSafe(sample) + ".bam";
    }

    /**
     * Returns the path to the latest BAM.
     * @param project Project
     * @param sample Sample
     * @return The path to the latest BAM.
     * @throws FileNotFoundException If a finished directory cannot be found for a sample.
     */
    public static String getSampleBam(String project, String sample) throws FileNotFoundException {
        return getSampleDir(project, sample) + IoUtil.makeFileNameSafe(sample) + ".bam";
    }

    /**
     * Returns the sample directory.
     * @param project Project
     * @param sample Sample
     * @param version Version
     * @return the sample directory.
     */
    public static String getSampleDir(String project, String sample, int version) {
        return PICARD_AGGREGATION_DIR + String.format("%s/%s/v%d/", IoUtil.makeFileNameSafe(project), IoUtil.makeFileNameSafe(sample), version);
    }

    /**
     * Returns the latest finished directory for this project sample.
     * @param project Project
     * @param sample Sample
     * @return The path to the latest finished directory.
     * @throws FileNotFoundException If a finished directory cannot be found for a sample.
     */
    public static String getSampleDir(String project, String sample) throws FileNotFoundException {
        int latestVersion = getLatestVersion(project, sample);
        if (latestVersion == 0)
            throw new FileNotFoundException("Unable to find a finished directory for project sample " + project + "/" + sample);
        return getSampleDir(project, sample, latestVersion);
    }

    /**
     * Returns the latest finished version directory
     * @param project Project
     * @param sample Sample
     * @return The highest finished version directory after startVersion
     */
    public static int getLatestVersion(String project, String sample) {
        File sampleDirectory = new File(PICARD_AGGREGATION_DIR + project + "/" + IoUtil.makeFileNameSafe(sample));

        if (!sampleDirectory.exists())
            return 0;

        final File[] versions = sampleDirectory.listFiles(new FileFilter() {
            @Override public boolean accept(final File file) {
                return file.getName().startsWith("v") && file.isDirectory();
            }
        });

        int latestVersion = 0;
        File latest = null;

        for (final File f : versions) {
            if (!aggregationFsUtil.isFinished(f)) continue;

            final int v = Integer.parseInt(f.getName().substring(1));
            if (latest == null || v > latestVersion) {
                latestVersion = v;
                latest = f;
            }
        }

        return latestVersion;
    }

    /**
     * Returns true if the project sample directory contains a finished.txt
     * @param project Project
     * @param sample Sample
     * @param version Version
     * @return true if the project sample directory contains a finished.txt
     */
    public static boolean isFinished(String project, String sample, int version) {
        return aggregationFsUtil.isFinished(new File(getSampleDir(project, sample, version)));
    }
}
