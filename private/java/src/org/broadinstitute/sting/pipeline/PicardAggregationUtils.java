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
import org.apache.commons.io.LineIterator;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.io.IOUtils;

import java.io.File;
import java.io.FileFilter;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

public class PicardAggregationUtils {
    private static final Logger log = Logger.getLogger(PicardAggregationUtils.class);

    private static final PicardAggregationFsUtil aggregationFsUtil = new PicardAggregationFsUtil();
    public static final String PICARD_AGGREGATION_DIR = aggregationFsUtil.AGGREGATION_DIRECTORY.getAbsolutePath() + "/";

    public static List<PicardSample> parseSamples(File tsv) {
        List<PicardSample> picardSamples = new ArrayList<PicardSample>();
        int errors = 0;

        LineIterator it = IOUtils.lineIterator(tsv);
        try {
            for (int lineNum = 1; it.hasNext(); lineNum++) {
                String line = it.nextLine();
                String[] tokens = line.split("\t");

                if (tokens.length != 2) {
                    log.error(String.format("Line %d: Does not contain two tab separated values a project/sample: %s", lineNum, line));
                    errors++;
                    continue;
                }

                String project = tokens[0];
                String sample = tokens[1];
                int version = PicardAggregationUtils.getLatestVersion(project, sample);

                if (version == 0) {
                    log.error(String.format("Line %d: Unable to find a latest version: %s", lineNum, line));
                    errors++;
                    continue;
                }

                picardSamples.add(new PicardSample(project, sample, version));
            }
        } finally {
            it.close();
        }

        if (errors > 0)
            throw new UserException.CouldNotReadInputFile(tsv, String.format("See logger errors for problematic lines."));

        return picardSamples;
    }

    public static PicardIntervals readAnalysisIntervals(List<PicardSample> picardSamples) {
        Set<PicardIntervals> seenIntervals = new LinkedHashSet<PicardIntervals>();

        for (PicardSample picardSample: picardSamples) {
            PicardAnalysisFiles analysis = new PicardAnalysisFiles(picardSample.getProject(), picardSample.getSample(), picardSample.getVersion());
            PicardIntervals picardIntervals = new PicardIntervals(analysis.getReferenceSequence(), analysis.getTargetIntervals());
            seenIntervals.add(picardIntervals);
        }

        if (seenIntervals.isEmpty())
            return null;

        if (seenIntervals.size() == 1)
            return seenIntervals.iterator().next();

        throw new UserException.BadArgumentValue("picardSamples", String.format("%d intervals found: %s", seenIntervals.size(), seenIntervals));
    }

    public static List<String> getSampleBams(List<PicardSample> picardSamples) {
        List<String> bams = new ArrayList<String>();
        for (PicardSample picardSample: picardSamples) {
            bams.add(PicardAggregationUtils.getSampleBam(picardSample.getProject(), picardSample.getSample(), picardSample.getVersion()));
        }
        return bams;
    }

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
     */
    public static String getSampleBam(String project, String sample) {
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
     */
    public static String getSampleDir(String project, String sample) {
        int latestVersion = getLatestVersion(project, sample);
        if (latestVersion == 0)
            throw new UserException.BadArgumentValue("project/sample", "Unable to find a finished directory for project sample " + project + "/" + sample);
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
