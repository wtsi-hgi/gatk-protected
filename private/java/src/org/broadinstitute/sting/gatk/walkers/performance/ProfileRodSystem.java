/*
 * Copyright (c) 2010, The Broad Institute
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
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrackBuilder;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.SimpleTimer;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.gcf.GCF;
import org.broadinstitute.sting.utils.gcf.GCFHeader;
import org.broadinstitute.sting.utils.gcf.GCFHeaderBuilder;
import org.broadinstitute.sting.utils.gcf.GCFWriter;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.List;

/**
 * Emits specific fields as dictated by the user from one or more VCF files.
 */
public class ProfileRodSystem extends RodWalker<Integer, Integer> {
    @Output(doc="File to which results should be written",required=true)
    protected PrintStream out;

    @Input(fullName="vcf", shortName = "vcf", doc="vcf", required=true)
    public RodBinding<VariantContext> vcf;

    @Input(fullName="indexTest", shortName = "indexTest", doc="vcf", required=false)
    public List<File> vcfsForIndexTest;

    @Argument(fullName="nIterations", shortName="N", doc="Number of raw reading iterations to perform", required=false)
    int nIterations = 1;

    @Argument(fullName="verbose", shortName="verbose", doc="Number of raw reading iterations to perform", required=false)
    boolean VERBOSE = false;

    @Argument(fullName="performanceTest", shortName="performanceTest", doc="Number of raw reading iterations to perform", required=false)
    boolean performanceTest = false;

    @Argument(fullName="maxRecords", shortName="M", doc="Max. number of records to process", required=false)
    int MAX_RECORDS = -1;

    @Argument(fullName="mode", shortName="mode", doc="What kind of profile should we do?", required=false)
    public ProfileType profileType = ProfileType.ALL;

    public enum ProfileType {
        /** Run all tests available */
        ALL,
        /** For profiling loading indices */
        JUST_LOAD_INDICES,
        /** Just test the low-level tribble I/O system */
        JUST_TRIBBLE,
        /** Just test decoding of the VCF file using the low-level tribble I/O system */
        JUST_TRIBBLE_DECODE,
        /** Just test the high-level GATK I/O system */
        JUST_GATK,
        JUST_GVCF
    }

    SimpleTimer timer = new SimpleTimer("myTimer");

    public void initialize() {
        if ( getToolkit().getIntervals() != null )
            throw new UserException.BadArgumentValue("intervals", "ProfileRodSystem cannot accept intervals");

        if ( profileType == ProfileType.JUST_LOAD_INDICES ) {
            RMDTrackBuilder builder = new RMDTrackBuilder(getToolkit().getReferenceDataSource().getReference().getSequenceDictionary(),
                    getToolkit().getGenomeLocParser(), ValidationExclusion.TYPE.ALL);
            int i = 0;
            for ( int iteration = 0; iteration < nIterations; iteration++ ) {
                for ( File x : vcfsForIndexTest ) {
                    try {
                        SimpleTimer gatkLoad = new SimpleTimer().start();
                        builder.loadIndex(x, new VCFCodec());
                        double gatkLoadTime = gatkLoad.getElapsedTime();

                        SimpleTimer tribbleLoad = new SimpleTimer().start();
                        File indexFile = Tribble.indexFile(x);
                        Index index = IndexFactory.loadIndex(indexFile.getAbsolutePath());

                        System.out.printf("GATK load index %d %.2f %.2f%n", ++i, gatkLoadTime, tribbleLoad.getElapsedTime());
                    } catch ( IOException e ) {
                        throw new RuntimeException(e);
                    }
                }
            }
            System.exit(0);
        }

        if ( profileType == ProfileType.JUST_GVCF ) {
            testGVCF();
            System.exit(0);
        }

        File rodFile = getRodFile();

        if ( !EnumSet.of(ProfileType.JUST_GATK).contains(profileType) ) {
            out.printf("# walltime is in seconds%n");
            out.printf("# file is %s%n", rodFile);
            out.printf("# file size is %d bytes%n", rodFile.length());
            out.printf("operation\titeration\twalltime%n");
        }

        for ( int i = 0; i < nIterations; i++ ) {
            if ( EnumSet.of(ProfileType.ALL, ProfileType.JUST_TRIBBLE).contains(profileType) ) {
                out.printf("read.bytes\t%d\t%.2f%n", i, readFile(rodFile, ReadMode.BY_BYTE));
                out.printf("read.line\t%d\t%.2f%n", i, readFile(rodFile, ReadMode.BY_LINE));
                out.printf("line.and.parts\t%d\t%.2f%n", i, readFile(rodFile, ReadMode.BY_PARTS));
                out.printf("decode.loc\t%d\t%.2f%n", i, readFile(rodFile, ReadMode.DECODE_LOC));
            }

            if ( EnumSet.of(ProfileType.ALL, ProfileType.JUST_TRIBBLE, ProfileType.JUST_TRIBBLE_DECODE).contains(profileType) ) {
                out.printf("full.decode\t%d\t%.2f%n", i, readFile(rodFile, ReadMode.DECODE));
            }

        }

        if ( EnumSet.of(ProfileType.JUST_TRIBBLE, ProfileType.JUST_TRIBBLE_DECODE).contains(profileType) )
            System.exit(0);

        timer.start(); // start up timer for map itself
    }

    private final void testGVCF() {
        try {
            final File vcfFile = getRodFile();
            final File gvcfFile = new File(vcfFile.getName() + ".gcf");
            final GCFWriter gcfWriter = new GCFWriter(gvcfFile, false, false);
            int counter = 0;
            final VCFCodec codec = new VCFCodec();
            final AsciiLineReader lineReader = new AsciiLineReader(new FileInputStream(vcfFile));
            VCFHeader header = (VCFHeader)codec.readHeader(lineReader);
            gcfWriter.writeHeader(header);

            final List<VariantContext> vcs = new ArrayList<VariantContext>();
            if ( performanceTest ) {
                logger.info("Beginning performance testing");
                SimpleTimer vcfTimer = new SimpleTimer();
                SimpleTimer gcfWriterTimer = new SimpleTimer();

                vcfTimer.start();
                while (counter++ < MAX_RECORDS || MAX_RECORDS == -1) {
                    String line = lineReader.readLine();
                    if ( line == null )
                        break;

                    VariantContext vc = (VariantContext)codec.decode(line);
                    vc.getNSamples(); // force parsing
                    vcs.add(vc);
                }
                vcfTimer.stop();

                gcfWriterTimer.start();
                for ( VariantContext vc : vcs ) {
                    // write GCF records
                    gcfWriter.add(vc);
                }
                gcfWriter.close();
                gcfWriterTimer.stop();

                logger.info("Read  " + counter + " VCF records in " + vcfTimer.getElapsedTime());
                logger.info("Wrote " + counter + " GCF records in " + gcfWriterTimer.getElapsedTime());
            } else {
                logger.info("Beginning size testing");
                while (counter++ < MAX_RECORDS || MAX_RECORDS == -1) {
                    String line = lineReader.readLine();
                    if ( line == null )
                        break;

                    VariantContext vc = (VariantContext)codec.decode(line);
                    gcfWriter.add(vc);
                }
                gcfWriter.close();
            }

            if ( performanceTest ) {
                for ( boolean skipGenotypes : Arrays.asList(false, true) ) {
                    final SimpleTimer gcfReaderTimer = new SimpleTimer().start();
                    readGVCF(gvcfFile, vcs, skipGenotypes);
                    logger.info("Read GVCF in " + gcfReaderTimer.getElapsedTime() + " skipGenotypes = " + skipGenotypes);
                }
            }
        } catch ( Exception e ) {
            throw new RuntimeException(e);
        }
    }

    private void readGVCF(File source, List<VariantContext> vcs, boolean skipGenotypes) throws IOException {
        FileInputStream fileInputStream = GCF.createFileInputStream(source);
        DataInputStream inputStream = GCF.createDataInputStream(fileInputStream);
        logger.info("Reading GVCF from " + source);
        GCFHeader GCFHeader = new GCFHeader(fileInputStream);

        try {
            if ( ! vcs.isEmpty() ) {
                for ( VariantContext vc : vcs ) {
                    if ( VERBOSE ) logger.info("Original VCF: " + vc);
                    GCF GCF = new GCF(inputStream, skipGenotypes);
                    VariantContext decoded = GCF.decode("gcf", GCFHeader);
                    //logger.info("GVCF        : " + gcf);
                    if ( VERBOSE ) logger.info("GVCF -> VCF : " + decoded);
                }
            } else {
                while ( true ) {
                    GCF GCF = new GCF(inputStream, skipGenotypes);
                    VariantContext decoded = GCF.decode("gcf", GCFHeader);
                    //logger.info("GVCF        : " + gcf);
                    if ( VERBOSE ) logger.info("GVCF -> VCF : " + decoded);
                }
            }
        } catch ( EOFException e ) {
            ; // done reading
        }

        inputStream.close();
    }

private enum ReadMode { BY_BYTE, BY_LINE, BY_PARTS, DECODE_LOC, DECODE };

private final double readFile(File f, ReadMode mode) {
        timer.start();

try {
        byte[] data = new byte[100000];
FileInputStream s = new FileInputStream(f);

if ( mode == ReadMode.BY_BYTE ) {
        while (true) {
        if ( s.read(data) == -1 )
        break;
}
        } else {
        int counter = 0;
VCFCodec codec = new VCFCodec();
String[] parts = new String[100000];
AsciiLineReader lineReader = new AsciiLineReader(s);

if ( mode == ReadMode.DECODE_LOC || mode == ReadMode.DECODE )
        codec.readHeader(lineReader);

while (counter++ < MAX_RECORDS || MAX_RECORDS == -1) {
        String line = lineReader.readLine();
if ( line == null )
        break;
else if ( mode == ReadMode.BY_PARTS ) {
        ParsingUtils.split(line, parts, VCFConstants.FIELD_SEPARATOR_CHAR);
}
        else if ( mode == ReadMode.DECODE_LOC ) {
        codec.decodeLoc(line);
}
        else if ( mode == ReadMode.DECODE ) {
        processOneVC((VariantContext)codec.decode(line));
}
        }
        }
        } catch ( Exception e ) {
        throw new RuntimeException(e);
}

        return timer.getElapsedTime();
}

private File getRodFile() {
        List<ReferenceOrderedDataSource> rods = this.getToolkit().getRodDataSources();
ReferenceOrderedDataSource rod = rods.get(0);
return rod.getFile();
}

public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
        return 0;

VariantContext vc = tracker.getFirstValue(vcf, context.getLocation());
if ( vc != null )
        processOneVC(vc);

return 0;
}

private static final void processOneVC(VariantContext vc) {
        vc.getNSamples(); // force us to parse the samples
}

public Integer reduceInit() {
        return 0;
}

public Integer reduce(Integer counter, Integer sum) {
        return counter + sum;
}

public void onTraversalDone(Integer sum) {
        out.printf("gatk.traversal\t%d\t%.2f%n", 0, timer.getElapsedTime());
}
        }