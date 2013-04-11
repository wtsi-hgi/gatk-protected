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

package org.broadinstitute.sting.gatk.walkers.performance;

import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.CloseableTribbleIterator;
import org.broad.tribble.FeatureReader;
import org.broad.tribble.Tribble;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.PositionalBufferedStream;
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
import org.broadinstitute.variant.bcf2.BCF2Codec;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFConstants;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;

import java.io.*;
import java.util.*;

/**
 * Emits specific fields as dictated by the user from one or more VCF files.
 */
public class ProfileRodSystem extends RodWalker<Integer, Integer> {
    @Output(doc="File to which results should be written")
    protected PrintStream out;

    @Input(fullName="vcf", shortName = "vcf", doc="vcf", required=true)
    public RodBinding<VariantContext> vcf;

    @Input(fullName="indexTest", shortName = "indexTest", doc="vcf", required=false)
    public List<File> vcfsForIndexTest = Collections.emptyList();

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
        /** Just test the VCF writing */
        JUST_OUTPUT,
        JUST_BCF2
    }

    SimpleTimer timer = new SimpleTimer("myTimer");

    public void initialize() {
        if ( getToolkit().getIntervals() != null )
            throw new UserException.BadArgumentValue("intervals", "ProfileRodSystem cannot accept intervals");

        if ( profileType == ProfileType.JUST_LOAD_INDICES ) {
            RMDTrackBuilder builder = new RMDTrackBuilder(getToolkit().getReferenceDataSource().getReference().getSequenceDictionary(),
                    getToolkit().getGenomeLocParser(), ValidationExclusion.TYPE.ALL, getToolkit().getArguments().disableAutoIndexCreationAndLockingWhenReadingRods);
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

        if ( profileType == ProfileType.JUST_BCF2 ) {
            testBCF2();
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

            if ( EnumSet.of(ProfileType.ALL, ProfileType.JUST_OUTPUT).contains(profileType) ) {
                out.printf("output.records\t%d\t%.2f%n", i, writeFile(rodFile));
            }
        }

        if ( EnumSet.of(ProfileType.JUST_TRIBBLE, ProfileType.JUST_TRIBBLE_DECODE, ProfileType.JUST_OUTPUT).contains(profileType) )
            System.exit(0);

        timer.start(); // start up timer for map itself
    }

    private void testBCF2() {
        try {
            final File vcfFile = getRodFile();
            final File bcf2File = new File(vcfFile.getName() + ".bcf");
            int counter = 0;
            FeatureReader<VariantContext> reader = AbstractFeatureReader.getFeatureReader(vcfFile.getAbsolutePath(), new VCFCodec(), false);
            FileOutputStream outputStream = new FileOutputStream(bcf2File);
            EnumSet<Options> options = EnumSet.of(Options.FORCE_BCF, Options.INDEX_ON_THE_FLY);
            final VariantContextWriter bcf2Writer = VariantContextWriterFactory.create(bcf2File, outputStream, getToolkit().getReferenceDataSource().getReference().getSequenceDictionary(), options);
            VCFHeader header = GATKVCFUtils.withUpdatedContigs((VCFHeader) reader.getHeader(), getToolkit());
            bcf2Writer.writeHeader(header);

            final List<VariantContext> vcs = new ArrayList<VariantContext>();
            Iterator<VariantContext> it = reader.iterator();
            if ( performanceTest ) {
                logger.info("Beginning performance testing");
                SimpleTimer vcfTimer = new SimpleTimer();
                SimpleTimer bcf2WriterTimer = new SimpleTimer();

                vcfTimer.start();
                while (it.hasNext() && (counter++ < MAX_RECORDS || MAX_RECORDS == -1)) {
                    VariantContext vc = it.next();
                    vc.getNSamples(); // force parsing
                    vcs.add(vc);
                }
                vcfTimer.stop();

                bcf2WriterTimer.start();
                for ( VariantContext vc : vcs ) {
                    // write BCF2 records
                    bcf2Writer.add(vc);
                }
                bcf2Writer.close();
                bcf2WriterTimer.stop();

                logger.info("Read  " + counter + " VCF records in " + vcfTimer.getElapsedTime());
                logger.info("Wrote " + counter + " BCF2 records in " + bcf2WriterTimer.getElapsedTime());
            } else {
                logger.info("Beginning size testing");
                while (it.hasNext() && (counter++ < MAX_RECORDS || MAX_RECORDS == -1)) {
                    VariantContext vc = it.next();
                    bcf2Writer.add(vc);
                }
                bcf2Writer.close();
            }

            if ( performanceTest ) {
                final SimpleTimer bcf2ReaderTimer = new SimpleTimer().start();
                readBCF2(bcf2File);
                logger.info("Read BCF2 in " + bcf2ReaderTimer.getElapsedTime());
            }
        } catch ( Exception e ) {
            throw new RuntimeException(e);
        }
    }

    private void readBCF2(File source) throws IOException {
        logger.info("Reading BCF2 from " + source);

        BCF2Codec codec = new BCF2Codec();
        AbstractFeatureReader<VariantContext> featureReader = AbstractFeatureReader.getFeatureReader(source.getAbsolutePath(), codec, false);

        int counter = 0;
        featureReader.getHeader();

        CloseableTribbleIterator<VariantContext> itor = featureReader.iterator();

        while (itor.hasNext() && (counter++ < MAX_RECORDS || MAX_RECORDS == -1)) {
            itor.next();
        }
        // Not so Closeable...
        //itor.close();
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
                AsciiLineReader lineReader = new AsciiLineReader(new PositionalBufferedStream(s));

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

    private final double writeFile(File f) {

        try {
            FeatureReader<VariantContext> reader = AbstractFeatureReader.getFeatureReader(f.getAbsolutePath(), new VCFCodec(), false);
            VCFHeader header = (VCFHeader)reader.getHeader();
            Iterator<VariantContext> it = reader.iterator();

            ArrayList<VariantContext> VCs = new ArrayList<VariantContext>(10000);

            int counter = 0;
            while ((counter++ < MAX_RECORDS || MAX_RECORDS == -1) && it.hasNext()) {
                VCs.add(it.next());
            }

            // now we start the timer
            timer.start();

            VariantContextWriter writer = VariantContextWriterFactory.create(new File(f.getAbsolutePath() + ".test"), getMasterSequenceDictionary());
            writer.writeHeader(header);
            for ( VariantContext vc : VCs )
                writer.add(vc);
            writer.close();

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