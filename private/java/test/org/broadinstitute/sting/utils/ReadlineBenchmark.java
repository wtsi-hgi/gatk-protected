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

package org.broadinstitute.sting.utils;

import com.google.caliper.Param;
import com.google.caliper.SimpleBenchmark;
import com.google.caliper.runner.CaliperMain;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.*;

/**
 * Caliper microbenchmark of reading lines from a VCF file
 */
public class ReadlineBenchmark extends SimpleBenchmark {

    private enum LineReaderAlgorithm {
        DEFAULT_ASCII_READER,
        ASCII_READER_WITH_ADAPTIVE_BUFFER_SIZE
    }

    private final String OmniGenotypesFile = "/humgen/gsa-hpprojects/GATK/bundle/current/b37/1000G_omni2.5.b37.vcf";
    private final int MAX_LINES = 1000;

    @Param({"1", "100", "1000"})
    int linesToRead; // set automatically by framework

    @Param({"512000", "1024000", "2048000"})
    int bufferSize; // set automatically by framework

    private String INPUT_STRING;

    @Override protected void setUp() {

        // read it into a String so that we don't try to benchmark IO issues
        try {
            FileInputStream s = new FileInputStream(new File(OmniGenotypesFile));
            AsciiLineReader lineReader = new AsciiLineReader(s);
            int counter = 0;
            StringBuffer sb = new StringBuffer();
            while (counter++ < MAX_LINES ) {
                String line = lineReader.readLine();
                if ( line == null )
                    break;
                sb.append(line + "\n");
            }
            s.close();
            INPUT_STRING = sb.toString();

        } catch (Exception e) {
            System.out.println("Benchmarking setup failure because of " + e.getMessage());
        }
    }

    private void run(int rep, LineReaderAlgorithm algorithm) {
        for ( int i = 0; i < rep; i++ )
            readLines(algorithm);
    }

    private void readLines(LineReaderAlgorithm algorithm) {
        try {
            InputStream is = new ByteArrayInputStream(INPUT_STRING.getBytes());
            AsciiLineReader lineReader;

            if ( algorithm == LineReaderAlgorithm.DEFAULT_ASCII_READER )
                lineReader = new AsciiLineReader(is, bufferSize);
            else // if ( algorithm == LineReaderAlgorithm.ASCII_READER_WITH_ADAPTIVE_BUFFER_SIZE )
                // TODO -- not implemented yet, so just use the default reader for now
                lineReader = new AsciiLineReader(is, bufferSize);

            int counter = 0;
            while (counter++ < linesToRead ) {
                String line = lineReader.readLine();
                if ( line == null )
                    break;
            }
        } catch (Exception e) {
            System.out.println("Benchmarking run failure because of " + e.getMessage());
        }
    }

    public void timeDefault(int rep) {
        run(rep, LineReaderAlgorithm.DEFAULT_ASCII_READER);
    }

    public void timeAdaptive(int rep) {
        run(rep, LineReaderAlgorithm.ASCII_READER_WITH_ADAPTIVE_BUFFER_SIZE);
    }

    public static void main(String[] args) {
        CaliperMain.main(ReadlineBenchmark.class, args);
    }
}
