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

package org.broadinstitute.sting.tools;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.util.SortingCollection;
import org.apache.log4j.BasicConfigurator;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.Feature;
import org.broad.tribble.bed.BEDCodec;
import org.broad.tribble.dbsnp.OldDbSNPCodec;
import org.broad.tribble.gelitext.GeliTextCodec;
import org.broad.tribble.readers.LineIterator;
import org.broad.tribble.readers.PositionalBufferedStream;
import org.broadinstitute.sting.gatk.features.maf.MafCodec;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.codecs.CGVarCodec;
import org.broadinstitute.sting.utils.codecs.SoapSNPCodec;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.variant.vcf.VCFCodec;

import java.io.*;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Jan 28, 2011
 * Time: 12:15:03 PM
 * To change this template use File | Settings | File Templates.
 */
public class SortROD {
    // setup the logging system, used by some codecs
    private static org.apache.log4j.Logger logger = org.apache.log4j.Logger.getRootLogger();

    /**
     * this class:
     *  1) checks to see that the feature file exists
     *  2) loads an index from disk, if one doesn't exist, it creates it and writes it to disk
     *  3) creates a FeatureSource
     *  4) iterates over the records, emitting a final tally for the number of features seen
     *
     * @param args a single parameter, the file name to load
     */
    public static void main(String[] args) throws IOException {
        BasicConfigurator.configure();
        logger.setLevel(org.apache.log4j.Level.INFO);
        // check yourself before you wreck yourself - we require one arg, the input file
        if (args.length != 3 )
            printUsage();


        String refarg = args[0];
        if ( ! refarg.endsWith(".fasta")) {
            System.err.println("Reference file name must end with .fasta");
            System.exit(1);
        }

        File refFile = new File(refarg);
        if ( ! refFile.exists() ) {
            System.err.println("Reference file "+refarg+" does not exist");
            System.exit(1);
        }

        String rodType = null;
        String inputArg;
        // our feature file
        int pos = args[1].indexOf(":");
        if ( pos == -1 ) {
            inputArg = args[1];
        } else {
            rodType = args[1].substring(0,pos);
            inputArg = args[1].substring(pos+1);
        }
        File featureFile = new File(inputArg);
        if (!featureFile.exists()) {
            System.err.println("File " + featureFile.getAbsolutePath() + " doesn't exist");
            printUsage();
        }

        BufferedWriter out = null;
        try {
            out = new BufferedWriter(new FileWriter(args[2]));
        } catch ( IOException e ) {
            System.err.println("Can not open output file "+args[2]+" for writing");
            System.exit(1);
        }

        // determine the codec
        AsciiFeatureCodec featureCodec = getFeatureCodec(featureFile,rodType);
        ReferenceSequenceFile ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(refFile);

        LineIterator source = null;
        try {
            source = featureCodec.makeSourceFromStream(new PositionalBufferedStream(new FileInputStream(featureFile)));
        } catch (FileNotFoundException e) {
            System.err.println("File "+featureFile.getAbsolutePath()+" doesn't exist");
            System.exit(1);
        }

        // read the headers
        featureCodec.readHeader(source);

        GenomeLocParser parser = new GenomeLocParser(ref.getSequenceDictionary());

        SortingCollection<String> sorter = SortingCollection.newInstance(String.class,
                new LineCodec(),
                new FeatureComparator(featureCodec,parser),200000);

        int nLines = 0;
        try {
            String currentLine = source.next();
            while ( currentLine != null ) {
                nLines++;

                // uncomment if null returns should be ignored
                //if ( featureCodec.decodeLoc(currentLine) != null )
                sorter.add(currentLine);

                currentLine = source.next();
            }
            
            for ( String s : sorter ) {
                out.write(s);
                out.write('\n');
            }
            out.close();
        } catch (IOException e) {
            System.err.println("Writing failed to the output file "+args[2]);
            System.exit(1);
        }

        logger.info("Sorting finished. Processed lines: "+nLines);
   //     runWithIndex(featureFile, codec, optimizeIndex);

    }


    /**
     * print usage information
     */
    public static void printUsage() {
        System.err.println("Usage: java -jar SortROD.jar <reference> [<rodType>:]<inputFile> <outputFile>");
        System.err.println("    Where input can be of type: VCF (ends in .vcf or .VCF)");
        System.err.println("                                Bed (ends in .bed or .bed)");
        System.err.println("                                DbSNP (ends in .snp or .rod)");
        System.err.println("                                MAF (ends in .maf)");
        System.err.println("    If input file has non-standard extension, rodType can be specified");
        System.err.println("    (rodType always takes precedence over file extension, even if the");
        System.err.println("    latter is otherwise recognizable). rodType can be vcf, bed, dbsnp, or maf");
        System.err.println("    Reference is what the input file needs to be sorted against");

        /**
         * you could add others here; also look in the GATK code-base for an example of a dynamic way
         * to load Tribble codecs.
         */
        System.exit(1);
    }


    public static AsciiFeatureCodec getFeatureCodec(File featureFile, String rodType) {
        // quickly determine the codec type
        if ( rodType != null ) {
            if (rodType.equals("vcf") ) return new VCFCodec();
            if (rodType.equals("bed") ) return new BEDCodec();
            if (rodType.equals("cgvar") || rodType.equals("CGVar") ) return new CGVarCodec();
            if (rodType.equals("snp") || rodType.equals("dbsnp") ) return new OldDbSNPCodec();
            if (rodType.equals("geli.calls") || rodType.equals("geli") ) return new GeliTextCodec();
            if (rodType.equals("txt") ) return new SoapSNPCodec();
            if (rodType.equals("maf") ) return new MafCodec();
            throw new StingException("Explicitly specified rod type "+rodType+" is not recognized");
        }
        if ( featureFile.getName().endsWith(".vcf") || featureFile.getName().endsWith(".VCF") )
            return new VCFCodec();
        if (featureFile.getName().endsWith(".bed") || featureFile.getName().endsWith(".BED") )
            return new BEDCodec();
        if ( featureFile.getName().endsWith(".tsv") || featureFile.getName().endsWith(".TSV") )
            return new CGVarCodec();
        if (featureFile.getName().endsWith(".snp") || featureFile.getName().endsWith(".rod") )
            return new OldDbSNPCodec();
        if (featureFile.getName().endsWith(".geli.calls") || featureFile.getName().endsWith(".geli") )
            return new GeliTextCodec();
        if (featureFile.getName().endsWith(".txt") || featureFile.getName().endsWith(".TXT") )
            return new SoapSNPCodec();
        if (featureFile.getName().endsWith(".maf") || featureFile.getName().endsWith(".MAF") )
            return new MafCodec();
        throw new IllegalArgumentException("Unable to determine correct file type based on the file name, for file -> " + featureFile);
    }

    static class LineCodec implements SortingCollection.Codec<String> {
        OutputStream os;
        InputStream is;

        public void setOutputStream(OutputStream outputStream) {
            os = outputStream;
        }

        public void setInputStream(InputStream inputStream) {
            is = inputStream;
        }

        public void encode(String s) {
            try {
                os.write(s.getBytes());
                os.write('\n');
            } catch (IOException e) {
                throw new StingException("SortingCollection: Write into temporary file failed",e);
            }
        }

        public String decode() {
            List<Byte> l = new ArrayList<Byte>(1024);
            try {
                int c = is.read();
                while ( c != -1 && c != '\n' ) {
                    l.add((byte)c);
                    c = is.read();
                }
            } catch (IOException e) {
                throw new StingException("SortingCollection: Read from temporary file failed",e);
            }
            return new String(toByteArray(l));  //To change body of implemented methods use File | Settings | File Templates.
        }

        public SortingCollection.Codec<String> clone() {
            LineCodec codec = new LineCodec();
            codec.setInputStream(is);
            codec.setOutputStream(os);
            return codec;  //To change body of implemented methods use File | Settings | File Templates.
        }

        private byte [] toByteArray(List<Byte> l) {
            byte[] ret = new byte[l.size()];
            for ( int i = 0 ; i < l.size() ; i++ ) ret[i] = l.get(i);
            return ret;
        }
    }

    static class FeatureComparator implements Comparator<String> {
        AsciiFeatureCodec codec ;
        GenomeLocParser parser;

        public FeatureComparator (AsciiFeatureCodec codec, GenomeLocParser parser) {
            this.codec = codec;
            this.parser=parser;
        }

        /**
         * Compares its two arguments for order.  Returns a negative integer,
         * zero, or a positive integer as the first argument is less than, equal
         * to, or greater than the second.<p>
         * <p/>
         * In the foregoing description, the notation
         * <tt>sgn(</tt><i>expression</i><tt>)</tt> designates the mathematical
         * <i>signum</i> function, which is defined to return one of <tt>-1</tt>,
         * <tt>0</tt>, or <tt>1</tt> according to whether the value of
         * <i>expression</i> is negative, zero or positive.<p>
         * <p/>
         * The implementor must ensure that <tt>sgn(compare(x, y)) ==
         * -sgn(compare(y, x))</tt> for all <tt>x</tt> and <tt>y</tt>.  (This
         * implies that <tt>compare(x, y)</tt> must throw an exception if and only
         * if <tt>compare(y, x)</tt> throws an exception.)<p>
         * <p/>
         * The implementor must also ensure that the relation is transitive:
         * <tt>((compare(x, y)&gt;0) &amp;&amp; (compare(y, z)&gt;0))</tt> implies
         * <tt>compare(x, z)&gt;0</tt>.<p>
         * <p/>
         * Finally, the implementor must ensure that <tt>compare(x, y)==0</tt>
         * implies that <tt>sgn(compare(x, z))==sgn(compare(y, z))</tt> for all
         * <tt>z</tt>.<p>
         * <p/>
         * It is generally the case, but <i>not</i> strictly required that
         * <tt>(compare(x, y)==0) == (x.equals(y))</tt>.  Generally speaking,
         * any comparator that violates this condition should clearly indicate
         * this fact.  The recommended language is "Note: this comparator
         * imposes orderings that are inconsistent with equals."
         *
         * @param o1 the first object to be compared.
         * @param o2 the second object to be compared.
         * @return a negative integer, zero, or a positive integer as the
         *         first argument is less than, equal to, or greater than the
         *         second.
         * @throws ClassCastException if the arguments' types prevent them from
         *                            being compared by this comparator.
         */
        public int compare(String o1, String o2) {
            Feature f1 = null, f2 = null;
            try {
                f1 = codec.decodeLoc(o1);
                f2 = codec.decodeLoc(o2);
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
            
            if ( f1 == null ) {
                if ( f2 == null ) return 0;
                else return -1; // null is less than non-null, this will hopefully push header strings up (but commented out lines will move up too!)
            }
            // f1 is not null
            if ( f2 == null ) return 1;

            GenomeLoc l1 = parser.createGenomeLoc(f1.getChr(),f1.getStart(),f1.getEnd());
            GenomeLoc l2 = parser.createGenomeLoc(f2.getChr(),f2.getStart(),f2.getEnd());
            return l1.compareTo(l2);  //To change body of implemented methods use File | Settings | File Templates.
        }
    }
}

