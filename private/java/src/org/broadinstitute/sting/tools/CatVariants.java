package org.broadinstitute.sting.tools;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.CloseableTribbleIterator;
import org.broad.tribble.FeatureReader;
import org.broadinstitute.sting.utils.codecs.bcf2.BCF2Codec;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.writer.Options;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriterFactory;
import scala.concurrent.LinkedListQueueCreator;

import java.io.*;
import java.util.*;


/**
 * Created with IntelliJ IDEA.
 * User: ami
 * Date: 11/1/12
 * Time: 7:49 AM
 * To change this template use File | Settings | File Templates.
 */
public class CatVariants {
    // setup the logging system, used by some codecs
    private static org.apache.log4j.Logger logger = org.apache.log4j.Logger.getRootLogger();


    public static void main(String[] args){
        BasicConfigurator.configure();
        logger.setLevel(Level.INFO);

        if (args.length < 3 )
            printUsage();
        boolean sorted = false;
        if (args.length == 4){
            if (args[3].equals("sorted"))
                sorted = true;
            else
                printUsage();
        }


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

        String inputFileName = args[1];
        if (!inputFileName.endsWith(".list")){
            System.err.println("Input File " + inputFileName + " should be <name>.list");
            printUsage();
        }

       // Comparator<Pair<Integer,FeatureReader<VariantContext>>> comparator = new PositionComparator();
        Comparator<Pair<Integer,String>> newComparator = new NewPositionComparator();


        //PriorityQueue<Pair<Integer,FeatureReader<VariantContext>>> queue =
        //        new PriorityQueue<Pair<Integer,FeatureReader<VariantContext>>>(2000, comparator);
        Queue<Pair<Integer,String>> priorityQueue;
        if(sorted)
            priorityQueue = new LinkedList<Pair<Integer,String>>();
        else
            priorityQueue = new PriorityQueue<Pair<Integer,String>>(10000, newComparator);

        try{
            FileInputStream fstream = new FileInputStream(inputFileName);
            DataInputStream in = new DataInputStream(fstream);
            BufferedReader br = new BufferedReader(new InputStreamReader(in));
            String fileName;
            while ((fileName = br.readLine()) != null)   {

                if (!(fileName.endsWith(".vcf") || fileName.endsWith(".VCF") || fileName.endsWith(".bcf") || fileName.endsWith(".BCF"))){
                    System.err.println("File " + fileName + " should be <name>.vcf or <name>.bcf");
                    printUsage();
                }
                if (sorted){
                    priorityQueue.add(new Pair<Integer, String>(0,fileName));
                }
                else{
                    File featureFile = new File(fileName);
                    if (!featureFile.exists()) {
                        System.err.println("File " + featureFile.getAbsolutePath() + " doesn't exist");
                    }
                    FeatureReader<VariantContext> reader;
                    boolean useVCF = (fileName.endsWith(".vcf") || fileName.endsWith(".VCF"));
                    if(useVCF)
                        reader = AbstractFeatureReader.getFeatureReader(featureFile.getAbsolutePath(), new VCFCodec(), false);
                    else
                        reader = AbstractFeatureReader.getFeatureReader(featureFile.getAbsolutePath(), new BCF2Codec(), false);
                    Iterator<VariantContext> it = reader.iterator();
                    if(!it.hasNext()){
                        System.err.println("File " + fileName + " is empty. This file will be ignored");
                        continue;
                    }
                    VariantContext vc = it.next();
                    int firstPosition = vc.getStart();
                    reader.close();
                    //queue.add(new Pair<Integer, FeatureReader<VariantContext>>(firstPosition,reader));
                    priorityQueue.add(new Pair<Integer, String>(firstPosition,fileName));
                }

            }
            in.close();
        }catch (Exception e){
            System.err.println("Error: " + e.getMessage());
            printUsage();
        }

        String outputName = args[2];
        if (!(outputName.endsWith(".vcf") || outputName.endsWith(".VCF"))){
            System.err.println("Output File " + outputName + " should be <name>.vcf");
            printUsage();
        }
        ReferenceSequenceFile ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(refFile);
        final File combinedFile = new File(outputName);

        int counter = 0;
        final int MAX_RECORDS = 10;
        try {
            FileOutputStream outputStream = new FileOutputStream(combinedFile);
            EnumSet<Options> options = EnumSet.of(Options.INDEX_ON_THE_FLY);
            final VariantContextWriter outputWriter = VariantContextWriterFactory.create(combinedFile, outputStream, ref.getSequenceDictionary(), options);

            boolean firstFile = true;
            int count =0;
            //while(!queue.isEmpty()){
            while(!priorityQueue.isEmpty() ){
                count++;
                //FeatureReader<VariantContext> reader = queue.remove().getSecond();
                String fileName = priorityQueue.remove().getSecond();
                File featureFile = new File(fileName);
                if (!featureFile.exists()) {
                    System.err.println("File " + featureFile.getAbsolutePath() + " doesn't exist");
                }
                FeatureReader<VariantContext> reader;
                boolean useVCF = (fileName.endsWith(".vcf") || fileName.endsWith(".VCF"));
                if(useVCF)
                    reader = AbstractFeatureReader.getFeatureReader(featureFile.getAbsolutePath(), new VCFCodec(), false);
                else
                    reader = AbstractFeatureReader.getFeatureReader(featureFile.getAbsolutePath(), new BCF2Codec(), false);

                if(count%10 ==0)
                    System.out.print(count);
                else
                    System.out.print(".");
                if (firstFile){
                    VCFHeader header = (VCFHeader)reader.getHeader();
                    outputWriter.writeHeader(header);
                    firstFile = false;
                }

                Iterator<VariantContext> it = reader.iterator();

                while (it.hasNext()){
                    VariantContext vc = it.next();
                    outputWriter.add(vc);
                    counter++;
                }

                reader.close();

            }
            System.out.println();

            outputStream.close();
            outputWriter.close();


        } catch ( Exception e ) {
            throw new RuntimeException(e);
        }
    }

    /**
     * print usage information
     */
    public static void printUsage() {
        System.err.println("Usage: java -cp dist/GenomeAnalysisTK.jar org.broadinstitute.sting.tools.AppendVariants <reference> <file with list of VCF or BCF files> <outputFile> [sorted (optional)]");
        System.err.println("    the list file must end with .list");
        System.err.println("    Where the files in the list can be of type: VCF (ends in .vcf or .VCF)");
        System.err.println("                                                BCF2 (ends in .bcf or .BCF)");
        System.err.println("    Output file must be vcf or bcf file (.vcf or .bcf)");
        System.err.println("    if the input files are already sorted, the last argument can indicate that");


        /**
         * you could add others here; also look in the GATK code-base for an example of a dynamic way
         * to load Tribble codecs.
         */
        System.exit(1);
    }

    public static AsciiFeatureCodec getFeatureCodec(File featureFile, String variantType) {
        // quickly determine the codec type
        if ( variantType != null ) {
            if (variantType.equals("vcf") ) return new VCFCodec();
            throw new StingException("Explicitly specified rod type "+variantType+" is not recognized");
        }
        if ( featureFile.getName().endsWith(".vcf") || featureFile.getName().endsWith(".VCF") )
            return new VCFCodec();
        // if ( featureFile.getName().endsWith(".bcf") || featureFile.getName().endsWith(".BCF") )
       //     return new BCF2Codec();
        throw new IllegalArgumentException("Unable to determine correct file type based on the file name, for file -> " + featureFile);
    }

    private void readBCF2(File source) throws IOException {
        logger.info("Reading BCF2 from " + source);

        BCF2Codec codec = new BCF2Codec();
        AbstractFeatureReader<VariantContext> featureReader = AbstractFeatureReader.getFeatureReader(source.getAbsolutePath(), codec, false);

        int counter = 0;
        featureReader.getHeader();

        CloseableTribbleIterator<VariantContext> itor = featureReader.iterator();

        //while (itor.hasNext() && (counter++ < MAX_RECORDS || MAX_RECORDS == -1)) {
        //    itor.next();
        //}

        // Not so Closeable...
        //itor.close();
    }

    static class PositionComparator implements Comparator<Pair<Integer,FeatureReader<VariantContext>>> {

        @Override
        public int compare(Pair<Integer,FeatureReader<VariantContext>> p1, Pair<Integer,FeatureReader<VariantContext>> p2) {
            int startPositionP1 = p1.getFirst();
            int startPositionP2 = p2.getFirst();
            if (startPositionP1  == startPositionP2)
                return 0;
            return startPositionP1 < startPositionP2 ? -1 : 1 ;
        }
    }
    static class NewPositionComparator implements Comparator<Pair<Integer,String>> {

        @Override
        public int compare(Pair<Integer,String> p1, Pair<Integer,String> p2) {
            int startPositionP1 = p1.getFirst();
            int startPositionP2 = p2.getFirst();
            if (startPositionP1  == startPositionP2)
                return 0;
            return startPositionP1 < startPositionP2 ? -1 : 1 ;
        }
    }

}
