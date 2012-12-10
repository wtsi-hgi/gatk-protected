package org.broadinstitute.sting.tools;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.CloseableTribbleIterator;
import org.broad.tribble.FeatureReader;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.utils.codecs.bcf2.BCF2Codec;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.writer.Options;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriterFactory;

import java.io.*;
import java.util.*;


/**
 * Created with IntelliJ IDEA.
 * User: ami
 * Date: 11/1/12
 * Time: 7:49 AM
 * To change this template use File | Settings | File Templates.
 */
public class CatVariants extends CommandLineProgram {
    // setup the logging system, used by some codecs
    private static org.apache.log4j.Logger logger = org.apache.log4j.Logger.getRootLogger();

    @Argument(fullName = "reference", shortName = "R", doc = "genome reference file <name>.fasta", required = true)
    private String refarg = null;

    @Argument(fullName = "inputFileName", shortName = "input", doc = "input file name <name>.list", required = true)
    private String inputFileName = null;

    @Argument(fullName = "outputName", shortName = "out", doc = "output file name <name>.vcf or <name>.bcf", required = true)
    private String outputName = null;

    @Argument(fullName = "emailAddress", shortName = "e", doc = "The user email address for which to generate a GATK key", required = true)
    private Boolean sorted = false;

    @Argument(fullName = "help", shortName = "help", doc = "print this hlp info", required = false)
    private Boolean help = false;

    protected int execute() throws Exception {
        if(help){
            printUsage();
            return 0;
        }

        BasicConfigurator.configure();
        logger.setLevel(Level.INFO);

        if ( ! refarg.endsWith(".fasta")) {
            throw new UserException("Reference file "+refarg+"name must end with .fasta");
        }

        File refFile = new File(refarg);
        if ( ! refFile.exists() ) {
            throw new UserException(String.format("Reference file %s does not exist", refFile.getAbsolutePath()));
        }

        if (!inputFileName.endsWith(".list")){
            throw new UserException(String.format("Input File %s should be <name>.list",inputFileName));
        }

        // Comparator<Pair<Integer,FeatureReader<VariantContext>>> comparator = new PositionComparator();
        Comparator<Pair<Integer,String>> newComparator = new PositionComparator();


        //PriorityQueue<Pair<Integer,FeatureReader<VariantContext>>> queue =
        //        new PriorityQueue<Pair<Integer,FeatureReader<VariantContext>>>(2000, comparator);
        Queue<Pair<Integer,String>> priorityQueue;
        if(sorted)
            priorityQueue = new LinkedList<Pair<Integer,String>>();
        else
            priorityQueue = new PriorityQueue<Pair<Integer,String>>(10000, newComparator);


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
                    throw new UserException(String.format("File %s doesn't exist",featureFile.getAbsolutePath()));
                }
                FeatureReader<VariantContext> reader;
                boolean useVCF = (fileName.endsWith(".vcf") || fileName.endsWith(".VCF"));
                if(useVCF)
                    reader = AbstractFeatureReader.getFeatureReader(featureFile.getAbsolutePath(), new VCFCodec(), false);
                else
                    reader = AbstractFeatureReader.getFeatureReader(featureFile.getAbsolutePath(), new BCF2Codec(), false);
                Iterator<VariantContext> it = reader.iterator();
                if(!it.hasNext()){
                    System.err.println(String.format("File %s is empty. This file will be ignored",featureFile.getAbsolutePath()));
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


        if (!(outputName.endsWith(".vcf") || outputName.endsWith(".VCF"))){
            throw new UserException(String.format("Output File ", outputName, " should be <name>.vcf"));
        }
        ReferenceSequenceFile ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(refFile);
        final File combinedFile = new File(outputName);

        FileOutputStream outputStream = new FileOutputStream(combinedFile);
        EnumSet<Options> options = EnumSet.of(Options.INDEX_ON_THE_FLY);
        final VariantContextWriter outputWriter = VariantContextWriterFactory.create(combinedFile, outputStream, ref.getSequenceDictionary(), options);

        boolean firstFile = true;
        int count =0;
        //while(!queue.isEmpty()){
        while(!priorityQueue.isEmpty() ){
            count++;
            //FeatureReader<VariantContext> reader = queue.remove().getSecond();
            fileName = priorityQueue.remove().getSecond();
            File featureFile = new File(fileName);
            if (!featureFile.exists()) {
                throw new UserException(String.format("File ", featureFile.getAbsolutePath(), " doesn't exist"));
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
            }

            reader.close();

        }
        System.out.println();

        outputStream.close();
        outputWriter.close();





        return 0;
    }


    public static void main(String[] args){
        try {
            CatVariants instance = new CatVariants();
            start(instance, args);
            System.exit(CommandLineProgram.result);
        } catch ( UserException e ) {
            printUsage();
            exitSystemWithUserError(e);
        } catch ( Exception e ) {
            exitSystemWithError(e);
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

    static class PositionComparator implements Comparator<Pair<Integer,String>> {

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
