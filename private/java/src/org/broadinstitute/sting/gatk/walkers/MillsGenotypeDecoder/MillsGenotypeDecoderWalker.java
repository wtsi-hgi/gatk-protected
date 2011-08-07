package org.broadinstitute.sting.gatk.walkers.MillsGenotypeDecoder;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.MendelianViolation;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by IntelliJ IDEA.
 * User: delangel
 * Date: 7/21/11
 * Time: 3:14 PM
 * To change this template use File | Settings | File Templates.
 */
@Requires(value={},referenceMetaData=@RMD(name="sites", type=VariantContext.class))

public class MillsGenotypeDecoderWalker  extends RodWalker<Integer, Integer> {
    @Output(doc="File to which variants should be written",required=true)
    protected VCFWriter vcfWriter = null;

    @Argument(fullName="gtFile", shortName="gtFile", doc="File with genotype data", required=true)
    private File GT_FILE = new File("");

    private final String variantRodName = "sites";
    private TreeSet<String> samples = new TreeSet<String>();
    private HashMap<String,HashMap<String,Integer>> recordData = new HashMap<String,HashMap<String,Integer>>();

    public void initialize() {
        Pattern samplePattern = Pattern.compile(".*_(NA\\d{5}).*");
        Pattern recordPattern = Pattern.compile("^(\\d+)_.*");

        String[] header;
        ArrayList<String> sampleNames = new ArrayList<String>();
        int numSamples = 0;
        try {
            for ( final String line : new XReadLines( GT_FILE ) ) {
                if (line.startsWith("#")) // headers
                    continue;
                if (line.startsWith("probe")) {
                    header = line.split("\t");
                    numSamples = header.length - 12;
                    for (int k=0; k < numSamples; k++) {
                        // get sample name
                        String fullSampleID = header[k+12];
                        Matcher m = samplePattern.matcher(fullSampleID);
                        if (m.matches()) {
                            String sample = m.group(1);

                            if (!samples.contains(sample))
                                samples.add((sample));

                            // keep also ordered list of names corresponding one to one to columns
                            sampleNames.add(sample);
                        }
                    }
                    final Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
                    hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));

                    vcfWriter.writeHeader(new VCFHeader(hInfo, samples));

                }  else {
                    //decode a record in file
                    String[] record = line.split("\t");
                    Matcher m = recordPattern.matcher(record[0]);

                    if (m.matches()) {
                        String recordID = m.group(1);
                        HashMap<String,Integer> gtMap = new HashMap<String,Integer>();
                        //System.out.println(recordID);
                        // For this Id, build a hashmap of the form Sample -> Genotype
                        for (int k=0; k < sampleNames.size(); k++) {
                            int gt = Integer.parseInt(record[k+12]);
                            // how to handle duplicate samples? ie duplicate columns in file
                            // get sample corresponding to this column
                            String sample = sampleNames.get(k);

                            if (gtMap.containsKey(sample)) {
                                // we've already seen this sample: try to merge genotypes.
                                Integer oldGT = gtMap.get(sample);
                                if (oldGT == -1) // no-call
                                    gtMap.put(sample,gt); // substitute no-call by new value

                            }  else
                                // new eample
                                gtMap.put(sample,gt);
                        }

                        // genotype hash map build: now add map to record map
                        recordData.put(recordID, gtMap);
                    }

                }


            }
        } catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(GT_FILE, e);
        }

    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        Collection<VariantContext> vcs = tracker.getVariantContexts(ref, variantRodName, null, context.getLocation(), true, false);

        if ( vcs == null || vcs.size() == 0) {
            return 0;
        }

        // get ID and see if we've seen the ID in genotype data
        for (VariantContext vc: vcs) {
            String id = vc.getID();
            Set<Allele> alleles = vc.getAlleles();
            if (recordData.containsKey(id)) {
                HashMap<String,Integer> data =  recordData.get(id);
                Map<String,Genotype> genotypes = new HashMap<String,Genotype> (vc.getGenotypes());

                for (String sample : data.keySet()) {
                    Allele a = null,b = null;
                    Integer g =  data.get(sample) ;
                    List<Allele> alleleList = new ArrayList<Allele>();

                    switch(g) {
                        case -1:
                            a = b= Allele.NO_CALL;
                            break;
                        case 0:
                            a = b= vc.getReference();

                            break;
                        case 1:
                            a =  vc.getReference();
                            b = vc.getAlternateAllele(0);
                            break;
                        case 2:
                            a = b= vc.getAlternateAllele(0);
                            break;
                    }
                    alleleList.add(a);
                    alleleList.add(b);
                    Genotype gt = new Genotype(sample,alleleList,0.0);
                    genotypes.put(sample,gt);
                }

                VariantContext vcnew = new VariantContext("GMIlls",vc.getChr(), vc.getStart(), vc.getEnd(), vc.getAlleles(),  genotypes, 99.0, null, null) ;
                vcfWriter.add(vcnew);
            }
        }
        return 1;
    }

    @Override
    public Integer reduceInit() { return 0; }

    @Override
    public Integer reduce(Integer value, Integer sum) { return value + sum; }

    public void onTraversalDone(Integer result) {
        logger.info(result + " records processed.");
    }


}
