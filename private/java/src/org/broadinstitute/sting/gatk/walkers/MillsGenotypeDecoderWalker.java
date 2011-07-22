package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.MendelianViolation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Collection;
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

    public void initialize() {
       // Pattern samplePattern = Pattern.compile("\w+_(\w+)\.\w+");

        try {
            for ( final String line : new XReadLines( GT_FILE ) ) {
                if (line.startsWith("#")) // headers
                    continue;
                if (line.startsWith("probe")) {
                    String[] header = line.split("\t");
                    int numSamples = header.length - 12;
                    for (int k=0; k < numSamples; k++) {
                        // get sample name
                        //Matcher m = samplePattern.matcher(header[k+12]);
                        String sample = header[k+12].replaceAll("\\w+_","");
                        System.out.format("%s %s\n", header[k+12], sample);
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
