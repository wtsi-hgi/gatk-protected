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

package org.broadinstitute.sting.gatk.walkers.concordance;

import org.broad.tribble.vcf.VCFGenotypeRecord;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.vcf.VCFRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.classloader.PackageUtils;
import org.broadinstitute.variant.utils.Pair;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.genotype.vcf.*;

import java.io.File;
import java.util.*;


/**
 * Determines the concordance between multiple VCF call sets at each position.
 * Users can specify which concordance tests should be run.
 */
@Requires(value={DataSource.REFERENCE})
@Reference(window=@Window(start=-20,stop=20))
public class CallsetConcordanceWalker extends RodWalker<Integer, Integer> {
    @Argument(fullName="concordance_output", shortName="CO", doc="VCF file to which output should be written", required=true)
    private File OUTPUT = null;
    @Argument(fullName="concordanceType", shortName="CT", doc="Concordance subset types to apply to given callsets.   Syntax: 'type[:key1=arg1,key2=arg2,...]'", required=false)
    private String[] TYPES = null;
    @Argument(fullName="list", shortName="ls", doc="List the available concordance types and exit", required=false)
    private Boolean LIST_ONLY = false;


    // the concordance tests to run
    private ArrayList<ConcordanceType> requestedTypes;

    // VCF writer for the output of the concordance tests
    private VCFWriter vcfWriter;

    // a map of rod name to uniquified sample name
    private HashMap<Pair<String, String>, String> rodNamesToSampleNames = new HashMap<Pair<String, String>, String>();


    /**
     * Prepare the output file and the list of available features.
     */
    public void initialize() {

        // get the possible concordance types
        List<Class<? extends ConcordanceType>> classes = PackageUtils.getClassesImplementingInterface(ConcordanceType.class);

        // print and exit if that's what was requested
        if ( LIST_ONLY ) {
            out.println("\nAvailable concordance types:");
            for (int i = 0; i < classes.size(); i++)
                out.println("\t" + classes.get(i).getSimpleName());
            out.println();
            System.exit(0);
        }

        // get the list of all sample names from the various input rods (they need to be uniquified in case there's overlap)
        HashSet<String> samples = new HashSet<String>();
        SampleUtils.getUniquifiedSamplesFromRods(getToolkit(), samples, rodNamesToSampleNames);

        for ( java.util.Map.Entry<Pair<String, String>, String> entry : rodNamesToSampleNames.entrySet() ) {
            logger.debug("Uniquified sample mapping: " + entry.getKey().first + "/" + entry.getKey().second + " -> " + entry.getValue());
        }

        // initialize requested concordance types
        requestedTypes = new ArrayList<ConcordanceType>();
        if (TYPES != null) {
            for ( String requestedTypeString : TYPES ) {
                String[] requestedPieces = requestedTypeString.split(":");
                String requestedType = requestedPieces[0];

                boolean foundClass = false;
                for ( Class type : classes ) {

                    if (requestedType.equalsIgnoreCase(type.getSimpleName())) {
                        foundClass = true;
                        try {
                            ConcordanceType concordance = (ConcordanceType)type.newInstance();
                            HashMap<String,String> requestedArgs = new HashMap<String,String>();
                            if ( requestedPieces.length == 2 ) {
                                String[] argStrings = requestedPieces[1].split(",");
                                for (int i = 0; i < argStrings.length; i++ ) {
                                    String[] arg = argStrings[i].split("=");
                                    if ( arg.length == 2 )
                                        requestedArgs.put(arg[0], arg[1]);
                                }
                            }

                            concordance.initialize(requestedArgs, samples);
                            requestedTypes.add(concordance);
                            break;
                        } catch (InstantiationException e) {
                            throw new StingException(String.format("Cannot instantiate concordance class '%s': must be concrete class", type.getSimpleName()));
                        } catch (IllegalAccessException e) {
                            throw new StingException(String.format("Cannot instantiate concordance class '%s': must have no-arg constructor", type.getSimpleName()));
                        }
                    }
                }

                if ( !foundClass )
                    throw new StingException("The requested concordance type (" + requestedType + ") isn't a valid concordance option");
            }
        }

        // set up the header fields
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("source", "CallsetConcordance"));
        hInfo.add(new VCFHeaderLine("note", "\"This file represents a concordance test of various call sets - NOT the output from a multi-sample caller\""));
        hInfo.addAll(getVCFAnnotationDescriptions(requestedTypes));

        vcfWriter = new VCFWriter(OUTPUT);
        vcfWriter.writeHeader(new VCFHeader(hInfo, samples));
    }

    public static Set<VCFHeaderLine> getVCFAnnotationDescriptions(Collection<ConcordanceType> types) {

        TreeSet<VCFHeaderLine> descriptions = new TreeSet<VCFHeaderLine>();
        for ( ConcordanceType type : types )
            descriptions.add(type.getInfoDescription());

        return descriptions;
    }

    public Integer map(RefMetaDataTracker rodData, ReferenceContext ref, AlignmentContext context) {
        if ( rodData == null ) // RodWalkers can make funky map calls
            return 0;

        // get all of the vcf rods at this locus
        Map<VCFRecord,String> vcfRods = new LinkedHashMap<VCFRecord,String>();
        Iterator<GATKFeature> rods = rodData.getAllRods().iterator();
        while (rods.hasNext()) {
            GATKFeature rod = rods.next();
            if ( rod.getUnderlyingObject() instanceof VCFRecord ) {
                if (vcfRods.containsKey(rod)) throw new StingException("Duplicate VCF's found");
                vcfRods.put((VCFRecord)rod.getUnderlyingObject(),rod.getName());
            }
        }

        if ( vcfRods.size() == 0 )
            return 0;

        // pull out all of the individual calls from the rods and insert into a map based on the
        // mapping from rod/sample to uniquified name
        HashMap<String, VCFGenotypeRecord> samplesToRecords = new HashMap<String, VCFGenotypeRecord>();
        for ( VCFRecord rod : vcfRods.keySet() ) {
            List<VCFGenotypeRecord> records = rod.getVCFGenotypeRecords();
            for ( VCFGenotypeRecord vcfRec : records ) {
                String uniquifiedSample = rodNamesToSampleNames.get(new Pair<String, String>(vcfRods.get(rod), vcfRec.getSampleName()));
                if ( uniquifiedSample == null )
                    throw new StingException("Unexpected sample encountered: " + vcfRec.getSampleName() + " in rod " + vcfRods.get(rod));

                samplesToRecords.put(uniquifiedSample, vcfRec);
            }
        }

        // create a merged record from all input VCFs
        VCFRecord record = VCFUtils.mergeRecords(vcfRods, rodNamesToSampleNames);

        // add in the info fields to the new record based on the results of each of the relevant concordance tests
        for ( ConcordanceType type : requestedTypes ) {
            String result = type.computeConcordance(samplesToRecords, ref);
            if ( result != null ) {
                record.addInfoField(type.getInfoName(), result);
            }
        }

        // emit the new record
        vcfWriter.addRecord(record);

        return 1;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public void onTraversalDone(Integer result) {
        vcfWriter.close();
        out.printf("Processed %d loci.\n", result);
    }
}
