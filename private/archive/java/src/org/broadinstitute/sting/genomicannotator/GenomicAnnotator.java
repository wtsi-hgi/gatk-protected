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

package org.broadinstitute.sting.gatk.walkers.annotator.genomicannotator;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.features.annotator.AnnotatorInputTableCodec;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.variant.vcf.VCFWriter;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.Map.Entry;

/**
 * Annotates variant calls with information from user-specified tabular files.
 *
 * For details, see:  http://www.broadinstitute.org/gsa/wiki/index.php/GenomicAnnotator
 */
@Requires(value={DataSource.REFERENCE},referenceMetaData=@RMD(name="variant",type=VariantContext.class))
@By(DataSource.REFERENCE)
public class GenomicAnnotator extends RodWalker<Integer, Integer> implements TreeReducible<Integer> {

    @Output(doc="File to which variants should be written",required=true)
    protected VCFWriter vcfWriter = null;

    @Argument(fullName="vcfOutput", shortName="vcf", doc="Please use --out instead", required=false)
    @Deprecated
    protected String oldOutArg;

    @Argument(fullName="sampleName", shortName="sample", doc="The sample (NA-ID) corresponding to the variant input (for non-VCF input only)", required=false)
    protected String sampleName = null;

    @Argument(fullName="select", shortName="s", doc="Optionally specifies which subset of columns from which -B inputs should be used for annotations. For example, -B:mydbsnp,AnnotatorInputTable /path/to/mydbsnp.txt -B:mytable,AnnotatorInputTable /path/mytable.txt -s mydbsnp.avHet,mydbsnp.name,mytable.column3 will cause annotations to only be generated from the 3 columns specified using -s.", required=false)
    protected String[] SELECT_COLUMNS = {};

    @Argument(fullName="join", shortName="J", doc="Optionally specifies a file and column within that file that should be LEFT-JOIN'ed to a column in a previously-specified file. The file provided to -J must be tab-delimited, with the first non-comment/non-empty line containing column names. (example: -B:name,AnnotatorInputTable /path/to/file1   -J name2,/path/to/file2,name.columnName=name2.columnName2  - this will join the table in file2 to the table in file1) ", required=false)
    protected String[] JOIN_ARGS = {};

    @Argument(fullName="oneToMany", shortName="m", doc="If more than one record from the same file matches a particular locus (for example, multiple dbSNP records with the same position), create multiple entries in the ouptut VCF file - one for each match. If a particular tabular file has J matches, and another tabular file has K matches for a given locus, then J*K output VCF records will be generated - one for each pair of K, J.   If this flag is not provided, the multiple records are still generated, but they are stored in the INFO field of a single output VCF record, with their annotation keys differentiated by appending '_i' with i varying from 1 to K*J. ", required=false)
    protected Boolean ONE_TO_MANY = false;

    @Argument(fullName="maxJoinTableSize", shortName="maxJoin", doc="The maximum allowed size (i.e. number of rows) for a table provided with the -J argument", required=false)
    protected Integer MAX_JOIN_TABLE_SIZE = 500000;

    @Argument(fullName="ignoreFilteredSites", shortName="noFilt", doc="If specified, don't annotate sites marked as filtered out")
    protected Boolean IGNORE_FILTERED_SITES = false;

    private VariantAnnotatorEngine engine;

    /**
     * Prepare the output file and the list of available features.
     */
    public void initialize() {

        //read all ROD file headers and construct a set of all column names to be used for validation of command-line args
        final Set<String> allFullyQualifiedColumnNames = new LinkedHashSet<String>();
        final Set<String> allBindingNames = new LinkedHashSet<String>();
            for(ReferenceOrderedDataSource ds : getToolkit().getRodDataSources()) {
                if(! ds.getType().equals(AnnotatorInputTableCodec.class)) {
                    continue; //skip all non-AnnotatorInputTable files.
                }
                final String bindingName = ds.getName();
                File file = ds.getFile();
                allBindingNames.add(bindingName);
                try {
                    final ArrayList<String> header = AnnotatorInputTableCodec.readHeader(file);
                    for(String columnName : header) {
                        allFullyQualifiedColumnNames.add(bindingName + "." + columnName);
                    }
                } catch(IOException e) {
                    throw new UserException.CouldNotReadInputFile(file, "Failed when attempting to read file header. ", e);
                }
            }

        //parse the JOIN_COLUMNS args, read in the specified files, and validate column names in the = relation. This end result of this loop is to populate the List of joinTables with one entry per -J arg.
        final List<JoinTable> joinTables = new LinkedList<JoinTable>();
        for(String joinArg : JOIN_ARGS) {

            //parse the tokens
            final String[] arg = joinArg.split(",");
            if(arg.length != 3) {
                throw new UserException.BadArgumentValue("-J", "The following -J arg: \"" + joinArg + "\" must contain 3 comma-separated values. (ex: -J name,/path/to/file,name.columnName=name2.columnName2)");
            }
            final String bindingName = arg[0];
            final String filename = arg[1];
            final String columnsToJoin = arg[2];

            if(allBindingNames.contains(bindingName)) {
                throw new UserException.BadArgumentValue("-J", "The name \"" + bindingName + "\" in the -J arg: \"" + joinArg + "\" has already been used in another binding.");
            }

            String[] splitOnEquals = columnsToJoin.split("=+");
            if(splitOnEquals.length != 2) {
                throw new UserException.BadArgumentValue("-J", "The -J arg: \"" + joinArg + "\" must specify the columns to join on. (ex: -J name,/path/to/file,name.columnName=name2.columnName2)");
            }

            String[] splitOnDot1 = splitOnEquals[0].split("\\.");
            String[] splitOnDot2 = splitOnEquals[1].split("\\.");
            if(splitOnDot1.length != 2 || splitOnDot2.length != 2) {
                throw new UserException.BadArgumentValue("-J", "The -J arg: \"" + joinArg + "\" must fully specify the columns to join on. (ex: -J name,/path/to/file,name.columnName=name2.columnName2)");
            }

            final String bindingName1 = splitOnDot1[0];
            final String columnName1 = splitOnDot1[1];
            final String bindingName2 = splitOnDot2[0];
            final String columnName2 = splitOnDot2[1];

            //figure out which of the 2 binding names within the = relation matches the -J bindingName
            final String localBindingName = bindingName; //alias
            final String localColumnName;
            final String externalBindingName;
            final String externalColumnName;
            if(bindingName1.equals(bindingName)) {
                localColumnName = columnName1;
                externalBindingName = bindingName2;
                externalColumnName = columnName2;
            } else if(bindingName2.equals(bindingName)) {
                localColumnName = columnName2;
                externalBindingName = bindingName1;
                externalColumnName = columnName1;
            } else {
                throw new UserException.BadArgumentValue("-J", "The name \"" + bindingName + "\" in the -J arg: \"" + joinArg + "\" must be specified in one the columns to join on. (ex: -J name,/path/to/file,name.columnName=name2.columnName2)");
            }

            //validate externalColumnName
            final String fullyQualifiedExternalColumnName = externalBindingName + '.' + externalColumnName;
            if( !allFullyQualifiedColumnNames.contains(fullyQualifiedExternalColumnName) ) {
                throw new UserException.BadArgumentValue("-J", "The -J arg: \"" + joinArg + "\" specifies an unknown column name: \"" + fullyQualifiedExternalColumnName + "\"");
            }

            //read in the file contents into a JoinTable object
            final JoinTable joinTable = new JoinTable(MAX_JOIN_TABLE_SIZE);
            joinTable.parseFromFile(filename, localBindingName, localColumnName, externalBindingName, externalColumnName);
            joinTables.add(joinTable);

            //validate localColumnName, and add all column names in this file to the list of allFullyQualifiedColumnNames so that they can be referenced from subsequent -J args.
            final List<String> columnNames = joinTable.getColumnNames();
            final List<String> fullyQualifiedColumnNames = new LinkedList<String>();
            boolean found = false;
            for ( String columnName : columnNames ) {
                 if ( columnName.equals(localColumnName) )
                     found = true;
                 fullyQualifiedColumnNames.add(localBindingName + '.' + columnName);
            }
            if ( !found )
                throw new UserException.BadArgumentValue("-J", "The -J arg: \"" + joinArg + "\" specifies an unknown column name: \"" + localColumnName + "\". It's not one of the column names in the header " + columnNames + " of the file: " + filename);

            allFullyQualifiedColumnNames.addAll(fullyQualifiedColumnNames);
        }

        //parse the SELECT_COLUMNS arg and validate the column names
        List<String> parsedSelectColumns = new LinkedList<String>();
        for ( String token : SELECT_COLUMNS )
            parsedSelectColumns.addAll(Arrays.asList(token.split(",")));
        SELECT_COLUMNS = parsedSelectColumns.toArray(SELECT_COLUMNS);

        for ( String columnName : SELECT_COLUMNS ) {
            if ( !allFullyQualifiedColumnNames.contains(columnName) )
                throw new UserException.BadArgumentValue("-s", "The column name '" + columnName + "' provided to -s doesn't match any of the column names in any of the -B files. Here is the list of available column names: " + allFullyQualifiedColumnNames);
        }

        //instantiate the VariantAnnotatorEngine
        ArrayList<String> annotationsToUse = new ArrayList<String>();
        annotationsToUse.add("GenomicAnnotation");
        engine = new VariantAnnotatorEngine(getToolkit(), new ArrayList<String>(), annotationsToUse);
        engine.setOneToMany(ONE_TO_MANY);
        engine.setRequestedColumns(SELECT_COLUMNS);
        engine.setJoinTables(joinTables);

        // set up the header fields
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(GATKVCFUtils.getHeaderFields(getToolkit(), Arrays.asList("variant")));
        hInfo.addAll(engine.getVCFAnnotationDescriptions());

        Set<String> rodName = new HashSet<String>();
        rodName.add("variant");
        Set<String> samples = SampleUtils.getUniqueSamplesFromRods(getToolkit(), rodName);
        VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter.writeHeader(vcfHeader);
    }

    /**
     * Initialize the number of loci processed to zero.
     *
     * @return 0
     */
    public Integer reduceInit() { return 0; }

    /**
     * We want reads that span deletions
     *
     * @return true
     */
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    /**
     * For each site of interest, annotate based on the requested annotation types
     *
     * @param tracker  the meta-data tracker
     * @param ref      the reference base
     * @param context  the context for the given locus
     * @return 1 if the locus was successfully processed, 0 if otherwise
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        Set<VariantContext> results = new LinkedHashSet<VariantContext>();
        for (VariantContext vc : tracker.getValues(VariantContext.class, "variant", context.getLocation())) {
            if ( (vc.isFiltered() && IGNORE_FILTERED_SITES) ||
                    (vc.isVariant() && !vc.isBiallelic()) ) {
                results.add(vc);
            } else {
                Map<String, AlignmentContext> stratifiedContexts = AlignmentContextUtils.splitContextBySampleName(context);
                if ( stratifiedContexts != null )
                    results.addAll(engine.annotateContext(tracker, ref, stratifiedContexts, vc));
                else
                    results.add(vc);
            }
        }

        for ( VariantContext vc : results )
            vcfWriter.add(vc ,ref.getBase());

        return 1;
    }

    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public Integer treeReduce(Integer lhs, Integer rhs) {
        return lhs + rhs;
    }

    public void onTraversalDone(Integer sum) {

        //out.printf("Generated %d annotated VCF records.\n", totalOutputVCFRecords);
        Map<String, Integer> inputTableHitCounter = engine.getInputTableHitCounter();
        for ( Entry<String, Integer> e : inputTableHitCounter.entrySet() ) {
            final String bindingName = e.getKey();
            final int counter = e.getValue();
            //final float percent = 100 * counter /(float) totalOutputVCFRecords;
            //out.printf(" %-6.1f%%   (%d) annotated with %s.\n", percent, counter, bindingName );
            System.out.printf(" %d annotated with %s.\n", counter, bindingName );
        }
    }
}

