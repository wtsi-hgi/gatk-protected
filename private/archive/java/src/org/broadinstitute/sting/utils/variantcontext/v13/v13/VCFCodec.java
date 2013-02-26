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

package org.broadinstitute.variant.variantcontext.v13;

import org.broad.tribble.TribbleException;
import org.broad.tribble.readers.LineReader;
import org.broad.tribble.util.ParsingUtils;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * A feature codec for the VCF 4 specification
 *
 * <p>
 * VCF is a text file format (most likely stored in a compressed manner). It contains meta-information lines, a
 * header line, and then data lines each containing information about a position in the genome.
 * </p>
 * <p>One of the main uses of next-generation sequencing is to discover variation amongst large populations
 * of related samples. Recently the format for storing next-generation read alignments has been
 * standardised by the SAM/BAM file format specification. This has significantly improved the
 * interoperability of next-generation tools for alignment, visualisation, and variant calling.
 * We propose the Variant Call Format (VCF) as a standarised format for storing the most prevalent
 * types of sequence variation, including SNPs, indels and larger structural variants, together
 * with rich annotations. VCF is usually stored in a compressed manner and can be indexed for
 * fast data retrieval of variants from a range of positions on the reference genome.
 * The format was developed for the 1000 Genomes Project, and has also been adopted by other projects
 * such as UK10K, dbSNP, or the NHLBI Exome Project. VCFtools is a software suite that implements
 * various utilities for processing VCF files, including validation, merging and comparing,
 * and also provides a general Perl and Python API.
 * The VCF specification and VCFtools are available from http://vcftools.sourceforge.net.</p>
 *
 * <p>
 * See also: @see <a href="http://vcftools.sourceforge.net/specs.html">VCF specification</a><br>
 * See also: @see <a href="http://www.ncbi.nlm.nih.gov/pubmed/21653522">VCF spec. publication</a>
 * </p>
 *
 * <h2>File format example</h2>
 * <pre>
 *     ##fileformat=VCFv4.0
 *     #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA12878
 *     chr1    109     .       A       T       0       PASS  AC=1    GT:AD:DP:GL:GQ  0/1:610,327:308:-316.30,-95.47,-803.03:99
 *     chr1    147     .       C       A       0       PASS  AC=1    GT:AD:DP:GL:GQ  0/1:294,49:118:-57.87,-34.96,-338.46:99
 * </pre>
 *
 * @author Mark DePristo
 * @since 2010
 */
public class VCFCodec extends AbstractVCFCodec {
    // Our aim is to read in the records and convert to VariantContext as quickly as possible, relying on VariantContext to do the validation of any contradictory (or malformed) record parameters.

    public final static String VCF4_MAGIC_HEADER = "##fileformat=VCFv4";

    /**
     * @param reader the line reader to take header lines from
     * @return the number of header lines
     */
    public Object readHeader(LineReader reader) {
        List<String> headerStrings = new ArrayList<String>();

        String line;
        try {
            boolean foundHeaderVersion = false;
            while ((line = reader.readLine()) != null) {
                lineNo++;
                if (line.startsWith(VCFHeader.METADATA_INDICATOR)) {
                    String[] lineFields = line.substring(2).split("=");
                    if (lineFields.length == 2 && VCFHeaderVersion.isFormatString(lineFields[0]) ) {
                        if ( !VCFHeaderVersion.isVersionString(lineFields[1]) )
                            throw new TribbleException.InvalidHeader(lineFields[1] + " is not a supported version");
                        foundHeaderVersion = true;
                        version = VCFHeaderVersion.toHeaderVersion(lineFields[1]);
                        if ( version == VCFHeaderVersion.VCF3_3 || version == VCFHeaderVersion.VCF3_2 )
                            throw new TribbleException.InvalidHeader("This codec is strictly for VCFv4; please use the VCF3 codec for " + lineFields[1]);
                        if ( version != VCFHeaderVersion.VCF4_0 && version != VCFHeaderVersion.VCF4_1 )
                            throw new TribbleException.InvalidHeader("This codec is strictly for VCFv4 and does not support " + lineFields[1]);
                    }
                    headerStrings.add(line);
                }
                else if (line.startsWith(VCFHeader.HEADER_INDICATOR)) {
                    if (!foundHeaderVersion) {
                        throw new TribbleException.InvalidHeader("We never saw a header line specifying VCF version");
                    }
                    return createHeader(headerStrings, line);
                }
                else {
                    throw new TribbleException.InvalidHeader("We never saw the required CHROM header line (starting with one #) for the input VCF file");
                }

            }
        } catch (IOException e) {
            throw new RuntimeException("IO Exception ", e);
        }
        throw new TribbleException.InvalidHeader("We never saw the required CHROM header line (starting with one #) for the input VCF file");
    }


    /**
     * parse the filter string, first checking to see if we already have parsed it in a previous attempt
     *
     * @param filterString the string to parse
     * @return a set of the filters applied or null if filters were not applied to the record (e.g. as per the missing value in a VCF)
     */
    protected Set<String> parseFilters(String filterString) {
        return parseFilters(filterHash, lineNo, filterString);
    }

    public static Set<String> parseFilters(final Map<String, LinkedHashSet<String>> cache, final int lineNo, final String filterString) {
        // null for unfiltered
        if ( filterString.equals(VCFConstants.UNFILTERED) )
            return null;

        if ( filterString.equals(VCFConstants.PASSES_FILTERS_v4) )
            return Collections.emptySet();
        if ( filterString.equals(VCFConstants.PASSES_FILTERS_v3) )
            generateException(VCFConstants.PASSES_FILTERS_v3 + " is an invalid filter name in vcf4", lineNo);
        if ( filterString.length() == 0 )
            generateException("The VCF specification requires a valid filter status: filter was " + filterString, lineNo);

        // do we have the filter string cached?
        if ( cache != null && cache.containsKey(filterString) )
            return Collections.unmodifiableSet(cache.get(filterString));

        // empty set for passes filters
        LinkedHashSet<String> fFields = new LinkedHashSet<String>();
        // otherwise we have to parse and cache the value
        if ( filterString.indexOf(VCFConstants.FILTER_CODE_SEPARATOR) == -1 )
            fFields.add(filterString);
        else
            fFields.addAll(Arrays.asList(filterString.split(VCFConstants.FILTER_CODE_SEPARATOR)));

        fFields = fFields;
        if ( cache != null ) cache.put(filterString, fFields);

        return Collections.unmodifiableSet(fFields);
    }


    /**
     * create a genotype map
     * @param str the string
     * @param alleles the list of alleles
     * @return a mapping of sample name to genotype object
     */
    public Map<String, Genotype> createGenotypeMap(String str, List<Allele> alleles, String chr, int pos) {
        if (genotypeParts == null)
            genotypeParts = new String[header.getColumnCount() - NUM_STANDARD_FIELDS];

        int nParts = ParsingUtils.split(str, genotypeParts, VCFConstants.FIELD_SEPARATOR_CHAR);

        Map<String, Genotype> genotypes = new LinkedHashMap<String, Genotype>(nParts);

        // get the format keys
        int nGTKeys = ParsingUtils.split(genotypeParts[0], genotypeKeyArray, VCFConstants.GENOTYPE_FIELD_SEPARATOR_CHAR);

        // cycle through the sample names
        Iterator<String> sampleNameIterator = header.getGenotypeSamples().iterator();

        // clear out our allele mapping
        alleleMap.clear();

        // cycle through the genotype strings
        for (int genotypeOffset = 1; genotypeOffset < nParts; genotypeOffset++) {
            int GTValueSplitSize = ParsingUtils.split(genotypeParts[genotypeOffset], GTValueArray, VCFConstants.GENOTYPE_FIELD_SEPARATOR_CHAR);

            double GTQual = VariantContext.NO_NEG_LOG_10PERROR;
            Set<String> genotypeFilters = null;
            Map<String, String> gtAttributes = null;
            String sampleName = sampleNameIterator.next();

            // check to see if the value list is longer than the key list, which is a problem
            if (nGTKeys < GTValueSplitSize)
                generateException("There are too many keys for the sample " + sampleName + ", keys = " + parts[8] + ", values = " + parts[genotypeOffset]);

            int genotypeAlleleLocation = -1;
            if (nGTKeys >= 1) {
                gtAttributes = new HashMap<String, String>(nGTKeys - 1);

                for (int i = 0; i < nGTKeys; i++) {
                    final String gtKey = new String(genotypeKeyArray[i]);
                    boolean missing = i >= GTValueSplitSize;

                    // todo -- all of these on the fly parsing of the missing value should be static constants
                    if (gtKey.equals(VCFConstants.GENOTYPE_KEY)) {
                        genotypeAlleleLocation = i;
                    } else if (gtKey.equals(VCFConstants.GENOTYPE_QUALITY_KEY)) {
                        GTQual = missing ? parseQual(VCFConstants.MISSING_VALUE_v4) : parseQual(GTValueArray[i]);
                    } else if (gtKey.equals(VCFConstants.GENOTYPE_FILTER_KEY)) {
                        genotypeFilters = missing ? parseFilters(VCFConstants.MISSING_VALUE_v4) : parseFilters(getCachedString(GTValueArray[i]));
                    } else if ( missing ) {
                        gtAttributes.put(gtKey, VCFConstants.MISSING_VALUE_v4);
                    } else {
                        gtAttributes.put(gtKey, new String(GTValueArray[i]));
                    }
                }
            }

            // check to make sure we found a genotype field if we are a VCF4.0 file
            if ( version == VCFHeaderVersion.VCF4_0 && genotypeAlleleLocation == -1 )
                generateException("Unable to find the GT field for the record; the GT field is required in VCF4.0");
            if ( genotypeAlleleLocation > 0 )
                generateException("Saw GT field at position " + genotypeAlleleLocation + ", but it must be at the first position for genotypes when present");

            List<Allele> GTalleles = (genotypeAlleleLocation == -1 ? null : parseGenotypeAlleles(GTValueArray[genotypeAlleleLocation], alleles, alleleMap));
            boolean phased = genotypeAlleleLocation != -1 && GTValueArray[genotypeAlleleLocation].indexOf(VCFConstants.PHASED) != -1;

            // add it to the list
            try {
                genotypes.put(sampleName,
                        new Genotype(sampleName,
                                GTalleles,
                                GTQual,
                                genotypeFilters,
                                gtAttributes,
                                phased));
            } catch (TribbleException e) {
                throw new TribbleException.InternalCodecException(e.getMessage() + ", at position " + chr+":"+pos);
            }
        }

        return genotypes;
    }

    @Override
    public boolean canDecode(final File potentialInput) {
        return canDecodeFile(potentialInput, VCF4_MAGIC_HEADER);
    }
}
