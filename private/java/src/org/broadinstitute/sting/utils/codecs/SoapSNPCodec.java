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

package org.broadinstitute.sting.utils.codecs;

import org.broad.tribble.*;
import org.broad.tribble.exception.CodecLineParsingException;
import org.broadinstitute.variant.variantcontext.*;

import java.util.*;

/**
 * A codec for parsing soapsnp files
 *
 * <p>
 * A simple text file format with the following whitespace separated fields of 17 columns:
 * <ol>
 *     <li>Chromosome ID</li>
 *     <li>Coordinate on chromosome, start from 1</li>
 *     <li>Reference genotype</li>
 *     <li>Consensus genotype</li>
 *     <li>Quality score of consensus genotype</li>
 *     <li>Best base</li>
 *     <li>Average quality score of best base</li>
 *     <li>Count of uniquely mapped best base</li>
 *     <li>Count of all mapped best base</li>
 *     <li>Second best bases</li>
 *     <li>Average quality score of second best base</li>
 *     <li>Count of uniquely mapped second best base</li>
 *     <li>Count of all mapped second best base</li>
 *     <li>Sequencing depth of the site</li>
 *     <li>Rank sum test p_value</li>
 *     <li>Average copy number of nearby region</li>
 *     <li>Whether the site is a dbSNP</li>
 * </ol>
 * Note this codec is for internal use only, and is not supported outside of GSA.
 * </p>
 *
 * <p>
 * See also: @see <a href="http://soap.genomics.org.cn/soapsnp.html#usage2">SOAPSNP usage page</a><br>
 * </p>

 * </p>
 *
 * <h2>File format example</h2>
 * <pre>
 *     chr1    205     A       C       2       C       18      2       2       A       0       0       0       2       2       1.00000 1.00000 0
 *     chr1    492     C       Y       19      T       34      2       2       C       34      1       1       3       3       0.666667        1.00000 1
 *     chr1    1540    G       C       3       C       34      2       2       G       0       0       1       3       3       1.00000 1.33333 0
 *     chr1    1555    A       C       3       C       33      2       2       A       0       0       0       2       2       1.00000 1.00000 0
 *     chr1    4770    A       G       14      G       33      6       8       A       0       0       3       11      11      1.00000 1.54545 0
 *     chr1    4793    A       G       17      G       33      7       8       A       0       0       3       11      11      1.00000 1.45455 0
 *     chr1    126137  C       S       0       G       34      1       1       C       0       0       1       2       2       0.00000 1.50000 0
 *     chr1    218136  G       R       36      G       34      32      378     A       34      7       58      436     436     0.507135        1.91055 0
 *     chr1    218178  G       R       54      G       33      44      655     A       34      9       185     841     841     0.504643        1.93698 0
 *     chr1    218326  G       S       20      G       33      100     1665    C       33      7       60      1727    1727    0.462359        1.93804 0
 * </pre>
 *
 * @author Mark DePristo
 * @since 2010
 */
public class SoapSNPCodec extends AsciiFeatureCodec<VariantContext> implements NameAwareCodec {
    private String[] parts;

    // we store a name to give to each of the variant contexts we emit
    private String name = "Unknown";

    public SoapSNPCodec() {
        super(VariantContext.class);
    }

    /**
     * Decode a line as a Feature.
     *
     * @param line
     *
     * @return Return the Feature encoded by the line,  or null if the line does not represent a feature (e.g. is
     *         a comment)
     */
    public VariantContext decode(String line) {
        try {
            // parse into lines
            parts = line.trim().split("\\s+");

            // check that we got the correct number of tokens in the split
            if (parts.length != 18)
                throw new CodecLineParsingException("Invalid SoapSNP row found -- incorrect element count.  Expected 18, got " + parts.length + " line = " + line);

            String contig = parts[0];
            long start = Long.valueOf(parts[1]);
            AlleleAndGenotype allelesAndGenotype = parseAlleles(parts[2], parts[3], line);

            double log10PError = Integer.valueOf(parts[4]) / -10.0;

            Map<String, Object> attributes = new HashMap<String, Object>();
            attributes.put("BestBaseQ", parts[6]);
            attributes.put("SecondBestBaseQ", parts[10]);
            attributes.put("RankSumP", parts[15]);
            // add info to keys

            //System.out.printf("Alleles  = " + allelesAndGenotype.alleles);
            //System.out.printf("genotype = " + allelesAndGenotype.genotype);
            
            VariantContext vc = new VariantContextBuilder(name, contig, start, start, allelesAndGenotype.alleles).genotypes(allelesAndGenotype.genotype).log10PError(log10PError).passFilters().attributes(attributes).make();

            //System.out.printf("line  = %s%n", line);
            //System.out.printf("vc    = %s%n", vc);

            return vc;
        } catch (CodecLineParsingException e) {
            throw new TribbleException("Unable to parse line " + line,e);
        } catch (NumberFormatException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            throw new TribbleException("Unable to parse line " + line,e);
        }
    }

    private static class AlleleAndGenotype {
        Collection<Allele> alleles;
        Collection<Genotype> genotype;

        public AlleleAndGenotype(Collection<Allele> alleles, Genotype genotype) {
            this.alleles = alleles;
            this.genotype = new HashSet<Genotype>();
            this.genotype.add(genotype);
        }
    }

    private AlleleAndGenotype parseAlleles(String ref, String consensusGenotype, String line) {
        /* A	Adenine
    C	Cytosine
    G	Guanine
    T (or U)	Thymine (or Uracil)
    R	A or G
    Y	C or T
    S	G or C
    W	A or T
    K	G or T
    M	A or C
    B	C or G or T
    D	A or G or T
    H	A or C or T
    V	A or C or G
    N	any base
    . or -	gap
    */
        if ( ref.equals(consensusGenotype) )
            throw new TribbleException.InternalCodecException("Ref base and consensus genotype are the same " + ref);

        Allele refAllele = Allele.create(ref, true);
        List<Allele> genotypeAlleles = null;

        char base = consensusGenotype.charAt(0);

        switch ( base ) {
            case 'A': case 'C': case 'G': case 'T':
                Allele a = Allele.create(consensusGenotype);
                genotypeAlleles = Arrays.asList(a, a);
                break;
            case 'R': case 'Y': case 'S': case 'W': case 'K': case 'M':
                genotypeAlleles = determineAlt(refAllele, ref.charAt(0), base);
                break;
            default:
                throw new TribbleException("Unexpected consensus genotype " + consensusGenotype + " at line = " + line);
        }


        Collection<Allele> alleles = new HashSet<Allele>(genotypeAlleles);
        alleles.add(refAllele);
        Genotype genotype = GenotypeBuilder.create("unknown", genotypeAlleles); // todo -- probably should include genotype quality

        return new AlleleAndGenotype( alleles, genotype );
    }

    private static final Map<Character, String> IUPAC_SNPS = new HashMap<Character, String>();
    static {
        IUPAC_SNPS.put('R', "AG");
        IUPAC_SNPS.put('Y', "CT");
        IUPAC_SNPS.put('S', "GC");
        IUPAC_SNPS.put('W', "AT");
        IUPAC_SNPS.put('K', "GT");
        IUPAC_SNPS.put('M', "AC");
    }

    private List<Allele> determineAlt(Allele ref, char refbase, char alt) {
        String alts = IUPAC_SNPS.get(alt);
        if ( alts == null )
            throw new IllegalStateException("BUG: unexpected consensus genotype " + alt);
            
        Allele a1 = alts.charAt(0) == refbase ? ref : Allele.create((byte)alts.charAt(0));
        Allele a2 = alts.charAt(1) == refbase ? ref : Allele.create((byte)alts.charAt(1));

        //if ( a1 != ref && a2 != ref )
        //    throw new IllegalStateException("BUG: unexpected consensus genotype " + alt + " does not contain the reference base " + ref);

        return Arrays.asList(a1, a2);
    }

    /**
     * get the name of this codec
     * @return our set name
     */
    public String getName() {
        return name;
    }

    /**
     * set the name of this codec
     * @param name new name
     */
    public void setName(String name) {
        this.name = name;
    }

    public static void main(String[] args) {
        System.out.printf("Testing " + args[0]);
    }
}