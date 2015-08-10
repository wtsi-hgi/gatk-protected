/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 (“BROAD”) and the LICENSEE and is effective at the date the downloading is completed (“EFFECTIVE DATE”).
* 
* WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
* WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
* NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
* 
* 1. DEFINITIONS
* 1.1 PROGRAM shall mean copyright in the object code and source code known as GATK3 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/gatk on the EFFECTIVE DATE.
* 
* 2. LICENSE
* 2.1 Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
* The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only. For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
* 2.2 No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD. LICENSEE shall ensure that all of its users agree to the terms of this Agreement. LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
* 2.3 License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
* 
* 3. PHONE-HOME FEATURE
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system (“PHONE-HOME”) which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE’S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2015 Broad Institute, Inc.
* Notice of attribution: The GATK3 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
* LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
* 
* 5. INDEMNIFICATION
* LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
* 
* 6. NO REPRESENTATIONS OR WARRANTIES
* THE PROGRAM IS DELIVERED AS IS. BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
* IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
* 
* 7. ASSIGNMENT
* This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
* 
* 8. MISCELLANEOUS
* 8.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
* 8.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
* 8.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
* 8.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested. All notices under this Agreement shall be deemed effective upon receipt.
* 8.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
* 8.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
* 8.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.gatk.engine.alignment;

import org.broadinstitute.gatk.engine.alignment.bwa.BWAConfiguration;
import org.broadinstitute.gatk.engine.alignment.bwa.BWTFiles;
import org.broadinstitute.gatk.engine.alignment.bwa.c.BWACAligner;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.ReadWalker;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.util.Iterator;

/**
 * Validates consistency of the aligner interface
 *
 * <p>Validates consistency of the aligner interface by taking reads already aligned by BWA in a BAM file, stripping them
 * of their alignment data, realigning them, and making sure one of the best resulting realignments matches the original
 * alignment from the input file.</p>
 *
 * <h3>Caveat</h3>
 * <p>This tool requires that BWA be available on the java path.</p>
 *
 * @author mhanna
 * @version 0.1
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
public class CheckAlignment extends ReadWalker<Integer,Integer> {
    /**
     * The supporting BWT index generated using BWT. 
     */
    @Argument(fullName="BWTPrefix",shortName="BWT",doc="Index files generated by bwa index -d bwtsw",required=false)
    private String prefix = null;

    /**
     * The instance used to generate alignments.
     */
    private BWACAligner aligner = null;

    /**
     * Create an aligner object.  The aligner object will load and hold the BWT until close() is called.
     */
    @Override
    public void initialize() {
        if(prefix == null)
            prefix = getToolkit().getArguments().referenceFile.getAbsolutePath();
        BWTFiles bwtFiles = new BWTFiles(prefix);
        BWAConfiguration configuration = new BWAConfiguration();
        aligner = new BWACAligner(bwtFiles,configuration);
    }

    /**
     * Aligns a read to the given reference.
     *
     * @param ref Reference over the read.  Read will most likely be unmapped, so ref will be null.
     * @param read Read to align.
     * @return Number of reads aligned by this map (aka 1).
     */
    @Override
    public Integer map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker metaDataTracker) {
        //logger.info(String.format("examining read %s", read.getReadName()));

        byte[] bases = read.getReadBases();
        if(read.getReadNegativeStrandFlag()) bases = BaseUtils.simpleReverseComplement(bases);

        boolean matches = true;
        Iterable<Alignment[]> alignments = aligner.getAllAlignments(bases);
        Iterator<Alignment[]> alignmentIterator = alignments.iterator();

        if(!alignmentIterator.hasNext()) {
            matches = read.getReadUnmappedFlag();
        }
        else {
            Alignment[] alignmentsOfBestQuality = alignmentIterator.next();
            for(Alignment alignment: alignmentsOfBestQuality) {
                matches = (alignment.getContigIndex() == read.getReferenceIndex());
                matches &= (alignment.getAlignmentStart() == read.getAlignmentStart());
                matches &= (alignment.isNegativeStrand() == read.getReadNegativeStrandFlag());
                matches &= (alignment.getCigar().equals(read.getCigar()));
                matches &= (alignment.getMappingQuality() == read.getMappingQuality());
                if(matches) break;
            }
        }

        if(!matches) {
            logger.error("Found mismatch!");
            logger.error(String.format("Read %s:",read.getReadName()));
            logger.error(String.format("    Contig index: %d",read.getReferenceIndex()));
            logger.error(String.format("    Alignment start: %d", read.getAlignmentStart()));
            logger.error(String.format("    Negative strand: %b", read.getReadNegativeStrandFlag()));
            logger.error(String.format("    Cigar: %s%n", read.getCigarString()));
            logger.error(String.format("    Mapping quality: %s%n", read.getMappingQuality()));
            for(Alignment[] alignmentsByScore: alignments) {
                for(int i = 0; i < alignmentsByScore.length; i++) {
                    logger.error(String.format("Alignment %d:",i));
                    logger.error(String.format("    Contig index: %d",alignmentsByScore[i].getContigIndex()));
                    logger.error(String.format("    Alignment start: %d", alignmentsByScore[i].getAlignmentStart()));
                    logger.error(String.format("    Negative strand: %b", alignmentsByScore[i].isNegativeStrand()));
                    logger.error(String.format("    Cigar: %s", alignmentsByScore[i].getCigarString()));
                    logger.error(String.format("    Mapping quality: %s%n", alignmentsByScore[i].getMappingQuality()));
                }
            }
            throw new ReviewedGATKException(String.format("Read %s mismatches!", read.getReadName()));
        }

        return 1;
    }

    /**
     * Initial value for reduce.  In this case, validated reads will be counted.
     * @return 0, indicating no reads yet validated.
     */
    @Override
    public Integer reduceInit() { return 0; }

    /**
     * Calculates the number of reads processed.
     * @param value Number of reads processed by this map.
     * @param sum Number of reads processed before this map.
     * @return Number of reads processed up to and including this map.
     */
    @Override
    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    /**
     * Cleanup.
     * @param result Number of reads processed.
     */
    @Override
    public void onTraversalDone(Integer result) {
        aligner.close();
        super.onTraversalDone(result);
    }

}
