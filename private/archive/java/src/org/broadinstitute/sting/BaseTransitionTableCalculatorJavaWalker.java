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

package org.broadinstitute.sting.walkers;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Oct 12, 2009
 * Time: 2:43:06 PM
 * To change this template use File | Settings | File Templates.
 */
@By(DataSource.REFERENCE)
@Reference(window = @Window(start = -3, stop = 3))
public class BaseTransitionTableCalculatorJavaWalker extends LocusWalker<Set<BaseTransitionTable>, Set<BaseTransitionTable>> implements TreeReducible<Set<BaseTransitionTable>> {
    @Output
    PrintStream out;

    @Argument(fullName = "usePreviousBases", doc = "Use previous bases of the reference as part of the calculation, uses the specified number, defaults to 0", required = false)
    int nPreviousBases = 0;
    @Argument(fullName = "useSecondaryBase", doc = "Use the secondary base of a read as part of the calculation", required = false)
    boolean useSecondaryBase = false;
    @Argument(fullName = "confidentRefThreshold", doc = "Set the lod score that defines confidence in ref, defaults to 4", required = false)
    int confidentRefThreshold = 5;
    @Argument(fullName = "maxNumMismatches", doc = "Set the maximum number of mismatches at a locus before choosing not to use it in calculation. Defaults to 1.", required = false)
    int maxNumMismatches = 1;
    @Argument(fullName = "minMappingQuality", doc = "Set the alignment quality below which to ignore reads; defaults to 30", required = false)
    int minMappingQuality = 30;
    @Argument(fullName = "minQualityScore", doc = "Set the base quality score below which to ignore bases in the pileup, defaults to 20", required = false)
    int minQualityScore = 20;
    @Argument(fullName = "usePileupMismatches", doc = "Use the number of mismatches in the pileup as a condition for the table", required = false)
    boolean usePileupMismatches = false;
    @Argument(fullName = "usePreviousReadBases", doc = "Use previous bases of the read as part of the calculation. Will ignore reads if there aren't this many previous bases. Uses the specified number. Defaults to 0", required = false)
    int nPreviousReadBases = 0;
    @Argument(fullName = "useReadGroup", doc = "Use the group number of the read as a condition of the table.", required = false)
    boolean useReadGroup = false;
    @Argument(fullName = "outputFile", shortName = "of", doc = "Output to this file rather than standard out. Must be used with -nt.", required = false)
    String outFilePath = null;
    @Argument(fullName = "forcePreviousReadBasesToMatchRef", doc = "Forces previous read bases to match the reference", required = false)
    boolean readBasesMustMatchRef = false;

    private UnifiedGenotyperEngine ug;
    // private ReferenceContextWindow refWindow;
    // private Set<BaseTransitionTable> conditionalTables;
    private List<Boolean> usePreviousBases;
    private List<GenomeLoc> previousBaseLoci;

    public void initialize() {
        if (nPreviousBases > 3 || (nPreviousReadBases > 3 && readBasesMustMatchRef)) {
            throw new UserException.CommandLineException("You have opted to use a number of previous bases in excess of 3. In order to do this you must change the reference window size in the walker itself.");
        }
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.baseModel = BaseMismatchModel.THREE_STATE;
        uac.ALL_BASES_MODE = true;
        ug = new UnifiedGenotyperEngine(getToolkit(), uac);
        // refWindow = new ReferenceContextWindow(nPreviousBases);
        usePreviousBases = new ArrayList<Boolean>();
        previousBaseLoci = new ArrayList<GenomeLoc>();

    }

    public Set<BaseTransitionTable> reduceInit() {
        return new TreeSet<BaseTransitionTable>();
    }

    public Set<BaseTransitionTable> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        ReadBackedPileup pileup = context.getBasePileup();
        Set<BaseTransitionTable> newCounts = null;
        //System.out.println(pileup.getBases());
        if (baseIsUsable(tracker, ref, pileup, context)) {
            //System.out.println("Pileup will be used");
            if (previousLociCanBeUsed(usePreviousBases, previousBaseLoci, context.getLocation())) {
                for (int r = 0; r < pileup.getReads().size(); r++) {
                    if (useRead(pileup.getReads().get(r), pileup.getOffsets().get(r), ref)) {
                        newCounts = updateTables(newCounts, pileup.getReads().get(r), pileup.getOffsets().get(r), ref, pileup);
                    }
                }
            }
            else {
                updatePreviousBases(usePreviousBases, true, previousBaseLoci, context.getLocation());
            }
        }
        else {
            updatePreviousBases(usePreviousBases, false, previousBaseLoci, context.getLocation());
        }

        return newCounts;
    }

    public Set<BaseTransitionTable> reduce(Set<BaseTransitionTable> map, Set<BaseTransitionTable> reduce) {
        if (map != null && !map.isEmpty()) {
            for (BaseTransitionTable t : map) {
                boolean add = true;
                for (BaseTransitionTable r : reduce) {
                    if (r.conditionsMatch(t)) {
                        r.incorporateTable(t);
                        add = false;
                        break;
                    }
                }
                if (add) {
                    reduce.add(t);
                }
            }
        }
        // System.out.println("Reduce: size of TransitionTable set is " + reduce.size() + " -- size of Map: " + (map != null ? map.size() : "null"));
        return reduce;
    }

    public Set<BaseTransitionTable> treeReduce(Set<BaseTransitionTable> reduce1, Set<BaseTransitionTable> reduce2) {
        // check to see if this is a truly tree-reducable calculation
        if (nPreviousBases >= 1) {
            String errMsg = "Parallelization cannot be used with UsePreviousBases due to the fact that internal walker data specifies whether a previous reference base is usable or not.";
            String errMsg2 = " This can cause cause concurrency issues and unpredictable behavior when used with parallelization. Either do not specify -nt, or try a the conjunction of ";
            String errMsg3 = "--usePreviousReadBases and --forcePreviousReadBasesToMatchRef.";
            throw new UserException.CommandLineException(errMsg + errMsg2 + errMsg3);
        }
        return reduce(reduce1, reduce2);
    }

    public void onTraversalDone(Set<BaseTransitionTable> conditionalTables) {
        PrintStream output;
        if (outFilePath == null) {
            output = out;
        }
        else {
            try {
                output = new PrintStream(outFilePath);
            } catch (FileNotFoundException e) {
                throw new UserException.CouldNotCreateOutputFile(new File(outFilePath), e);
            }
        }
        output.print(createHeaderFromConditions());
        for (BaseTransitionTable t : conditionalTables)
            t.print(output);
    }

    public void updatePreviousBases(List<Boolean> usage, boolean canUse, List<GenomeLoc> loci, GenomeLoc locus) {
        // early return
        if (nPreviousBases < 1) {
            return;
        }

        if (usage.size() <= nPreviousBases) {
            usage.add(canUse);
            loci.add(locus);
        }
        else {
            usage.remove(0);
            usage.add(canUse);
            loci.remove(0);
            loci.add(locus);
        }
    }

    public boolean previousLociCanBeUsed(List<Boolean> canUse, List<GenomeLoc> loci, GenomeLoc locus) {
        if (nPreviousBases < 1) {
            return true;
        }

        boolean use = true;
        for (boolean b : canUse) {
            use = use && b;
        }

        if (use) {
            use = use && (loci.get(0).distance(locus) == 1); // truly is PREVIOUS base
        }

        return use;
    }

    public Set<BaseTransitionTable> updateTables(Set<BaseTransitionTable> tables, SAMRecord read, int offset, ReferenceContext ref, ReadBackedPileup pileup) {
        List<Comparable> readConditions = buildConditions(read, offset, ref, pileup);
        // System.out.println("Updating table with pileup: "+pileup.getBases()+ ( read.getReadNegativeStrandFlag() ? "-" : "+" ) + "  Quality: "+read.getBaseQualities()[offset] + "  MapQ:  "+read.getMappingQuality());

        if (tables == null) {
            tables = new TreeSet<BaseTransitionTable>();
        }

        boolean createNewTable = true;

        for (BaseTransitionTable t : tables) {
            if (t.conditionsMatch(readConditions)) {
                updateTable(t, read, offset, ref);
                createNewTable = false;
                break;
            }
        }

        if (createNewTable) {
            BaseTransitionTable t = new BaseTransitionTable(readConditions);
            updateTable(t, read, offset, ref);
            tables.add(t);
        }

        return tables;
    }

    public void updateTable(BaseTransitionTable t, SAMRecord r, int o, ReferenceContext ref) {
        // System.out.println("Update Table");
        if (r.getReadNegativeStrandFlag()) {
            t.update((byte) BaseUtils.simpleComplement((char) r.getReadBases()[o]), (byte) BaseUtils.simpleComplement(ref.getBaseAsChar()));
        }
        else {
            t.update(r.getReadBases()[o], ref.getBase());
        }
    }

    public boolean useRead(SAMRecord read, int offset, ReferenceContext ref) {

        if (Character.toUpperCase(read.getReadBases()[offset]) == Character.toUpperCase(ref.getBase())) {
            return false;
        }
        else if (read.getMappingQuality() <= minMappingQuality) {
            return false;
        }
        else if (!BaseUtils.isRegularBase(read.getReadBases()[offset])) {
            return false;
        }
        else if (read.getBaseQualities()[offset] <= minQualityScore) {
            return false;
        }
        else if (useSecondaryBase && read.getAttribute("SQ") == null) {
            return false;
        }
        else if (nPreviousBases >= 1 && previousReadBasesMismatchRef(read, offset, ref)) {
            return false;
        }
        else if (nPreviousReadBases >= 1 && readLacksPreviousBases(read, offset, nPreviousReadBases)) {
            return false;
        }
        else if (nPreviousReadBases >= 1 && readBasesMustMatchRef && previousReadBasesMismatchRef(read, offset, ref)) {
            return false;
        }
        else {
            return true;
        }
    }

    public boolean previousReadBasesMismatchRef(SAMRecord read, int offset, ReferenceContext ref) {
        int c = read.getReadNegativeStrandFlag() ? 1 : -1;
        if (offset + nPreviousBases * c < 0) {
            return true;
        }
        else if (offset + nPreviousBases * c > read.getReadLength()) {
            return true;
        }

        for (int prevBase = 1; prevBase <= nPreviousBases; prevBase++) {
            if (Character.toUpperCase(read.getReadBases()[offset + prevBase * c]) != Character.toUpperCase(ref.getBases()[nPreviousBases + 1 + prevBase * c]) || !BaseUtils.isRegularBase(ref.getBases()[nPreviousBases + 1 + prevBase * c])) {
                return true;
            }
        }

        return false;
    }

    public boolean readLacksPreviousBases(SAMRecord read, int offset, int prevBases) {
        if (!read.getReadNegativeStrandFlag()) {
            return offset - prevBases < 0;
        }
        else {
            return offset + prevBases + 1 >= read.getReadLength();
        }
    }

    public List<Comparable> buildConditions(SAMRecord read, int offset, ReferenceContext ref, ReadBackedPileup pileup) {
        ArrayList<Comparable> conditions = new ArrayList<Comparable>();

        if (nPreviousBases > 0) {
            conditions.add(buildRefString(ref, nPreviousBases, !read.getReadNegativeStrandFlag()));

        }

        if (useSecondaryBase) {
            conditions.add(getSecondaryBase(read, offset));
        }

        if (nPreviousReadBases > 0) {
            conditions.add(buildReadString(read, offset, nPreviousReadBases));
        }

        if (usePileupMismatches) {
            conditions.add(countMismatches(ref.getBase(), pileup));
        }

        if (useReadGroup) {
            conditions.add(read.getReadGroup().getReadGroupId());
        }

        return conditions;
    }

    public String buildRefString(ReferenceContext ref, int bases, boolean forwardRead) {
        if (forwardRead) {
            return (new String(ref.getBases())).substring(0, nPreviousBases - 1);
        }
        else {
            return BaseUtils.simpleReverseComplement((new String(ref.getBases())).substring(nPreviousBases + 1));
        }
    }

    public String buildReadString(SAMRecord read, int offset, int nPreviousReadBases) {
        if (!read.getReadNegativeStrandFlag()) {
            return read.getReadString().substring(offset - nPreviousReadBases, offset);
        }
        else {
            return BaseUtils.simpleReverseComplement(read.getReadString().substring(offset + 1, offset + nPreviousReadBases + 1));
        }
    }

    public String createHeaderFromConditions() {
        String header = "Observed_base\tTrue_base";

        if (nPreviousBases > 0) {
            header = header + "\tPrevious_" + nPreviousBases + "_bases";
        }

        if (useSecondaryBase) {
            header = header + "\tSecondary_base";
        }

        if (nPreviousReadBases > 0) {
            header = header + "\tPrevious_" + nPreviousReadBases + "_read_bases";
        }

        if (usePileupMismatches) {
            header = header + "\tNumber_of_pileup_mismatches";
        }

        if (useReadGroup) {
            header = header + "\tRead_group";
        }

        return String.format("%s\t%s%n", header, "Counts");
    }

    public int countMismatches(byte ref, ReadBackedPileup p) {
        int refM = p.getBaseCounts()[BaseUtils.simpleBaseToBaseIndex(ref)];
        return p.size() - refM;
    }

    public char getSecondaryBase(SAMRecord read, int offset) {
        return BaseUtils.baseIndexToSimpleBaseAsChar(QualityUtils.compressedQualityToBaseIndex(((byte[]) read.getAttribute("SQ"))[offset]));
    }

    public boolean baseIsUsable(RefMetaDataTracker tracker, ReferenceContext ref, ReadBackedPileup pileup, AlignmentContext context) {
        return pileupContainsNoNs(pileup) && baseIsConfidentRef(tracker, ref, context) && pileupBelowMismatchThreshold(ref, pileup);
    }

    public boolean pileupBelowMismatchThreshold(ReferenceContext ref, ReadBackedPileup pileup) {
        return countMismatches(ref.getBase(), pileup) <= maxNumMismatches;
    }

    public boolean pileupContainsNoNs(ReadBackedPileup pileup) {
        for (byte c : pileup.getBases()) {
            if (c == 'N') {
                return false;
            }
        }

        return true;
    }

    public boolean baseIsConfidentRef(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (!BaseUtils.isRegularBase(ref.getBase()))
            return false;
        VariantCallContext calls = ug.calculateLikelihoodsAndGenotypes(tracker, ref, context);
        if (calls == null || calls.vc == null)
            return false;
        return (calls.vc.getNSamples() > 0 && calls.vc.getGenotype(0).isHomRef() && calls.vc.getGenotype(0).getNegLog10PError() > confidentRefThreshold);

    }

}

class BaseTransitionTable implements Comparable {

    /*
     * no direct manipulation of these objects ever
     */
    private int[][] table;
    private List<Comparable> conditions;

    public BaseTransitionTable(List<Comparable> conditions) {
        table = new int[BaseUtils.BASES.length][BaseUtils.BASES.length];
        for (int i = 0; i < BaseUtils.BASES.length; i++) {
            for (int j = 0; j < BaseUtils.BASES.length; j++) {
                table[i][j] = 0;
            }
        }

        this.conditions = conditions;
    }

    public boolean conditionsMatch(Object obj) {
        if (obj == null) {
            return false;
        }
        else if (obj instanceof BaseTransitionTable) {
            return ((BaseTransitionTable) obj).conditionsMatch(conditions);
        }
        else if (!(obj instanceof List)) {

            return false;
        }
        else if (this.numConditions() != ((List) obj).size()) {
            return false;
        }
        else {
            boolean eq = true;
            ListIterator thisIter = this.getConditionIterator();
            ListIterator thatIter = ((List) obj).listIterator();

            while (thisIter.hasNext()) {
                eq = eq && thisIter.next().equals(thatIter.next());
            }

            return eq;
        }
    }

    public int compareTo(Object obj) {
        if (!(obj instanceof BaseTransitionTable)) {
            return -1;
        }
        else {
            BaseTransitionTable t = (BaseTransitionTable) obj;
            if (this.conditionsMatch(t.conditions)) {
                return 0;
            }
            else {
                if (this.numConditions() == t.numConditions()) {
                    ListIterator<Comparable> thisIter = this.conditions.listIterator();
                    ListIterator<Comparable> thatIter = t.conditions.listIterator();
                    int g = 0;
                    do {
                        g = thisIter.next().compareTo(thatIter.next());
                    } while (g == 0);

                    return g;

                }
                else {
                    return (this.numConditions() > t.numConditions()) ? 1 : -1;
                }
            }
        }

    }

    public void print(PrintStream out) {
        StringBuilder s = new StringBuilder();
        for (byte observedBase : BaseUtils.BASES) {
            for (byte refBase : BaseUtils.BASES) {
                s.append(String.format("%s\t%s", (char) observedBase, (char) refBase));
                for (Comparable c : conditions) {
                    s.append(String.format("\t%s", c.toString()));
                }
                s.append(String.format("\t%d%n", table[BaseUtils.simpleBaseToBaseIndex(observedBase)][BaseUtils.simpleBaseToBaseIndex(refBase)]));
            }
        }

        out.print(s.toString());
    }

    public void update(byte observedBase, byte refBase) {
        //if ( observedBase == refBase ) {
        //    throw new StingException("BaseTransitionTable received equal observed and reference bases, which should not happen.");
        //}
        // System.out.println("Table updating: Observed Base: "+observedBase+" Ref base: "+refBase);
        table[BaseUtils.simpleBaseToBaseIndex(observedBase)][BaseUtils.simpleBaseToBaseIndex(refBase)]++;
    }

    public int numConditions() {
        return conditions.size();
    }

    private Comparable getCondition(int offset) {
        return conditions.get(offset);
    }

    private ListIterator getConditionIterator() {
        return conditions.listIterator();
    }

    public void incorporateTable(BaseTransitionTable t) {
        for (int i = 0; i < BaseUtils.BASES.length; i++) {
            for (int j = 0; j < BaseUtils.BASES.length; j++) {
                table[i][j] += t.observationsOf(i, j);
            }
        }
    }

    public int observationsOf(int observedBaseIndex, int referenceBaseIndex) {
        return table[observedBaseIndex][referenceBaseIndex];
    }

}