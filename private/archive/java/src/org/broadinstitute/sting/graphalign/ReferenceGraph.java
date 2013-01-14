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

package org.broadinstitute.sting.gatk.walkers.graphalign;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleDirectedGraph;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.*;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.util.IntervalTree;
import net.sf.samtools.util.StringUtil;

import java.util.*;
import java.io.Serializable;
import java.io.IOException;
import java.io.ObjectInputStream;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Oct 22, 2009
 * Time: 2:15:54 PM
 * To change this template use File | Settings | File Templates.
 */
class ReferenceGraph extends SimpleDirectedGraph<Fragment, DefaultEdge> implements Serializable {
    final private static boolean USE_IT = true;
    final private static boolean THROW_ERRORS_ON_BAD_INPUTS = false;

    private boolean DEBUG = false;
    int nSkippedIndels = 0;
    int nSkippedBecauseOfContinguousVariation = 0;
    int nBadPolymorphisms = 0;
    int nMultiStateAlleles = 0;

    GenomeLoc initialLoc = null;

    private transient IntervalTree<Fragment> loc2Fragment = new IntervalTree<Fragment>();

    public ReferenceGraph(boolean printDebuggingInfo) {
        super(DefaultEdge.class);
        DEBUG = printDebuggingInfo;
    }

    public ReferenceGraph() {
        this(false);
    }

    public void setDebugPrinting(boolean enable) {
        this.DEBUG = enable;
    }

    public void bindRefenceSequence(ReferenceSequence seq) {
        GenomeLoc refSeqLoc = GenomeLocParser.createGenomeLoc(seq.getContigIndex(), 1, seq.length());
        String refString = StringUtil.bytesToString(seq.getBases()).toUpperCase(); 
        Fragment frag = new Fragment(refSeqLoc, 1, StringUtil.stringToBytes(refString));
        addFragment(frag);
        initialLoc = refSeqLoc;
    }

    private void addFragmentToIntervalTree(Fragment frag) {
        loc2Fragment.put((int)frag.getLocation().getStart(), (int)frag.getLocation().getStop(), frag);
    }

    private void addFragment(Fragment frag) {
        addFragmentToIntervalTree(frag);
        addVertex(frag);
    }

    private void removeFragment(Fragment frag) {
        loc2Fragment.remove((int)frag.getLocation().getStart(), (int)frag.getLocation().getStop());
        removeVertex(frag);
    }

    public void validateGraph() {
        for ( Fragment v : this.vertexSet() ) {
            if ( this.inDegreeOf(v) == 0 && v.getLocation().getStart() != initialLoc.getStart() ) {
                throw new StingException(String.format("Fragment %s has no incoming edges but isn't at the start of the contig %s", v, initialLoc));
            }
            if ( this.outDegreeOf(v) == 0 && v.getLocation().getStop() != initialLoc.getStop() ) {
                throw new StingException(String.format("Fragment %s has no outgoing edges but isn't at the end of the contig %s", v, initialLoc));
            }
        }
        
        //System.out.printf("Passed validation: %s%n", this.toBriefString());
    }

    private void rebuildIntervalTree() {
        if ( DEBUG ) System.out.printf("rebuilding IntervalTree()%n");
        for ( Fragment v : this.vertexSet() ) {
            if ( DEBUG ) System.out.printf("  adding interval tree: %s%n", v);
            addFragmentToIntervalTree(v);
        }
    }

    private boolean allelesAreInExcisedFragment(Fragment cut, List<String> alleles) {
        boolean foundRef = false;
        for ( String allele : alleles ) {
            if ( allele.equals(cut.getBases()) ) 
                foundRef = true;
        }

        if ( ! foundRef && THROW_ERRORS_ON_BAD_INPUTS )
            throw new StingException(String.format("Polymorphic alleles %s do not contain the reference sequence %s", alleles, cut.getBases()));
        
        return foundRef;
    }

    public void addVariation(VariantContext variant, GenomeLoc loc, List<String> alleles) {
        if ( DEBUG ) System.out.printf("addVariation(%s, %s)%n", loc, alleles);
        //validateGraph();

        if ( variant.isSNP() ) {
            Fragment frag = getContainingFragment(loc);

            if ( frag == null ) {
                nMultiStateAlleles++;
                return;
            }

            if ( ! allelesAreInExcisedFragment(subFragment(frag, loc, 1), alleles)) {
                nBadPolymorphisms++;
                return;
            }

            List<Fragment> split = exciseFragment(frag, loc);
            if ( split != null ) {
                Fragment left = split.get(0);
                Fragment cut = split.get(1);
                Fragment right = split.get(2);

                if ( DEBUG ) System.out.printf("  cutFrag(%s, %s)%n", loc, cut);

                for ( String allele : alleles ) {
                    byte[] bases = StringUtil.stringToBytes(allele);
                    double freq = 1.0 / alleles.size();
                    Fragment alleleFrag = new Fragment(loc, freq, 0, bases.length, bases);
                    if ( DEBUG ) System.out.printf("  Creating allele fragment %s%n", alleleFrag);
                    addFragment(alleleFrag);
                    if ( left != null ) addEdge(left, alleleFrag);
                    if ( right != null ) addEdge(alleleFrag, right);
                }
            } else {
                nSkippedBecauseOfContinguousVariation++;
            }
        } else {
            nSkippedIndels++;
        }
    }


    private Fragment subFragment(Fragment frag, GenomeLoc loc, double freq ) {
        return new Fragment(loc, 1, frag.getUnderlyingBases());
    }

    public List<Fragment> exciseFragment(Fragment frag, GenomeLoc loc) {
        if ( DEBUG ) System.out.printf("  exciseFragment(%s, %s)%n", frag, loc);
        GenomeLoc fragLoc = frag.getLocation();

        Fragment cut = subFragment(frag, loc, 1);

        Set<DefaultEdge> inToFrag = incomingEdgesOf(frag);
        Set<DefaultEdge> outOfFrag = outgoingEdgesOf(frag);

        Fragment left = null;
        if ( fragLoc.getStart() == loc.getStart() ) {
            if ( ! inToFrag.isEmpty() ) {
                if ( THROW_ERRORS_ON_BAD_INPUTS )
                    throw new StingException(String.format("Attempting to create a variation at the start of a fragment %s %s", frag, loc));
                return null;
            }
        } else {
            GenomeLoc leftLoc = GenomeLocParser.createGenomeLoc(fragLoc.getContigIndex(), fragLoc.getStart(), loc.getStart()-1);
            left = new Fragment(leftLoc, 1, frag.getUnderlyingBases());
            addFragment(left);

            for ( DefaultEdge e : inToFrag ) {
                addEdge(getEdgeSource(e), left);
            }

            removeAllEdges(inToFrag);
        }

        Fragment right = null;
        if ( fragLoc.getStop() == loc.getStop() ) {
            if ( ! outOfFrag.isEmpty() ) {
                throw new StingException(String.format("Attempting to create a variation at the end of a fragment %s %s", frag, loc));
            }
        } else {
            GenomeLoc rightLoc = GenomeLocParser.createGenomeLoc(fragLoc.getContigIndex(), loc.getStop()+1, fragLoc.getStop());
            right = new Fragment(rightLoc, 1, frag.getUnderlyingBases());
            addFragment(right);

            for ( DefaultEdge e : outOfFrag ) {
                addEdge(right, getEdgeTarget(e));
            }

            removeAllEdges(outOfFrag);
        }

        if ( DEBUG ) System.out.printf("    removing %s%n", frag);
        removeFragment(frag);
        if ( DEBUG ) System.out.printf("    returning left=%s right=%s%n", left, right);
        return Arrays.asList(left, cut, right);
    }

    public Fragment getContainingFragment(GenomeLoc loc) {
        Fragment frag = USE_IT ? getContainingFragmentIT(loc) : getContainingFragmentG(loc);


        if ( frag == null ) {
            if ( THROW_ERRORS_ON_BAD_INPUTS )
                throw new StingException("No spanning fragment was found for " + loc);
            else
                return null;
        }
        else if ( frag.getLocation().getStart() > loc.getStart() || frag.getLocation().getStop() < loc.getStop() )
            throw new StingException("BUG: bad spanning fragment found for " + loc + " was " + frag.getLocation() );
        else
            return frag;
    }

    public Fragment getContainingFragmentG(GenomeLoc loc) {
        for ( Fragment v : this.vertexSet() ) {
            if ( v.getLocation().containsP(loc) ) {
                return v;
            }
        }

        return null;
    }

    public Fragment getContainingFragmentIT(GenomeLoc loc) {
        IntervalTree.Node<Fragment> node = loc2Fragment.minOverlapper((int)loc.getStart(), (int)loc.getStop());
        if ( node == null )
            return null;
        else
            return node.getValue();
    }

    public Collection<Fragment> getStartingFragment(GenomeLoc loc) {
        Collection<Fragment> frags = USE_IT ? getStartingFragmentIT(loc) : getStartingFragmentG(loc);
        //Collection<Fragment> frags = getStartingFragmentTest(loc);

        if ( frags == null || frags.size() == 0 )
            throw new StingException("No fragment contains location start of " + loc);
        if ( frags.size() == 1 && MathUtils.compareDoubles(frags.iterator().next().getFrequency(), 1.0) != 0 ) {
            Fragment bad = frags.iterator().next();
            throw new StingException(String.format("Only one fragment was found but it's frequency < 1 %s with %e", bad, bad.getFrequency()));
        }
        else
            return frags;
    }

    public Collection<Fragment> getStartingFragmentTest(GenomeLoc loc) {
        Collection<Fragment> fragsFromIT = getStartingFragmentIT(loc);
        Collection<Fragment> fragsFromG = getStartingFragmentG(loc);

        if ( fragsFromIT.size() != fragsFromG.size() ) {
            throw new StingException(String.format("Fragment sizes differ %d from IntervalTree, %d from graph", fragsFromIT.size(), fragsFromG.size()));
        }

        return USE_IT && false ? fragsFromIT : fragsFromG;
    }

    public Collection<Fragment> getStartingFragmentIT(GenomeLoc loc) {
        Collection<Fragment> frags = new HashSet<Fragment>();

        Iterator<IntervalTree.Node<Fragment>> it = loc2Fragment.overlappers((int)loc.getStart(), (int)loc.getStart());
        IntervalTree.Node<Fragment> node = null;
        while ( it.hasNext() ) {
            node = it.next();
            frags.add(node.getValue());
        }

        // todo -- painful bug work around -- should be removed
        if ( frags.size() == 1 && MathUtils.compareDoubles(node.getValue().getFrequency(), 1.0) != 0 ) {
            System.out.printf(">>> Using IT workaround at %s <<<%n", loc);
            return getStartingFragmentG(loc);
        }
        
        return frags;
//        IntervalTree.Node<Fragment> node = loc2Fragment.minOverlapper((int)loc.getStart(), (int)loc.getStart());
//        if ( node == null )
//            return null;
//        else
//            return node.getValue();
    }

    public Collection<Fragment> getStartingFragmentG(GenomeLoc loc) {
        Collection<Fragment> frags = new HashSet<Fragment>();
        for ( Fragment v : this.vertexSet() ) {
            //if ( v.getLocation().getStart() < loc.getStop() )
            //    System.out.printf("Checking %s vs. %s%n", loc, v.getLocation());
            if ( v.getLocation().containsStartPosition(loc.getStart()) ) {
            //    System.out.printf(" Adding %s%n", v.getLocation());
                frags.add(v);
            }
        }

        return frags;
    }

    public Set<Fragment> outgoingFragments( Fragment frag ) {
        Set<Fragment> outgoingFrags = new HashSet<Fragment>();
        
        for ( DefaultEdge e : outgoingEdgesOf(frag) ) {
            outgoingFrags.add(getEdgeTarget(e));
        }

        if ( outgoingFrags.size() == 0 && frag.getLocation().getStop() != initialLoc.getStop() ) {

        }

        return outgoingFrags;
    }

    public String toString() {
        StringBuilder s = new StringBuilder();

        for ( Fragment v : this.vertexSet() ) {
            s.append(String.format("Fragment: %s%n", v.toString()));
            for ( DefaultEdge e : this.incomingEdgesOf(v) ) {
                s.append(String.format("  [IN FROM] %s%n", this.getEdgeSource(e)));
            }
            for ( DefaultEdge e : this.outgoingEdgesOf(v) ) {
                s.append(String.format("  [OUT TO ] %s%n", this.getEdgeTarget(e)));
            }
        }

        return s.toString();
    }

    public String toBriefString() {
        return String.format("GraphRef: %d fragments, %d edges, skipped %d contingous variants, %d indels, %d polymorphisms w/o ref allele, %d multi-state",
                this.vertexSet().size(), this.edgeSet().size(), nSkippedBecauseOfContinguousVariation, nSkippedIndels, nBadPolymorphisms, nMultiStateAlleles);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // serialization
    //
    // --------------------------------------------------------------------------------------------------------------
    private void readObject(ObjectInputStream stream) throws IOException, ClassNotFoundException {
        //always perform the default de-serialization first
        stream.defaultReadObject();
        loc2Fragment = new IntervalTree<Fragment>();
        rebuildIntervalTree();
    }
}

class Fragment implements Serializable {
    GenomeLoc loc = null;   // Index position of this fragment into the reference
    int offset, stop;
    double freq = -1;
    byte[] bases = null;

    public Fragment( GenomeLoc loc, double freq, int offset, int stop, byte[] bases ) {
        this.loc = loc;
        this.freq = freq;
        this.bases = bases;
        this.offset = offset;
        this.stop = stop;
    }

    public Fragment( GenomeLoc loc, double freq, byte[] bases ) {
        this(loc, freq, (int)loc.getStart()-1, (int)loc.getStop(), bases);
    }

    public String toString() {
        //return String.format("%s:%.2f:%s", loc.toString(), getFrequency(), getBases());
        return String.format("%s:%.2f", loc.toString(), getFrequency());
    }

    public GenomeLoc getLocation() {
        return loc;
    }

    public double getFrequency() {
        return freq;
    }

    public int getUnderlyingOffset() { return offset; }
    public int getStop() { return stop; }
    public int getLength() { return getStop() - getUnderlyingOffset(); }

    public byte[] getUnderlyingBases() {
        return bases;
    }

    /**
     * how many bases over in the fragment are we over in this fragment?
     *
     * @param loc
     * @return
     */
    public int getFragOffsetFrom(GenomeLoc loc) {
        // todo -- ignores contigs -- can we fix this?
        if ( getLocation().getStart() > loc.getStart() )
            throw new StingException("BUG: Request for offset from " + loc + " in frag at " + getLocation() + " but this is beyond the location of the fragment");
        return (int)(loc.getStart() - getLocation().getStart());
    }

    public int getBaseLengthFrom( int fragOffset, int maxLength ) {
        int fragRemaining = getLength() - fragOffset;

        if ( fragRemaining < 0 )
            throw new StingException("BUG: Request for length from offset " + fragOffset + " but this is longer than the fragment itself");
        
        return Math.min(fragRemaining, maxLength);
    }

    public byte getBase(int fragOffset) {
        return bases[getUnderlyingOffset() + fragOffset];
    }

    public String getBases() {
        return StringUtil.bytesToString(getUnderlyingBases(), getUnderlyingOffset(), getLength());
    }
}
