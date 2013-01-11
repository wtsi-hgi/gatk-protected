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

package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broad.tribble.readers.AsciiLineReader;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.MappingQualityZeroFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.io.*;
import java.util.*;

/**
 * Walks along all variant ROD loci and verifies the phasing from the reads for user-defined pairs of sites.
 */
@Allows(value = {DataSource.READS, DataSource.REFERENCE})
@Requires(value = {DataSource.READS, DataSource.REFERENCE})
@By(DataSource.READS)

@ReadFilters({MappingQualityZeroFilter.class})
// Filter out all reads with zero mapping quality

public class ReadBasedPhasingValidationWalker extends RodWalker<Integer, Integer> {
    @Input(fullName="variant", shortName = "V", doc="Validate variants from this VCF file", required=true)
    public RodBinding<VariantContext> variants;

    @Argument(fullName = "sitePairsFile", shortName = "sitePairsFile", doc = "File of pairs of variants for which phasing in ROD should be assessed using input reads", required = true)
    protected File sitePairsFile = null;

    @Output
    protected PrintStream out;

    private Set<SitePair> sitePairs = null;
    private String sampleName = null;

    SiteGenotypeAndReads prevSiteAndReads = null;

    private final static int NUM_IN_PAIR = 2; // trivial

    // enable deletions in the pileup
    public boolean includeReadsWithDeletionAtLoci() {
        return true;
    }

    public void initialize() {
        sitePairs = new TreeSet<SitePair>();
        GenomeLocParser locParser = getToolkit().getGenomeLocParser();

        InputStream sitePairsStream = null;
        try {
            sitePairsStream = new FileInputStream(sitePairsFile);
        } catch (FileNotFoundException fnfe) {
            fnfe.printStackTrace();
            throw new UserException("Problem opening file: " + sitePairsFile);
        }

        AsciiLineReader sitePairsReader = new AsciiLineReader(sitePairsStream);
        while (true) {
            String line = null;
            try {
                line = sitePairsReader.readLine();
            } catch (IOException ioe) {
                ioe.printStackTrace();
                throw new UserException("Problem reading file: " + sitePairsFile);
            }
            if (line == null)
                break; // reached end of file

            String[] twoSites = line.split("\t");
            if (twoSites.length != 2)
                throw new UserException("Must have PAIRS of sites in line " + line + " of " + sitePairsFile);

            SitePair sp = new SitePair(locParser.parseGenomeLoc(twoSites[0]), locParser.parseGenomeLoc(twoSites[1]));
            sitePairs.add(sp);
        }
    }

    public boolean generateExtendedEvents() {
        return false;
    }

    public Integer reduceInit() {
        return 0;
    }

    /**
     * @param tracker the meta-data tracker
     * @param ref     the reference base
     * @param context the context for the given locus
     * @return statistics of and list of all phased VariantContexts and their base pileup that have gone out of cacheWindow range.
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return null;

        boolean relevantSitePair = false;
        SitePair sp = null;
        if (prevSiteAndReads != null) {
            // all vc's below start at ref.getLocus() [due to requireStartHere = true]:
            sp = new SitePair(prevSiteAndReads.site, ref.getLocus());
            relevantSitePair = sitePairs.contains(sp);
        }

        if (context == null || !context.hasBasePileup())
            return null;
        ReadBackedPileup pileup = context.getBasePileup();
        String nextName = null;

        Collection<String> sampNames = pileup.getSamples();
        if (sampNames.size() != 1)
            throw new UserException("Reads must be for exactly one sample [not multi-sample]");
        nextName = sampNames.iterator().next();
        if (nextName == null)
            throw new UserException("Reads must be for exactly one sample");

        if (sampleName == null)
            sampleName = nextName;
        else if (!nextName.equals(sampleName))
            throw new UserException("Reads must have a single consistent sample name");

        pileup = pileup.getPileupForSample(sampleName);

        ReadBasesAtPosition readBases = new ReadBasesAtPosition();
        for (PileupElement p : pileup)
            readBases.putReadBase(p);

        ReadCounter rdCounts = null;
        if (relevantSitePair) { // otherwise, processed the reads for their possible use in the future:
            PhasingReadList buildReads = new PhasingReadList(NUM_IN_PAIR);
            buildReads.updateBases(0, prevSiteAndReads.readBases);
            buildReads.updateBases(1, readBases);

            List<PhasingRead> reads = new LinkedList<PhasingRead>();
            for (Map.Entry<String, PhasingRead> readEntry : buildReads.entrySet()) {
                PhasingRead rd = readEntry.getValue();
                if (rd.getNonNullIndices().length == NUM_IN_PAIR) { // only want reads with BOTH bases called [possibly as deleted ("D")]
                    reads.add(rd);
                    logger.debug("Read name: " + readEntry.getKey() + "\trd: " + rd);
                }
            }

            // Count the occurence of each "haplotype":
            rdCounts = new ReadCounter();
            for (PhasingRead rd : reads)
                rdCounts.incrementCount(rd);
        }


        // Now, read the ROD and note the genotypes and their phase to be validated:
        Set<Haplotype> calledHaplotypes = null;
        List<Haplotype> allPossibleHaplotypes = null;

        for (VariantContext vc : Arrays.asList(tracker.getFirstValue(variants, context.getLocation()))) {
            if (vc.isFiltered() || !vc.isSNP())
                continue;

            if (vc.getNSamples() != 1)
                throw new UserException("ROD file must have exactly one sample [not multi-sample]");
            nextName = vc.getSampleNames().iterator().next();
            if (sampleName == null)
                sampleName = nextName;
            else if (!nextName.equals(sampleName))
                throw new UserException("ROD must have a single consistent sample name");

            Genotype gt = vc.getGenotype(sampleName);

            if (relevantSitePair) {
                Genotype prevGt = prevSiteAndReads.gt;
                List<Allele> prevAlleles = prevGt.getAlleles();
                List<Allele> curAlleles = gt.getAlleles();

                calledHaplotypes = new TreeSet<Haplotype>(); // implemented Haplotype.compareTo()
                if (gt.isPhased()) {
                    if (gt.getPloidy() != prevGt.getPloidy())
                        throw new UserException("Invalid ROD file: cannot be phased AND have different ploidys!");

                    // Consider only the haplotypes called to be phased
                    Iterator<Allele> curAllIt = curAlleles.iterator();
                    for (Allele prevAll : prevAlleles) {
                        Allele curAll = curAllIt.next();
                        calledHaplotypes.add(successiveAllelesToHaplotype(prevAll, curAll));
                    }
                }

                // Consider EVERY combination of alleles as haplotypes [IF PHASED, this will give the contingency table in the CORRECT order]:
                allPossibleHaplotypes = new LinkedList<Haplotype>();
                for (Allele prevAll : prevAlleles) {
                    for (Allele curAll : curAlleles) {
                        allPossibleHaplotypes.add(successiveAllelesToHaplotype(prevAll, curAll));
                    }
                }
            }

            prevSiteAndReads = new SiteGenotypeAndReads(ref.getLocus(), gt, readBases);
        }

        int processedPairs = 0;
        if (relevantSitePair) {
            Map<Haplotype, Integer> haplotypeCounts = new TreeMap<Haplotype, Integer>(); // implemented Haplotype.compareTo()

            processedPairs = 1;
            int totalCount = rdCounts.totalCount();
            System.out.println("\nPair: " + sp + " [# reads = " + totalCount + "]");

            int matchCount = 0;
            for (Map.Entry<PhasingRead, Integer> rdEntry : rdCounts.entrySet()) {
                PhasingRead read = rdEntry.getKey();
                int count = rdEntry.getValue();

                Haplotype readsHaplotype = new Haplotype(read);
                haplotypeCounts.put(readsHaplotype, count);

                boolean readMatchesCalledHaplotype = calledHaplotypes != null && calledHaplotypes.contains(readsHaplotype);
                if (readMatchesCalledHaplotype)
                    matchCount += count;

                System.out.println("read" + ": " + read + (readMatchesCalledHaplotype ? "*" : "") + "\tcount: " + count);
            }

            double percentMatchingReads = 100 * (matchCount / (double) totalCount);
            System.out.println("% MATCHING reads: " + percentMatchingReads + " [of " + totalCount + " TOTAL reads]");

            out.print(sp);
            if (allPossibleHaplotypes != null) {
                for (Haplotype hap : allPossibleHaplotypes) {
                    Integer count = haplotypeCounts.get(hap);
                    if (count == null) // haplotype may not have been observed in ANY reads
                        count = 0;

                    out.print("\t" + count);
                }
            }
            out.println();
        }

        return processedPairs;
    }

    private Haplotype successiveAllelesToHaplotype(Allele prevAll, Allele curAll) {
        byte prevBase = SNPallelePair.getSingleBase(prevAll);
        byte curBase = SNPallelePair.getSingleBase(curAll);

        byte[] hapBases = new byte[NUM_IN_PAIR];
        hapBases[0] = prevBase;
        hapBases[1] = curBase;
        return new Haplotype(hapBases);
    }

    public Integer reduce(Integer addIn, Integer runningCount) {
        if (addIn == null)
            addIn = 0;

        return runningCount + addIn;
    }

    /**
     * @param result the number of reads and VariantContexts seen.
     */
    public void onTraversalDone(Integer result) {
        System.out.println("Validated " + result + " pairs of sites.");
    }
}

class SitePair implements Comparable<SitePair> {
    public GenomeLoc site1;
    public GenomeLoc site2;

    public SitePair(GenomeLoc site1, GenomeLoc site2) {
        if (site1.size() > 1 || site2.size() > 1)
            throw new UserException("Must give pairs of SINGLE-LOCUS record start sites");

        this.site1 = site1;
        this.site2 = site2;
    }

    public String toString() {
        return site1.toString() + "\t" + site2.toString();
    }

    public int compareTo(SitePair other) {
        int comp1 = site1.compareTo(other.site1);
        if (comp1 != 0)
            return comp1;

        return site2.compareTo(other.site2);
    }
}

class SiteGenotypeAndReads {
    public GenomeLoc site;
    public Genotype gt;
    public ReadBasesAtPosition readBases;

    public SiteGenotypeAndReads(GenomeLoc site, Genotype gt, ReadBasesAtPosition readBases) {
        this.site = site;
        this.gt = gt;
        this.readBases = readBases;
    }
}

class PhasingReadList {
    private Map<String, PhasingRead> readsAtSites = null;
    private int numSites;

    public PhasingReadList(int numSites) {
        this.readsAtSites = new HashMap<String, PhasingRead>();
        this.numSites = numSites;
    }

    public void updateBases(int index, ReadBasesAtPosition readBases) {
        if (readBases == null)
            return;

        for (ReadBase rb : readBases) {
            String readName = rb.readName;

            PhasingRead rd = readsAtSites.get(readName);
            if (rd == null) {
                rd = new PhasingRead(numSites, rb.mappingQual);
                readsAtSites.put(readName, rd);
            }

            // Arbitrarily updates to the last base observed for this sample and read (rb.base):
            rd.updateBaseAndQuality(index, rb.base, rb.baseQual);
        }
    }

    public Set<Map.Entry<String, PhasingRead>> entrySet() {
        return readsAtSites.entrySet();
    }

    public int size() {
        return readsAtSites.size();
    }
}

class ReadCounter {
    private Map<PhasingRead, Integer> counts;
    private int totalCount;

    public ReadCounter() {
        this.counts = new TreeMap<PhasingRead, Integer>(); // implemented PhasingRead.compareTo()
    }

    public void incrementCount(PhasingRead rd) {
        Integer cnt = counts.get(rd);
        if (cnt == null)
            cnt = 0;

        counts.put(rd, cnt + 1);
        totalCount++;
    }

    public Set<Map.Entry<PhasingRead, Integer>> entrySet() {
        return counts.entrySet();
    }

    public int totalCount() {
        return totalCount;
    }
}