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

package org.broadinstitute.sting.gatk.features.maf;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
import org.broadinstitute.variant.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.variant.variantcontext.*;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Jan 24, 2011
 * Time: 12:01:12 PM
 * To change this template use File | Settings | File Templates.
 */
public class MafFeature implements Feature {
    private String contig;                      // our contig location
    private int start;                          // our starting location, zero based
    private int stop;                           // our stopping location

    private String refAllele = ".";             // the reference allele
    private String[] observedTumAlleles = null;          // The sequences of the observed alleles in tumor
    private String[] observedNormAlleles = null;          // The sequences of the observed alleles in normal
    private String tumorSampleId = null;
    private String normalSampleId = null;
    private String hugoSymbol = null;
    private Classification classification = null;

    public enum Type {
        UNKNOWN,SNP,MNP,INS,DEL
    };

    public enum Classification {
        Unclassified, Intergenic,Intron,Noncoding_transcript,UTR3,UTR5,Flank5,Silent,Missense, Nonsense, Splice_site, miRNA,
        Frameshift, Inframe, Stop_deletion, Promoter,De_novo_start, De_novo_start_out_of_frame
    }

    private Type type = Type.UNKNOWN;

    /**
     * create the dbSNP feature, given the following information:
     *
     * @param contig the contig rsID
     * @param start  the start position, one based
     * @param stop   the stop position, one based
     */
    MafFeature(String contig,
                 int start,
                 int stop) {
        this.contig = contig;
        this.start = start;
        this.stop = stop;
    }

    public void setVariantType(String t) {
        type=Type.valueOf(t);
    }

    public void setObservedTumor(String[] obs) {
        observedTumAlleles = obs;
    }

    public void setObservedTumor(String allele1, String allele2) {
        observedTumAlleles = new String[2];
        observedTumAlleles[0] = allele1;
        observedTumAlleles[1] = allele2;
    }

    public void setRefAllele(String ref) {
        this.refAllele = ref;
    }

    public void setTumorSample(String sampleId) {
        this.tumorSampleId = sampleId;
    }

    public void setNormalSample(String sampleId) {
        this.normalSampleId = sampleId;
    }

    public String getRefBases() {
        return refAllele;
    }

    public String getHugoGeneSymbol() {
        return hugoSymbol;
    }

    public void setHugoGeneSymbol(String genename) {
        int pos = genename.indexOf('|');
        if ( pos < 0 ) {
            hugoSymbol = genename;
        } else {
            hugoSymbol = genename.substring(0,pos);
        }
    }

    /**
     * Returns list of alleles (represented as strings) observed in Tumor. Returned alleles
     * could be redundant (e.g. if we have homozygous non-ref at ploidy 2+).
     * @return
     */
    public List<String> getObservedTumorAlleleList() {
        return Arrays.asList(observedTumAlleles);
    }

    /**
     * Returns list of alleles (represented as strings) observed in Tumor. Returned alleles
     * could be redundant (e.g. if we have homozygous non-ref at ploidy 2+).
     * @return
     */
    public List<String> getObservedNormalAlleleList() {
        if ( observedNormAlleles == null ) {
            // if we got no ref allele observations recorded in the maf, we assume its ref/ref (somatic event)
            List<String> l = new ArrayList<String>(2);
            l.add(refAllele);
            l.add(refAllele);
            return l;
        }
        else return Arrays.asList(observedTumAlleles);
    }

    /** Returns a (non-redundant) list of all distinct alleles
     * observed at the site, plus a reference allele (whether it
     * was actually observed or not). The reference allele is always returned as first
     * element of the list.
     * @return
     */
    public List<String> getAllAlleleList() {
        List<String> l = new ArrayList<String>();
        l.add(refAllele);
        for ( String a : observedTumAlleles ) {
            if ( l.contains(a) ) continue;
            l.add(a);
        }
        if ( observedNormAlleles != null ) {
            for ( String a : observedNormAlleles ) {
                if ( l.contains(a) ) continue;      // already have this allele
                l.add(a);
            }
        }
        return l;
    }

    /** Returns a (non-redundant) list of all distinct non-reference alleles
     * observed at the site
     * @return
     */
    public List<String> getAllNonRefAlleleList() {
        List<String> l = new ArrayList<String>();
        for ( String a : observedTumAlleles ) {
            if ( l.contains(a) ) continue;      // already have this allele
            if ( a.equals(refAllele)) continue; // allele is ref, we do not need it
            l.add(a); 
        }
        if ( observedNormAlleles != null ) {
            for ( String a : observedNormAlleles ) {
                if ( l.contains(a) ) continue;      // already have this allele
                if ( a.equals(refAllele)) continue; // allele is ref, we do not need it
                l.add(a);
            }
        }
        return l;
    }

    public String getTumorSampleId() { return tumorSampleId; }
    public String getNormalSampleId() { return normalSampleId; }

    public boolean isRefAllele(String a) { return refAllele.equals(a); }

    public Type getType() { return type; }

    public int lengthOnRef() {
        switch ( type ) {
            case SNP:
            case MNP:
            case DEL:
                return refAllele.length();
            case INS:
                return 0;
            default:
                throw new StingException("Unrecognized event type in Maf record: "+type);
        }
    }

    public boolean isSomatic() {
        if ( observedTumAlleles[0].equals(refAllele) && observedTumAlleles[1].equals(refAllele) ) return false; // tumor is ref
        // we get here only if tumor is non-ref
        if ( observedNormAlleles == null ) return true; // norm alleles are omitted from maf only if they are all ref
        if ( observedNormAlleles[0].equals(refAllele) && observedNormAlleles[1].equals(refAllele) ) return true;
        return false;
    }

    public void setVariantClassification(String s) {
        if ( s.equals("IGR") ) { classification = Classification.Intergenic ; return; }
        if ( s.equals("Intron") ) { classification = Classification.Intron ; return; }
        if ( s.equals("3'UTR") || s.equals("3'-UTR")) { classification = Classification.UTR3 ; return; }
        if ( s.equals("5'UTR") || s.equals("5'-UTR")) { classification = Classification.UTR5 ; return; }
        if ( s.equals("5'-Flank") ) { classification = Classification.Flank5 ; return; }
        if ( s.equals("Silent") || s.equals("Synonymous")) { classification = Classification.Silent ; return; }
        if ( s.equals("Non-coding_Transcript")) { classification = Classification.Noncoding_transcript; return; }
        if ( s.equals("Missense") || s.equals("Missense_Mutation") ) { classification = Classification.Missense ; return; }
        if ( s.equals("Nonsense_Mutation") || s.equals("Nonsense") ) { classification = Classification.Nonsense ; return; }
        if ( s.equals("Splice_Site") ) { classification = Classification.Splice_site ; return; }
        if ( s.equals("miRNA") ) { classification = Classification.miRNA ; return; }
        if ( s.equals("Frame_Shift_Ins") ) { classification = Classification.Frameshift ; return; }
        if ( s.equals("Frame_Shift_Del") ) { classification = Classification.Frameshift ; return; }
        if ( s.equals("In_Frame_Ins") ) { classification = Classification.Inframe ; return; }
        if ( s.equals("In_Frame_Del") ) { classification = Classification.Inframe ; return; }
        if ( s.equals("Stop_Codon_Del") ) { classification = Classification.Stop_deletion ; return; }
        if ( s.equals("Splice_Site_Del") ) { classification = Classification.Splice_site ; return; }
        if ( s.equals("Splice_Site_Ins") ) { classification = Classification.Splice_site ; return; }
        if ( s.equals("Splice_Site_SNP") ) { classification = Classification.Splice_site ; return; }
        if ( s.equals("Promoter") ) { classification = Classification.Promoter ; return; }
        if ( s.equals("De_novo_Start") ) { classification = Classification.De_novo_start ; return; }
        if ( s.equals("De_novo_Start_OutOfFrame") ) { classification = Classification.De_novo_start_out_of_frame ; return; }
        if ( s.equals("TX-REF-MISMATCH") ) { classification = Classification.Unclassified ; return; }
        throw new UserException.MalformedFile("Unknown variant classification: " + s);
    }

    public Classification getVariantClassification() {
        return classification;
    }

   /*
     * the required getting and setter methods
     */

    public String getChr() {
        return contig;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return stop;
    }
    
}

class MafAdaptor implements VariantContextAdaptors.VCAdaptor {
    /**
     * Converts Maf features to VariantContext.
     * @return MafFeature.
     */
    @Override
    public Class<? extends Feature> getAdaptableFeatureType() { return MafFeature.class; }

    /**
     * convert to a Variant Context, given:
     * @param name the name of the ROD
     * @param input the Rod object, in this case a MafFeature
     * @return a VariantContext object
     */
//        VariantContext convert(String name, Object input) {
//            return convert(name, input, null);
//        }

    /**
     * convert to a Variant Context, given:
     * @param name  the name of the ROD
     * @param input the Rod object, in this case a MafFeature
     * @param ref   the reference context
     * @return a VariantContext object
     */
    @Override
    public VariantContext convert(String name, Object input, ReferenceContext ref) {

        if ( ref == null )
            throw new UnsupportedOperationException("Conversion from MAF to VariantContext requires a reference context, null received");

        MafFeature maf = (MafFeature)input;
        if ( ! Allele.acceptableAlleleBases(maf.getRefBases()) )
            return null;

        List<Allele> alleles = new ArrayList<Allele>();

        Allele refAllele = Allele.create(maf.getRefBases(), true);
        // add the reference allele:
        alleles.add(refAllele);

        // add all of the alt alleles
        for ( String alt : maf.getAllNonRefAlleleList() ) {
            if ( ! Allele.acceptableAlleleBases(alt) ) {
                //System.out.printf("Excluding dbsnp record %s%n", dbsnp);
                return null;
            }
            alleles.add(Allele.create(alt, false));
        }

        // make a mapping from sample to genotype

        String normalSample = maf.getNormalSampleId();
        String tumorSample = maf.getTumorSampleId();

//                String[] genotypeStrings = hapmap.getGenotypes();

        GenotypesContext genotypes = GenotypesContext.create(2);

        addGenotype(genotypes, normalSample, maf.getObservedNormalAlleleList(),maf.getRefBases());
        addGenotype(genotypes,tumorSample,maf.getObservedTumorAlleleList(),maf.getRefBases());


        HashMap<String, Object> attrs = new HashMap<String, Object>(10);
        // fill attributes:
        if ( maf.getHugoGeneSymbol() != null && ! maf.getHugoGeneSymbol().equals("Unknown"))
            attrs.put("Gene",maf.getHugoGeneSymbol());

        if ( maf.isSomatic() ) {
            attrs.put(VCFConstants.SOMATIC_KEY,true);
            attrs.put("SS","Somatic");
        } else {
            attrs.put("SS","Germline");
        }

        if ( maf.getVariantClassification() != null ) {
            switch(maf.getVariantClassification()) {
                case Intergenic: attrs.put("VC","Genomic"); break;
                case Intron: attrs.put("VC","Intron"); break;
                case Noncoding_transcript: attrs.put("VC","Noncoding_transcript"); break;
                case UTR3: attrs.put("VC","3'UTR"); break;
                case UTR5: attrs.put("VC","5'UTR"); break;
                case Flank5: attrs.put("VC","5'flank"); break;
                case Promoter: attrs.put("VC","5'flank"); break;
                case De_novo_start: attrs.put("VC","De_novo_start"); break;
                case De_novo_start_out_of_frame: attrs.put("VC","De_novo_start_out_of_frame"); break;
                case Silent: attrs.put("VC","Silent"); break;
                case Missense: attrs.put("VC","Missense"); break;
                case Nonsense: attrs.put("VC","Nonsense"); break;
                case Splice_site: attrs.put("VC","Splice_site"); break;
                case miRNA: attrs.put("VC","miRNA"); break;
                case Frameshift: attrs.put("VC","Frameshift"); break;
                case Inframe: attrs.put("VC","Inframe"); break;
                case Stop_deletion: attrs.put("VC","Stop_codon_deletion");
                case Unclassified: attrs.put("VC","Unclassified");
                default:
            }
        }

        attrs.put("VT",maf.getType());

        int end = maf.getEnd();
        VariantContext vc = new VariantContextBuilder(name, maf.getChr(), maf.getStart(), end, alleles)
                .genotypes(genotypes).attributes(attrs).make();
        return vc;
    }

    private void addGenotype(GenotypesContext dest, String sampleId, List<String> alleles, String refAllele) {
        List<Allele> myAlleles = new ArrayList<Allele>(2);

        boolean success = true;

        for ( String a : alleles ) {
            if ( a.isEmpty() || a.contains("N") || a.contains(".")) return; // bad allele found
            myAlleles.add(Allele.create(a,refAllele.equals(a)));
        }
        dest.add(GenotypeBuilder.create(sampleId,myAlleles));
    }

}