package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.codecs.table.TableFeature;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;

/**
 * A walker to convert table-formatted tracks into VCF files. Please see the wiki
 * for how to format your table file for this conversion.
 */
@Reference(window=@Window(start=-50,stop=50))
public class TableToVCF extends RodWalker<VariantContext,Integer> {
    @Input(shortName = "-t", fullName = "--table", doc = "The input table we will convert to VCF")
    RodBinding<TableFeature> table;

    //final private String CNV_HEADER = "HEADER,loc,size,type";
    /**
     * Todo -- when ##ALT= header lines can be written into VCF, enable CNV tables to be used
     */
    final private String VAR_HEADER = "HEADER,loc,size,ref,alt,type";

    private Set<Pair<GenomeLoc,String>> active;

    @Output
    VCFWriter vcfWriter = null;

    public void initialize() {
        active = new HashSet<Pair<GenomeLoc,String>>();
        vcfWriter.writeHeader(new VCFHeader(new HashSet<VCFHeaderLine>()));
    }

    public Integer reduceInit() { return 0; }

    public VariantContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) { return null; }

        VariantContext vcToPrint = null;
        for ( TableFeature tFeature : tracker.getValues(table) ) {
            if ( ! Utils.join(",",tFeature.getHeader()).equals(VAR_HEADER) ) {
                throw new UserException("Invalid Header Format");
            }
            Pair<GenomeLoc,String> rep = new Pair<GenomeLoc,String>(tFeature.getLocation(),String.format("%s.%s",table.getName(),Utils.join(",",tFeature.getAllValues())));
            if ( ! active.contains(rep) ) {
                active.add(rep);
                logger.debug(String.format("%s\t%s%n", rep.getFirst(), rep.getSecond()));
                if ( Utils.join(",",tFeature.getHeader()).equals(VAR_HEADER) ) {
                    vcToPrint = parseVar(tFeature.getAllValues(),ref,tFeature.getLocation());
                }
            }
        }

        HashSet<Pair<GenomeLoc,String>> toRem = new HashSet<Pair<GenomeLoc, String>>();
        for ( Pair<GenomeLoc,String> act : active ) {
            if ( act.getFirst().isBefore(ref.getLocus()) ) {
                toRem.add(act);
            }
        }

        active.removeAll(toRem);
        vcfWriter.add(vcToPrint);
        return vcToPrint;
    }

    public Integer reduce(VariantContext map, Integer reduce) {
        return (map == null ? 0 : 1) + reduce;
    }

    private VariantContext parseVar(List<String> features, ReferenceContext reference, GenomeLoc featureLoc) {
        VariantType type = VariantType.match(features.get(4));
        Allele ref;
        Allele alt;
        int end;
        if ( type.equals(VariantType.SNP) ) {
            alt = Allele.create(features.get(3), false);
            ref = Allele.create(reference.getBase(),true);
            if ( ! ( features.get(2).equals("<FILL>") || features.get(2).equals(Character.toString((char)reference.getBase())) ) ) {
                throw new UserException("Allelic mismatch, ref "+((char) reference.getBase())+" and given "+features.get(2));
            }
            end = 1+featureLoc.getStart();
        } else if ( type.equals(VariantType.DELETION) ) {
            int size = Integer.parseInt(features.get(1));
            byte[] refBase = new byte[size];
            for ( int idx = 1; idx <= size; idx++ ) {
                refBase[idx-1] = reference.getBases()[50+idx];
            }
            ref = Allele.create(refBase,true);
            alt = Allele.create(Allele.NULL_ALLELE_STRING,false);
            end = featureLoc.getStart()+size    ;
        } else if ( type.equals(VariantType.INSERTION) ) {
            alt = Allele.create(features.get(3).getBytes());
            ref = Allele.create(Allele.NULL_ALLELE_STRING,true);
            end = featureLoc.getStart();
        } else { // MNP
            int size = Integer.parseInt(features.get(1));
            byte[] refBase = new byte[size];
            for ( int idx = 1; idx <= size; idx ++ ) {
                refBase[idx-1] = reference.getBases()[50+idx];
            }
            ref = Allele.create(refBase,true);
            alt = Allele.create(features.get(3),false);
            end = featureLoc.getStop();
        }

        // String source, String contig, long start, long stop, Collection<Allele> alleles, Map<String, Genotype> genotypes,
        // double negLog10PError, Set<String> filters, Map<String, ?> attributes


        logger.debug(ref);
        logger.debug(alt);

        return new VariantContext("Table2VCF",featureLoc.getContig(),featureLoc.getStart(),end, Arrays.asList(ref,alt),
                new HashMap<String,Genotype>(), 1.0, new HashSet<String>(), new HashMap<String,Object>());
    }

    protected enum VariantType {
        SNP("SNP", VariantContext.Type.SNP),
        INSERTION("INS", VariantContext.Type.INDEL),
        DELETION("DEL", VariantContext.Type.INDEL),
        MNP("MNP", VariantContext.Type.MNP);

        String name;
        VariantContext.Type type;

        VariantType(String name, VariantContext.Type type) {
            this.name = name;
            this.type = type;
        }

        public String getName() { return name; }
        public VariantContext.Type getVariantContextType() { return type; }

        public static VariantType match(String given) {
            for ( VariantType type : VariantType.values() ) {
                if ( given.equals(type.getName()) ) {
                    return type;
                }
            }

            return null;
        }
    }
}
