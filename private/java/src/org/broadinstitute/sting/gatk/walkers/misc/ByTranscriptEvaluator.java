package org.broadinstitute.sting.gatk.walkers.misc;

import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 6/1/12
 * Time: 2:50 PM
 * To change this template use File | Settings | File Templates.
 */
public class ByTranscriptEvaluator extends RodWalker<VariantContext,ByTranscriptEvaluator.EvalContext>  implements TreeReducible<ByTranscriptEvaluator.EvalContext> {
    // todo -- some way to specify how to properly parse String annotations from the header that aren't really strings (but instead float/integer)
    @Input(doc="eval file",shortName="eval",fullName="eval",required=true)
    RodBinding<VariantContext> eval;

    @Argument(doc="Additional keys added to the info field by the annotation engine. Must have a corresponding info header line.",required=false,fullName="keys",shortName="k")
    List<String> additionalKeys = new ArrayList<String>();

    @Argument(doc="Ignore transcripts which are flagged as non-coding for a SNP",required=false,shortName="inc",fullName="ignoreNonCoding")
    boolean ignoreNonCoding = false;

    @Argument(doc="Only use CCDS transcripts",required=false,shortName="ccds",fullName="ccdsOnly")
    boolean ccdsOnly = false;

    @Argument(doc="Minimum allele frequency",required=false,shortName="minAAF",fullName="minAAF")
    double minAAF = 0.0;

    @Argument(doc="Maximum allele frequency",required=false,shortName="maxAAF",fullName="maxAAF")
    double maxAAF = 1.0;

    @Hidden
    @Argument(doc="Allele frequency key",required=false,fullName="afKey")
    String afkey = "AF";

    @Output
    PrintStream out;

    TranscriptInfoParser transcriptInfoParser;

    private static final String TRANSCRIPT_INFO_KEY = "TranscriptInfo";
    private static final String TRANSCRIPT_NAME_KEY = "Feature";
    private static final String GENE_NAME_KEY = "Gene";
    private static final String CONSEQ_NAME_KEY = "Consequence";
    private static final String SIFT_KEY = "SIFT";
    private static final String POLYPHEN_KEY = "PolyPhen";
    private static final String CCDS_KEY = "CCDS";

    private final List<String> REQUESTED_FIELDS = Arrays.asList(new String[]{TRANSCRIPT_NAME_KEY,GENE_NAME_KEY,CONSEQ_NAME_KEY,SIFT_KEY,POLYPHEN_KEY,CCDS_KEY});

    public void initialize() {

        VCFHeader header = VCFUtils.getVCFHeadersFromRods(getToolkit(), Arrays.asList(eval)).get(eval.getName());
        VCFInfoHeaderLine csqFormat = header.getInfoHeaderLine("CSQ");
        Map<String,VCFInfoHeaderLine> additionalFormats = new HashMap<String,VCFInfoHeaderLine>(additionalKeys.size());
        for ( String k : additionalKeys ) {
            additionalFormats.put(k,header.getInfoHeaderLine(k));
        }

        transcriptInfoParser = new TranscriptInfoParser(csqFormat,additionalFormats,REQUESTED_FIELDS);
    }

    public VariantContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext alignment) {
        if ( tracker != null && tracker.hasValues(eval) ) {
            VariantContext myEval = tracker.getFirstValue(eval);
            // todo -- update the filtering clause here to take command line values
            if ( ! myEval.isFiltered() && myEval.isBiallelic() && myEval.isSNP() && myEval.isPolymorphicInSamples() && matchesFrequency(myEval)) {
                logger.debug(ref.getLocus());
                return transcriptInfoParser.parseTranscriptInfo(myEval);
            }
        }

        return null;
    }

    private boolean matchesFrequency(VariantContext vc) {
        if ( ! vc.hasAttribute(afkey)  )
            return false;
        double af = Double.parseDouble(vc.getAttribute(afkey).toString());
        return (af > minAAF && af < maxAAF);
    }

    public EvalContext reduce(VariantContext map, EvalContext prevReduce) {
        if ( map == null ) {
            return prevReduce;
        }


        prevReduce.addContext((Map<String,Set<VariantTranscriptContext>>) map.getAttribute(TRANSCRIPT_INFO_KEY));

        return prevReduce;
    }

    public EvalContext treeReduce(EvalContext left, EvalContext right) {
        left.merge(right);

        return left;
    }

    public EvalContext reduceInit() {
        return new EvalContext();
    }

    public void onTraversalDone(EvalContext context) {
        List<String> colNames = new ArrayList<String>(25);
        colNames.addAll(Arrays.asList(new String[]{"Gene","Transcript","nVariants"}));
        // todo -- logic to dynamically determine what requested statistics to output
        colNames.addAll(Arrays.asList(new String[]{"SYN","NONSYN","STOP","SPLICE","PolyphenDamaging","SiftDeleterious"}));
        GATKReport report = GATKReport.newSimpleReport("TranscriptReport",colNames);
        ArrayList<Object> values = new ArrayList<Object>(colNames.size());
        for ( GeneInfo geneInfo : context.getGenes() ) {
            Set<TranscriptInfo> transcripts = new HashSet<TranscriptInfo>(geneInfo.getTranscripts());
            transcripts.add(geneInfo.worstAnnotationTranscript);
            for ( TranscriptInfo transcriptInfo : transcripts ) {
                // todo -- logic to hook up user-requested statistics to access pattern
                values.clear();
                values.add(geneInfo.geneName);
                values.add(transcriptInfo.transcriptName);
                values.add(transcriptInfo.numVariants);
                values.add(transcriptInfo.getCounts(ConsequenceType.SYNONYMOUS_CODING));
                values.add(transcriptInfo.getCounts(ConsequenceType.NONSYNONYMOUS_CODING));
                values.add(transcriptInfo.getCounts(ConsequenceType.STOP_GAIN)+transcriptInfo.getCounts(ConsequenceType.STOP_LOSS));
                values.add(transcriptInfo.getCounts(ConsequenceType.SPLICE_SITE));
                values.add(transcriptInfo.getPolyphenAbove(0.85));
                values.add(transcriptInfo.getSiftBelow(0.05));
                report.addRowList(values);
            }
        }
        report.print(out);
    }

    class TranscriptInfoParser {

        Map<String,Integer> fieldOffset;

        public TranscriptInfoParser(VCFInfoHeaderLine csqFormat, Map<String,VCFInfoHeaderLine> extraFields, List<String> requestedFields) {
            String fieldStr = csqFormat.getDescription().replaceFirst("Consequence type as predicted by VEP. Format: ","");
            // fieldStr should be of the format KEY1|KEY2|KEY3|KEY4|...
            String[] fields = fieldStr.split("\\|");
            fieldOffset = new HashMap<String,Integer>(requestedFields.size());
            // todo -- this is n^2 but nobody cares since it's only done once and the fields aren't really that big, but n^2 is kinda ugly
            for ( String rf : requestedFields ) {
                fieldOffset.put(rf,ArrayUtils.indexOf(fields,rf));
            }
        }

        public Map<String,Set<VariantTranscriptContext>> parse(String CSQvalue) {
            String[] tVals = CSQvalue.split(",");
            Map<String,Set<VariantTranscriptContext>> toRet = new HashMap<String,Set<VariantTranscriptContext>>(tVals.length);
            for ( String tval : tVals ) {
                String[] fields = tval.split("\\|");
                VariantTranscriptContext vtc = new VariantTranscriptContext();
                vtc.setGeneName(fields[fieldOffset.get(GENE_NAME_KEY)]);
                if ( vtc.getGeneName().equals("") )
                    continue;
                vtc.setTranscriptName(fields[fieldOffset.get(TRANSCRIPT_NAME_KEY)]);
                String ccds = "";
                if ( fieldOffset.get(CCDS_KEY) < fields.length)
                    ccds = fields[fieldOffset.get(CCDS_KEY)];
                if ( ccds.equals("") && ccdsOnly ) {
                    logger.debug(ccds);
                    continue;
                }

                // if specified, ignore genes that are flagged as not coding
                Set<ConsequenceType> consequences = ConsequenceType.decode(fields[fieldOffset.get(CONSEQ_NAME_KEY)]);
                if ( consequences.contains(ConsequenceType.GENE_NOT_CODING) && ignoreNonCoding ) {
                    continue;
                }
                vtc.setConsequences(consequences);

                if ( fields.length > Math.max(fieldOffset.get(SIFT_KEY),fieldOffset.get(POLYPHEN_KEY))) {
                    String sift = fields[fieldOffset.get(SIFT_KEY)];
                    String polyphen = fields[fieldOffset.get(POLYPHEN_KEY)];
                    if ( ! sift.equals("") ) {
                        vtc.setSiftScore(parseSift(sift));
                    }
                    if ( ! polyphen.equals("") ) {
                        vtc.setPolyphenScore(parseSift(polyphen));
                    }
                }
                if ( ! toRet.containsKey(vtc.getGeneName()) ) {
                    toRet.put(vtc.getGeneName(),new HashSet<VariantTranscriptContext>(12));
                }
                toRet.get(vtc.getGeneName()).add(vtc);
            }

            return toRet;
        }

        public VariantContext parseTranscriptInfo(VariantContext context) {
            VariantContextBuilder builder = new VariantContextBuilder(context);
            builder.attribute(TRANSCRIPT_INFO_KEY,this.parse(context.getAttribute("CSQ").toString()));
            return builder.make();
        }

        private double parseSift(String sift) {
            // string of form EFFECT(number)
            return Double.parseDouble(sift.split("\\(")[1].replace(")",""));
        }
    }

    class VariantTranscriptContext {

        private String transcriptName;
        private String geneName;
        private Set<ConsequenceType> consequences;
        private Double siftScore;
        private Double polyphenScore;

        public String getTranscriptName() {
            return transcriptName;
        }

        public String getGeneName() {
            return geneName;
        }

        public void setTranscriptName(String name) {
            transcriptName = name;
        }

        public void setGeneName(String name) {
            geneName = name;
        }

        public void setConsequences(Set<ConsequenceType> types) {
            consequences = types;
        }

        public Set<ConsequenceType> getConsequences() {
            return consequences;
        }

        public void setSiftScore(double sift) {
            siftScore = sift;
        }

        public void setPolyphenScore(double polyphen) {
            polyphenScore = polyphen;
        }

        public boolean hasSiftScore() {
            return siftScore != null;
        }

        public boolean hasPolyphenScore() {
            return polyphenScore != null;
        }

        public double getSiftScore() {
            return siftScore;
        }

        public double getPolyphenScore() {
            return polyphenScore;
        }
    }

    class GeneInfo {
        private String geneName;
        private Map<String,TranscriptInfo> transcripts;
        private TranscriptInfo worstAnnotationTranscript;
        private long numVariants;

        public GeneInfo(String name) {
            geneName = name;
            numVariants = 0l;
            transcripts = new HashMap<String,TranscriptInfo>(16);
            worstAnnotationTranscript = new TranscriptInfo(geneName+"_worst");
        }

        public void addContexts(Set<VariantTranscriptContext> contexts) {
            ConsequenceType worstType = ConsequenceType.INTERGENIC;
            double polyPhen = Double.MIN_VALUE;
            double sift = Double.MAX_VALUE;
            for ( VariantTranscriptContext context : contexts ) {
                if ( ! transcripts.containsKey(context.getTranscriptName()) ) {
                    addTranscript(context.getTranscriptName());
                }

                transcripts.get(context.getTranscriptName()).addContext(context);
                for ( ConsequenceType type : context.getConsequences() ) {
                    if ( type.compareTo(worstType) < 0 ) {
                        worstType = type;
                    }
                }
                if ( context.hasPolyphenScore() && context.getPolyphenScore() > polyPhen ) {
                    polyPhen = context.getPolyphenScore();
                }
                if ( context.hasSiftScore() && context.getSiftScore() < sift ) {
                    sift = context.getSiftScore();
                }
            }

            worstAnnotationTranscript.numVariants++;
            if ( worstType.compareTo(ConsequenceType.INTRON) <= 0 ) {
                worstAnnotationTranscript.addConsequence(worstType);
                if ( polyPhen > -1.0 && sift > -1.0 ) {
                    worstAnnotationTranscript.polyphenScores.add(polyPhen);
                    worstAnnotationTranscript.siftScores.add(sift);
                }
            }
            //logger.debug(numVariants);
            numVariants++;
        }

        private void addTranscript(String tName) {
            transcripts.put(tName,new TranscriptInfo(tName));
        }

        public void merge(GeneInfo other) {
            if ( ! other.geneName.equals(this.geneName) ) {
                throw new IllegalStateException("Gene info objects referencing different genes can not be merged");
            }
            for ( Map.Entry<String,TranscriptInfo> info : other.transcripts.entrySet() ) {
                if ( transcripts.containsKey(info.getKey()) ) {
                    transcripts.get(info.getKey()).merge(info.getValue());
                } else {
                    transcripts.put(info.getKey(),info.getValue());
                }
            }
            numVariants += other.numVariants;
        }

        public Collection<TranscriptInfo> getTranscripts() {
            return transcripts.values();
        }
    }

    class TranscriptInfo {

        private String transcriptName;
        private long numVariants;
        private Map<ConsequenceType,Integer> consequenceCounts;
        private List<Double> polyphenScores;
        private List<Double> siftScores;

        public TranscriptInfo(String name) {
            transcriptName = name;
            numVariants = 0;
            consequenceCounts = new HashMap<ConsequenceType,Integer>(ConsequenceType.values().length);
            polyphenScores = new ArrayList<Double>(16);
            siftScores = new ArrayList<Double>(16);
        }

        public void merge(TranscriptInfo other) {
            // todo -- right now do nothing (names will be the same)
            this.numVariants += other.numVariants;
            this.polyphenScores.addAll(other.polyphenScores);
            this.siftScores.addAll(other.siftScores);
            for ( Map.Entry<ConsequenceType,Integer> otherEtry : other.consequenceCounts.entrySet() ) {
                if ( this.consequenceCounts.containsKey(otherEtry.getKey()) ) {
                    this.consequenceCounts.put(otherEtry.getKey(),this.consequenceCounts.get(otherEtry.getKey())+otherEtry.getValue());
                } else {
                    this.consequenceCounts.put(otherEtry.getKey(),otherEtry.getValue());
                }
            }
        }

        public void addContext(VariantTranscriptContext context) {
            numVariants++;
            for ( ConsequenceType type : context.getConsequences() ) {
                if ( ! type.describesTranscript )
                    addConsequence(type);
            }
            if ( context.hasPolyphenScore() ) {
                polyphenScores.add(context.getPolyphenScore());
            }
            if (context.hasSiftScore()) {
                siftScores.add(context.getSiftScore());
            }
        }

        private void addConsequence(ConsequenceType type) {
            if ( ! consequenceCounts.containsKey(type) ) {
                consequenceCounts.put(type,0);
            }
            consequenceCounts.put(type,consequenceCounts.get(type)+1);
        }

        public Integer getPolyphenAbove(double score) {
            int ct = 0;
            for ( double d : polyphenScores ) {
                ct += d >= score ? 1 : 0;
            }

            return ct;
        }

        public Integer getSiftBelow(double score) {
            int ct = 0;
            for ( double d : siftScores ) {
                ct += d <= score ? 1 : 0;
            }

            return ct;
        }

        public Integer getCounts(ConsequenceType type) {
            return consequenceCounts.containsKey(type) ? consequenceCounts.get(type) : 0;
        }

        public String toString() {
            StringBuilder builder = new StringBuilder();
            builder.append(transcriptName);
            builder.append("\t");
            builder.append("numVar:");
            builder.append(numVariants);
            for ( ConsequenceType type : ConsequenceType.values() ) {
                if ( type.describesTranscript )
                    continue;
                builder.append("\t");
                builder.append(type.ensembleTerm);
                builder.append(":");
                Integer num = consequenceCounts.get(type);
                if ( num == null ) {
                    num = 0;
                }
                builder.append(num);
            }
            builder.append("\t");
            builder.append("polyphen: ");
            builder.append(Utils.join(",",polyphenScores));
            builder.append("\tsift: ");
            builder.append(Utils.join(",",siftScores));
            return builder.toString();
        }
    }

    class EvalContext {
        private Map<String,GeneInfo> geneInfoMap;

        public EvalContext() {
            geneInfoMap = new HashMap<String,GeneInfo>(1024);
        }

        public boolean hasGene(String name) {
            return geneInfoMap.keySet().contains(name);
        }

        public void addContext(Map<String,Set<VariantTranscriptContext>> context) {
            for ( Map.Entry<String,Set<VariantTranscriptContext>> geneVTC : context.entrySet() ) {
                if ( ! geneInfoMap.containsKey(geneVTC.getKey()) ) {
                    geneInfoMap.put(geneVTC.getKey(),new GeneInfo(geneVTC.getKey()));
                }
                geneInfoMap.get(geneVTC.getKey()).addContexts(geneVTC.getValue());
            }
        }

        public void merge(EvalContext other) {
            for ( Map.Entry<String,GeneInfo> infoEntry : other.geneInfoMap.entrySet() ) {
                if ( geneInfoMap.containsKey(infoEntry.getKey()) ) {
                    geneInfoMap.get(infoEntry.getKey()).merge(infoEntry.getValue());
                } else {
                    geneInfoMap.put(infoEntry.getKey(),infoEntry.getValue());
                }
            }
        }

        public Collection<GeneInfo> getGenes() {
            return geneInfoMap.values();
        }
    }

    enum ConsequenceType {


        STOP_GAIN("Stop gained","STOP_GAINED"),
        STOP_LOSS("Stop lost","STOP_LOST"),
        FRAMESHIFT("Frameshift coding","FRAMESHIFT_CODING"),
        ESSENTIAL_SPLICE("Essential splice site - both donor and acceptor","ESSENTIAL_SPLICE_SITE"),
        SPLICE_SITE("1-3 bps into an exon or 3-8 bps into an intron","SPLICE_SITE"),
        NONSYNONYMOUS_CODING("Nonsynonymous coding. Includes codon change, codon loss, and codon gain.","NON_SYNONYMOUS_CODING"),
        TF_BINDING("Transcription factor binding motif","TRANSCRIPTION_FACTOR_BINDING_MOTIF"),
        REGULATORY("Regulatory region","REGULATORY_REGION"),
        SYNONYMOUS_CODING("Synonymous coding - includes both stop and other codons","SYNONYMOUS_CODING"),
        INTRON("Intronic","INTRONIC"),
        PRIME_5("5 prime UTR","5PRIME_UTR"),
        PRIME_3("3 prime UTR","3PRIME_UTR"),
        UPSTREAM("upstream - within 5KB","UPSTREAM"),
        DOWNSTREAM("downstream - within 5KB","DOWNSTREAM"),
        COMPLEX("Complex in/del","COMPLEX_INDEL"),
        PARTIAL_CODON("Partial codon","PARTIAL_CODON"),
        CODING_UNKNOWN("Coding unknown","CODING_UNKNOWN"),
        MIRNA("Within mature miRNA","WITHIN_MATURE_miRNA"),
        NMD("NMD transcript","NMD_TRANSCRIPT",true),
        GENE_NOT_CODING("Within non-coding gene","WITHIN_NON_CODING_GENE",true),
        INTERGENIC("Intergenic","INTERGENIC");


        static Map<String,ConsequenceType> parsingMap = new HashMap<String,ConsequenceType>(ConsequenceType.values().length);
        // i know this is dumb, and you should use this.valueOf, but you can't have an enum named "5_PRIME_UTR" as it starts with a number.
        // so this is the only way. Annoying.

        private String description;
        private String ensembleTerm;
        private boolean describesTranscript;

        ConsequenceType(String desc, String term, boolean describesT) {
            description = desc;
            ensembleTerm = term;
            describesTranscript = describesT;
        }

        ConsequenceType(String desc, String term) {
            this(desc,term,false);
        }

        public static Set<ConsequenceType> decode(String fromVC) {
            if ( parsingMap.size() == 0 ) {
                for ( ConsequenceType type : ConsequenceType.values() ) {
                    parsingMap.put(type.ensembleTerm,type);
                }
            }
            String[] consequenceTypes = fromVC.split("\\&");
            Set<ConsequenceType> matching = new HashSet<ConsequenceType>(consequenceTypes.length);
            for ( String ct : consequenceTypes ) {
                matching.add(parsingMap.get(ct));
            }
            return matching;
        }

    }

}
