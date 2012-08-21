package org.broadinstitute.sting.gatk.walkers.misc;

import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.EvaluationContext;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

import java.io.File;
import java.io.IOException;
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

    @Input(doc="bootstrap file",shortName="boot",fullName = "boot",required=false)
    RodBinding<VariantContext> bootBinding = null;

    @Argument(doc="bootSample",shortName="bs",fullName="bootSample",required=false)
    String bootSam = null;

    @Argument(doc="Additional keys added to the info field by the annotation engine. Must have a corresponding info header line.",required=false,fullName="keys",shortName="k")
    List<String> additionalKeys = new ArrayList<String>();

    @Argument(doc="Ignore transcripts which are flagged as non-coding for a SNP",required=false,shortName="inc",fullName="ignoreNonCoding")
    boolean ignoreNonCoding = false;

    @Argument(doc="Only use CCDS transcripts",required=false,shortName="ccds",fullName="ccdsOnly")
    boolean ccdsOnly = false;

    @Argument(doc="Only use Refseq transcripts",required=false,shortName="refseq",fullName="refseqOnly")
    boolean refseqOnly = false;

    @Argument(doc="Only use ENSEMBL transcripts",required=false,shortName="ensembl",fullName="ensemblOnly")
    boolean ensemblOnly = false;

    @Argument(doc="Minimum allele frequency",required=false,shortName="minAAF",fullName="minAAF")
    double minAAF = 0.0;

    @Argument(doc="Maximum allele frequency",required=false,shortName="maxAAF",fullName="maxAAF")
    double maxAAF = 1.0;

    @Argument(doc="Ignore these transcripts (e.g. for being spurious/noncoding)",required=false,shortName="xt",fullName="excludeTranscripts")
    File ignoreTranscriptFile = null;

    @Hidden
    @Argument(doc="Allele frequency key",required=false,fullName="afKey")
    List<String> afkey = Arrays.asList(new String[]{"AF"});

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

    private long nSeen = 0l;
    private long nProcessed = 0l;
    public void initialize() {
        assertCodeIsWorking();
        assertInputsAreGood();
        transcriptInfoParser = initializeTranscriptParser();
    }

    public VariantContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext alignment) {
        if ( tracker != null && tracker.hasValues(eval) ) {
            VariantContext myEval = tracker.getFirstValue(eval);
            // todo -- update the filtering clause here to take command line values
            if ( ! bootBinding.isBound() ) {
                if ( ! myEval.isFiltered() && myEval.isSNP() ) {
                    nSeen++;
                    return transcriptInfoParser.addTranscriptInformationToVC(myEval);
                }
            } else {
                if ( getBootstrapAC(tracker,bootBinding,bootSam) > 0 ) {
                    nSeen++;
                    return transcriptInfoParser.addTranscriptInformationToVC(myEval);
                } else {
                    return null;
                }
            }
        }

        return null;
    }

    public EvalContext reduce(VariantContext map, EvalContext prevReduce) {
        if ( map == null ) {
            return prevReduce;
        }
        Map<String,Set<VariantTranscriptContext>> contextsByGene = (Map<String,Set<VariantTranscriptContext>>) map.getAttribute(TRANSCRIPT_INFO_KEY);
        if ( contextsByGene.size() > 0 ) {
            prevReduce.addContext(contextsByGene);
            nProcessed++;
        }

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
        colNames.addAll(Arrays.asList(new String[]{"Gene","Transcript","nVariants","nCoding"}));
        // todo -- logic to dynamically determine what requested statistics to output
        colNames.addAll(Arrays.asList(new String[]{"SYN","NONSYN","STOP","SPLICE","PolyphenDamaging","SiftDeleterious"}));
        GATKReport transcriptReport = GATKReport.newSimpleReport("TranscriptReport",colNames);
        for ( GeneInfo geneInfo : context.getGenes() ) {

            for ( TranscriptInfo transcriptInfo : geneInfo.getTranscripts() ) {
                // todo -- logic to hook up user-requested statistics to access pattern
                ArrayList<Object> transValues = getTranscriptOutputValues(geneInfo,transcriptInfo,colNames);
                transcriptReport.addRowList(transValues);
            }

            if ( geneInfo.hasMultipleTranscripts() ) {
                ArrayList<Object> mostDamagingValues = getTranscriptOutputValues(geneInfo,geneInfo.mostDeleteriousAnnotationTranscript,colNames);
                ArrayList<Object> leastDamagingValues = getTranscriptOutputValues(geneInfo,geneInfo.leastDeleteriousAnnotationTranscript,colNames);
                transcriptReport.addRowList(mostDamagingValues);
                transcriptReport.addRowList(leastDamagingValues);
            }
        }
        transcriptReport.print(out);
        out.printf("%n");
        GATKReport variantReport = GATKReport.newSimpleReport("VariantReport",Arrays.asList(new String[]{"NumberOfTranscripts","NumCodingVariants"}));
        Map<Integer,Integer> variantCountsByNTranscript = countVariantsByNTranscript(context);
        for (Map.Entry<Integer,Integer> countByNTranscript : variantCountsByNTranscript.entrySet() ) {
            variantReport.addRow(countByNTranscript.getKey(),countByNTranscript.getValue());
        }
        variantReport.print(out);

        logger.info(String.format("Traversal Result: Num_Seen=%d and Num_Processed=%d",nSeen,nProcessed));
    }

    private void assertInputsAreGood() {
        if ( refseqOnly && ccdsOnly || refseqOnly && ensemblOnly || ensemblOnly && ccdsOnly) {
            throw new UserException("Cannot exclusively use multiple of specific transcripts (ensembl,refseq,ccds). Choose exactly one and not multiple.");
        }

        if ( bootSam == null && bootBinding.isBound() ) {
            throw new UserException("Please provide a specific bootstrap sample");
        }
    }

    public int getBootstrapAC(RefMetaDataTracker tracker, RodBinding<VariantContext> binding, String sample) {
        if ( ! tracker.hasValues(bootBinding) ) {
            return -1;
        }
        VariantContext bootstrap = tracker.getFirstValue(binding);
        Genotype bootGeno = bootstrap.getGenotype(sample);
        Object bac = bootGeno.getAnyAttribute("BAC");
        return parseBAC(bac);
    }

    private int parseBAC(Object bac) {
        String[] sp = bac.toString().split(",");
        int b = 0;
        for ( String s : sp ) {
            int t = Integer.parseInt(s);
            if ( t > b ) {
                b = t;
            }
        }

        return b;
    }

    private Map<Integer,Integer> countVariantsByNTranscript(EvalContext evalContext) {
        Map<Integer,Integer> countsByNTranscript = new HashMap<Integer,Integer>(16);
        for ( String genesOverlappingVariantStr : evalContext.variantGeneList ) {
            logger.debug(genesOverlappingVariantStr);
            String[] genesOverlappingVariant = genesOverlappingVariantStr.split(",");
            int maxNTranscript = 0;
            for ( String gene : genesOverlappingVariant ) {
                int nTranscript = evalContext.geneInfoMap.get(gene).getTranscripts().size();
                if ( nTranscript > maxNTranscript )
                    maxNTranscript = nTranscript;
            }

            // cutoff at 10+ transcripts
            if ( maxNTranscript > 10 )
                maxNTranscript = 10;

            if (!  countsByNTranscript.containsKey(maxNTranscript) ) {
                countsByNTranscript.put(maxNTranscript,0);
            }

            int counts = countsByNTranscript.get(maxNTranscript);
            countsByNTranscript.put(maxNTranscript,1+counts);
        }

        return countsByNTranscript;
    }

    private TranscriptInfoParser initializeTranscriptParser() {
        VCFHeader header = VCFUtils.getVCFHeadersFromRods(getToolkit(), Arrays.asList(eval)).get(eval.getName());
        VCFInfoHeaderLine csqFormat = header.getInfoHeaderLine("CSQ");
        Map<String,VCFInfoHeaderLine> additionalFormats = new HashMap<String,VCFInfoHeaderLine>(additionalKeys.size());
        for ( String k : additionalKeys ) {
            additionalFormats.put(k,header.getInfoHeaderLine(k));
        }

        Set<String> ignoreTranscripts = readTranscriptsToIgnore(ignoreTranscriptFile);

        return new TranscriptInfoParser(csqFormat,additionalFormats,REQUESTED_FIELDS,ignoreTranscripts);
    }

    private Set<String> readTranscriptsToIgnore(File ignoreFile) {
        if ( ignoreFile == null )
            return new HashSet<String>(0);

        Set<String> badTranscripts = new HashSet<String>(100000);

        try {
            for ( String transcript : new XReadLines(ignoreFile) ) {
                badTranscripts.add(transcript);
            }
        } catch (IOException ioException ) {
            throw new UserException("Error opening file for reading: "+ignoreFile.getAbsolutePath());
        }

        return badTranscripts;
    }

    private ArrayList<Object> getTranscriptOutputValues(GeneInfo geneInfo, TranscriptInfo transcriptInfo, List<String> columnNames) {
        ArrayList<Object> values = new ArrayList<Object>(columnNames.size());
        values.add(geneInfo.geneName);
        values.add(transcriptInfo.transcriptName);
        values.add(transcriptInfo.numVariants);
        values.add(getNumCodingVariants(transcriptInfo));
        values.add(transcriptInfo.getCounts(ConsequenceType.SYNONYMOUS_CODING));
        values.add(transcriptInfo.getCounts(ConsequenceType.NONSYNONYMOUS_CODING));
        values.add(transcriptInfo.getCounts(ConsequenceType.STOP_GAIN)+transcriptInfo.getCounts(ConsequenceType.STOP_LOSS));
        values.add(transcriptInfo.getCounts(ConsequenceType.SPLICE_SITE)+transcriptInfo.getCounts(ConsequenceType.ESSENTIAL_SPLICE));
        values.add(transcriptInfo.getPolyphenAbove(0.75));
        values.add(transcriptInfo.getSiftBelow(0.10));
        return values;
    }

    private int getNumCodingVariants(TranscriptInfo transcriptInfo) {
        int nCoding = 0;
        for ( Map.Entry<ConsequenceType,Integer> consequenceCounts : transcriptInfo.consequenceCounts.entrySet() ) {
            if ( ConsequenceType.isCoding(consequenceCounts.getKey()) ) {
                nCoding += consequenceCounts.getValue();
            }
        }
        return nCoding;
    }

    private void assertCodeIsWorking() {
        if ( ConsequenceType.isCoding(ConsequenceType.GENE_NOT_CODING) )
            throw new StingException("The ConsequenceType Enum is busted. GENE_NOT_CODING should not be coding.");

        if ( ! ConsequenceType.isCoding(ConsequenceType.SYNONYMOUS_CODING) )
            throw new StingException("The ConsequenceType Enum is busted! Syn_Coding should be coding.");
    }

/*    // note: from previous incarnation that did not use bootstraps
    @Deprecated
    private boolean matchesFrequency(VariantContext vc) {
        String afk = null;
        for ( String k : afkey ) {
            if ( vc.hasAttribute(k) ) {
                afk = k;
                break;
            }
        }
        if ( afk == null  )
            return false;

        if ( ! vc.isBiallelic() ) {
            Object afObject = vc.getAttribute(afk);
            double  af = 0.0;
            if ( afObject instanceof  String) {
                String[] afstring = ((String) afObject).split(",");
                for ( String s : afstring ) {
                    af += Double.parseDouble(s);
                }
            } else {
                List<Object> afList = (List<Object>) afObject;
                for ( Object o : afList ) {
                    if (o instanceof  String ) {
                        af += Double.parseDouble((String) o );
                    } else {
                        af += (Double) o;
                    }
                }

            }

            return (af > minAAF && af < maxAAF);
        }
        double af = Double.parseDouble(vc.getAttribute(afk).toString());
        return (af > minAAF && af < maxAAF);
    }*/

    class TranscriptInfoParser {

        Map<String,Integer> fieldOffset;
        final Set<String> ignoreTranscripts;

        public TranscriptInfoParser(VCFInfoHeaderLine csqFormat, Map<String,VCFInfoHeaderLine> extraFields, List<String> requestedFields, Set<String> ignoreTrans) {
            String fieldStr = csqFormat.getDescription().replaceFirst("Consequence type as predicted by VEP. Format: ","");
            // fieldStr should be of the format KEY1|KEY2|KEY3|KEY4|...
            String[] fields = fieldStr.split("\\|");
            fieldOffset = new HashMap<String,Integer>(requestedFields.size());
            for ( String rf : requestedFields ) {
                fieldOffset.put(rf,ArrayUtils.indexOf(fields,rf));
            }

            ignoreTranscripts = ignoreTrans;
        }

        private String getGeneName(String[] transcriptFields) {
            return transcriptFields[fieldOffset.get(GENE_NAME_KEY)];
        }

        private String getTranscriptName(String[] transcriptFields) {
            return transcriptFields[fieldOffset.get(TRANSCRIPT_NAME_KEY)];
        }

        private String getCCDS_ID(String[] transcriptFields) {
            if ( fieldOffset.get(CCDS_KEY) >= transcriptFields.length ) {
                return null;
            }

            String ccdsID = transcriptFields[fieldOffset.get(CCDS_KEY)];

            if ( ccdsID.equals("") )
                return null;

            return ccdsID;
        }

        private Set<ConsequenceType> decodeConsequences(String[] transcriptFields) {
            return ConsequenceType.decode(transcriptFields[fieldOffset.get(CONSEQ_NAME_KEY)]);
        }

        private boolean isRefseqTranscript(VariantTranscriptContext context) {
            return context.getTranscriptName().startsWith("NM_");
        }

        private boolean isEnsemblTranscript(VariantTranscriptContext context) {
            return context.getTranscriptName().startsWith("ENS");
        }

        private boolean isCCDSTranscript(VariantTranscriptContext context) {
            return context.getCCDS_ID() != null;
        }

        private VariantTranscriptContext setTranscriptIDs(String[] transcriptFields,VariantTranscriptContext context) {
            context.setGeneName(getGeneName(transcriptFields));
            context.setTranscriptName(getTranscriptName(transcriptFields));
            context.setCCDS_ID(getCCDS_ID(transcriptFields));
            return context;
        }

        private VariantTranscriptContext setSIFTscore(String[] transFields, VariantTranscriptContext context) {
            if ( transFields.length > fieldOffset.get(SIFT_KEY) ) {
                String siftStr = transFields[fieldOffset.get(SIFT_KEY)];
                if ( ! siftStr.equals("") ) {
                    double sift = parseSiftOrPolyphen(siftStr);
                    context.setSiftScore(sift);
                }
            }

            return context;
        }

        private VariantTranscriptContext setPolyphenScore(String[] transFields, VariantTranscriptContext context) {
            if ( transFields.length > fieldOffset.get(POLYPHEN_KEY) ) {
                String polyStr = transFields[fieldOffset.get(POLYPHEN_KEY)];
                if ( ! polyStr.equals("") ) {
                    double polyphen = parseSiftOrPolyphen(polyStr);
                    context.setPolyphenScore(polyphen);
                }
            }

            return context;
        }

        private boolean isInSpecifiedDB(VariantTranscriptContext context) {
            return  ! ( context.geneName == null || context.geneName.equals("") ||
                    ( refseqOnly && ! isRefseqTranscript(context) ) ||
                    ( ensemblOnly && ! isEnsemblTranscript(context) ) ||
                    ( ccdsOnly && ! isCCDSTranscript(context) ) );
        }

        private VariantTranscriptContext setConsequences(String[] transcriptFields, VariantTranscriptContext context) {
            Set<ConsequenceType> consequences = decodeConsequences(transcriptFields);
            context.setConsequences(consequences);
            return context;
        }

        private VariantTranscriptContext setConservationScores(String[] transcriptFields, VariantTranscriptContext context) {
            VariantTranscriptContext vtc = setSIFTscore(transcriptFields,context);
            vtc = setPolyphenScore(transcriptFields,vtc);
            return vtc;
        }

        private boolean transcriptOrGeneNotCoding(VariantTranscriptContext context) {
            return context.getConsequences().contains(ConsequenceType.GENE_NOT_CODING);
        }

        private boolean ignoreTranscript(VariantTranscriptContext context) {
            return ignoreTranscripts.contains(context.getTranscriptName());
        }

        public Map<String,Set<VariantTranscriptContext>> parse(String CSQvalue) {
            String[] transInfoByTranscript = CSQvalue.split(",");
            Map<String,Set<VariantTranscriptContext>> contextsByTranscriptName = new HashMap<String,Set<VariantTranscriptContext>>(transInfoByTranscript.length);
            for ( String tval : transInfoByTranscript ) {
                String[] fields = tval.split("\\|");
                VariantTranscriptContext vtc = new VariantTranscriptContext();
                vtc = setTranscriptIDs(fields,vtc);

                if ( ! isInSpecifiedDB(vtc) )
                    continue;

                if ( ignoreTranscript(vtc) )
                    continue;

                vtc = setConsequences(fields,vtc);

                if ( transcriptOrGeneNotCoding(vtc) && ignoreNonCoding )
                    continue;

                vtc = setConservationScores(fields,vtc);

                if ( ! contextsByTranscriptName.containsKey(vtc.getGeneName()) ) {
                    contextsByTranscriptName.put(vtc.getGeneName(),new HashSet<VariantTranscriptContext>(12));
                }
                contextsByTranscriptName.get(vtc.getGeneName()).add(vtc);
            }

            return contextsByTranscriptName;
        }

        public VariantContext addTranscriptInformationToVC(VariantContext context) {
            VariantContextBuilder builder = new VariantContextBuilder(context);
            builder.attribute(TRANSCRIPT_INFO_KEY,this.parse(context.getAttribute("CSQ").toString()));
            return builder.make();
        }

        private double parseSiftOrPolyphen(String sift) {
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
        private String CCDSid;

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

        public void setCCDS_ID(String id) {
            CCDSid = id;
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

        public String getCCDS_ID() {
            return CCDSid;
        }
    }

    class GeneInfo {
        private String geneName;
        private Map<String,TranscriptInfo> transcripts;
        private TranscriptInfo mostDeleteriousAnnotationTranscript;
        private TranscriptInfo leastDeleteriousAnnotationTranscript;
        private long numVariants;

        public GeneInfo(String name) {
            geneName = name;
            numVariants = 0l;
            transcripts = new HashMap<String,TranscriptInfo>(16);
            mostDeleteriousAnnotationTranscript = new TranscriptInfo(geneName+"_mostDeleterious");
            leastDeleteriousAnnotationTranscript = new TranscriptInfo(geneName+"_leastDeleterious");
        }

        public void addContexts(Set<VariantTranscriptContext> contexts) {
            ConsequenceType mostDeleterious = null;
            ConsequenceType leastDeleterious = null;
            double polyPhen = Double.MIN_VALUE;
            double sift = Double.MAX_VALUE;

            for ( VariantTranscriptContext context : contexts ) {
                if ( ! transcripts.containsKey(context.getTranscriptName()) ) {
                    addTranscript(context.getTranscriptName());
                }

                transcripts.get(context.getTranscriptName()).addContext(context);

                for ( ConsequenceType type : context.getConsequences() ) {

                    if ( ConsequenceType.inCodingRegion(type) && ( mostDeleterious == null || type.compareTo(mostDeleterious) < 0 )) {
                        mostDeleterious = type;
                    }

                    if ( ConsequenceType.inCodingRegion(type) && (leastDeleterious == null || type.compareTo(leastDeleterious) > 0)) {
                        leastDeleterious = type;
                    }
                }

                if ( context.hasPolyphenScore() && context.getPolyphenScore() > polyPhen ) {
                    polyPhen = context.getPolyphenScore();
                }

                if ( context.hasSiftScore() && context.getSiftScore() < sift ) {
                    sift = context.getSiftScore();
                }
            }

            updateConsensusTranscripts(mostDeleterious,leastDeleterious,polyPhen,sift);

            numVariants++;
        }

        private void updateConsensusTranscripts(ConsequenceType mostDeleterious, ConsequenceType leastDeleterious, Double polyPhen, Double sift) {
            if ( mostDeleterious != null && ConsequenceType.inCodingRegion(mostDeleterious) )
                updateConsensus(mostDeleteriousAnnotationTranscript, mostDeleterious, polyPhen,sift);
            if ( leastDeleterious != null && ConsequenceType.inCodingRegion(leastDeleterious) )
                updateConsensus(leastDeleteriousAnnotationTranscript,leastDeleterious,polyPhen,sift);
        }

        private void updateConsensus(TranscriptInfo consensus,ConsequenceType consequence, Double polyPhen, Double sift) {
            consensus.numVariants++;
            consensus.addConsequence(consequence);
            if ( polyPhen > -1.0 && sift > -1.0 ) {
                consensus.polyphenScores.add(polyPhen);
                consensus.siftScores.add(sift);
            }
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

            mostDeleteriousAnnotationTranscript.merge(other.mostDeleteriousAnnotationTranscript);
            leastDeleteriousAnnotationTranscript.merge(other.leastDeleteriousAnnotationTranscript);
            numVariants += other.numVariants;
        }

        public Collection<TranscriptInfo> getTranscripts() {
            return transcripts.values();
        }

        public boolean hasMultipleTranscripts() {
            return getTranscripts().size() > 1;
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
        private ArrayList<String> variantGeneList;

        public EvalContext() {
            geneInfoMap = new HashMap<String,GeneInfo>(4096);
            variantGeneList = new ArrayList<String>(4194304);
        }

        public boolean hasGene(String name) {
            return geneInfoMap.keySet().contains(name);
        }

        public void addContext(Map<String,Set<VariantTranscriptContext>> contextByGene) {
            Set<String> genesVariantIsCodingIn = new HashSet<String>(8);
            for ( Map.Entry<String,Set<VariantTranscriptContext>> geneVTC : contextByGene.entrySet() ) {
                if ( ! geneInfoMap.containsKey(geneVTC.getKey()) ) {
                    geneInfoMap.put(geneVTC.getKey(),new GeneInfo(geneVTC.getKey()));
                }
                geneInfoMap.get(geneVTC.getKey()).addContexts(geneVTC.getValue());
                if ( isCodingInSomeTranscript(geneVTC.getValue()) )
                    genesVariantIsCodingIn.add(geneVTC.getKey());
            }
            if ( genesVariantIsCodingIn.size() > 0 )
                variantGeneList.add(Utils.join(",",genesVariantIsCodingIn));
        }

        private boolean isCodingInSomeTranscript(Set<VariantTranscriptContext> vtContexts) {
            for ( VariantTranscriptContext vtContext : vtContexts ) {
                for ( ConsequenceType consequence : vtContext.getConsequences() ) {
                    if ( ConsequenceType.isCoding(consequence) )
                        return true;
                }
            }

            return false;
        }

        public void merge(EvalContext other) {
            for ( Map.Entry<String,GeneInfo> infoEntry : other.geneInfoMap.entrySet() ) {
                if ( geneInfoMap.containsKey(infoEntry.getKey()) ) {
                    geneInfoMap.get(infoEntry.getKey()).merge(infoEntry.getValue());
                } else {
                    geneInfoMap.put(infoEntry.getKey(),infoEntry.getValue());
                }
            }
            variantGeneList.addAll(other.variantGeneList);
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
        SYNONYMOUS_CODING("Synonymous coding - includes both stop and other codons","SYNONYMOUS_CODING"),
        INTRON("Intronic","INTRONIC"),
        TF_BINDING("Transcription factor binding motif","TRANSCRIPTION_FACTOR_BINDING_MOTIF"),
        PRIME_5("5 prime UTR","5PRIME_UTR"),
        PRIME_3("3 prime UTR","3PRIME_UTR"),
        REGULATORY("Regulatory region","REGULATORY_REGION"),
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
        // i know this is dumb, and you should use this.valueOf, but you can't have an enum named "5PRIME_UTR" as it starts with a number.
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

        public static boolean isCoding(ConsequenceType type) {
            return type.compareTo(INTRON) < 0;
        }

        public static boolean inCodingRegion(ConsequenceType type) {
            return type.compareTo(INTERGENIC) < 0;
        }

    }

}
