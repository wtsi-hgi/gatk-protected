/*
 * Copyright (c) 2011, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.pipeline;

import net.sf.picard.io.IoUtil;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.*;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

/**
 * Standalone program to create an XML index of the GSA pipeline directory tree suitable for loading in IGV.
 *
 * Non-empty directories are added to the XML index as "Category" elements, with the exception that
 * directories in the *root* of the GSA pipeline directory tree must have five-character names starting
 * with "P" in order to be indexed.
 *
 * The following are added to the XML index as "Resource" elements:
 *
 *     -bam files for Picard projects/samples listed in tsv files found within the GSA pipeline directory tree,
 *      provided the bams and their corresponding bais actually exist within the appropriate Picard directory.
 *      The most current version of each Picard bam is always chosen. Symlinks are created to the actual Picard
 *      bam/bai files in subdirectories named BAM_LINKS_SUBDIRECTORY within the directory in which each tsv file
 *      was found.
 *
 *     -"Project" vcf files -- that is, vcf files that have a tsv file of the same name present in the same directory.
 *
 *     -Files that end with one of the RESOURCE_FILE_EXTENSIONS
 *
 * Files/directories beginning with one of the prefixes in FILENAME_BLACKLIST are always excluded from indexing,
 * even if they fall into one of the categories above.
 *
 * All elements are sorted in the output XML file, with "Category" elements always coming before "Resource"
 * elements, and elements of the same type being sorted lexicographically by their "name" attributes.
 *
 * USAGE EXAMPLES:
 *
 * -Run with default settings:
 *
 *     java -classpath dist/GenomeAnalysisTK.jar org.broadinstitute.sting.pipeline.GSAPipelineIndexer
 *
 * This produces an index of the DEFAULT_ROOT_DIRECTORY, replacing the path to that directory with
 * DEFAULT_WEB_ROOT in all URLs generated, and writing the output to XML file DEFAULT_OUTPUT_FILE.
 *
 * -Run with custom settings:
 *
 *     java -classpath dist/GenomeAnalysisTK.jar org.broadinstitute.sting.pipeline.GSAPipelineIndexer \
 *          --rootDirectory directory_to_index \
 *          --webRoot web_alias_for_directory_to_index \
 *          --outputFile path_to/output.xml
 *
 * This produces an index of directory_to_index, replacing the path to that directory with
 * web_alias_for_directory_to_index in all URLs generated, and writing the output to XML file
 * path_to/output.xml
 *
 * @author David Roazen
 */
public class GSAPipelineIndexer extends CommandLineProgram {

    private static Logger logger = Logger.getLogger(GSAPipelineIndexer.class);

    @Argument(fullName = "rootDirectory", shortName = "r", doc = "Root of the directory tree to index", required = false)
    private File rootDirectory = new File(DEFAULT_ROOT_DIRECTORY);

    @Argument(fullName = "webRoot", shortName = "w", doc = "Alias for the root directory path in the URLs we generate", required = false)
    private String webRoot = DEFAULT_WEB_ROOT;

    @Output(fullName = "outputFile", shortName = "o", doc = "File to which the XML index should be written", required = false)
    private File outputFile = new File(DEFAULT_OUTPUT_FILE);

    @Argument(fullName = "useDebugSamples", shortName = "d", doc = "If set, index only samples in the debug directory (" + DEBUG_DIRECTORY_NAME + ")", required = false)
    private boolean useDebugSamples = false;

    @Argument(fullName = "printPaths", shortName = "p", doc = "If set, print the path to each TSV file as it's processed", required = false)
    private boolean printPaths = false;

    @Argument(fullName = "noSymlinks", shortName = "nsl", doc = "If set, do not create symlinks to Picard bams (for debugging/dry-run purposes only)", required = false)
    private boolean noSymlinks = false;

    @Argument(fullName = "debugSymlinkDir", shortName = "ds", doc = "If given, create all symlinks in the specified directory instead of in their " +
              "\"live\" locations within the directory tree being indexed (for debugging/dry-run purposes only)", required = false)
    private File debugSymlinkDir = null;

    @Argument(fullName = "noSort", shortName = "ns", doc = "If set, don't sort the XML elements", required = false)
    private boolean noSort = false;

    public static final String DEFAULT_ROOT_DIRECTORY = "/humgen/gsa-pipeline/";
    public static final String DEFAULT_WEB_ROOT = "gsa_pipeline_output";
    public static final String DEFAULT_OUTPUT_FILE = "index.xml";

    public static final String BAM_LINKS_SUBDIRECTORY = ".bamlinks";
    public static final String IGV_HYPERLINK = "https://iwww.broadinstitute.org/gsa/wiki/index.php/GSA_Firehose";
    public static final String HOST_PREFIX = "http://gsa-igv-web.broadinstitute.org/";
    public static final String PICARD_ROOT_DIRECTORY = "/seq/picard_aggregation/";
    public static final String DEBUG_DIRECTORY_NAME = "P0010";

    public static final String ROOT_ELEMENT_TAG_NAME = "Global";
    public static final String CATEGORY_ELEMENT_TAG_NAME = "Category";
    public static final String RESOURCE_ELEMENT_TAG_NAME = "Resource";
    public static final String NAME_ATTRIBUTE = "name";
    public static final String PATH_ATTRIBUTE = "path";
    public static final String HYPERLINK_ATTRIBUTE = "hyperlink";
    public static final String VERSION_ATTRIBUTE = "version";

    public static final String ROOT_ELEMENT_NAME_ATTRIBUTE = "GSA IGV Web";
    public static final String ROOT_ELEMENT_VERSION_ATTRIBUTE = "1.0";

    // Files ending in one of these extensions will be indexed as resources, provided they are not blacklisted:
    public static final String[] RESOURCE_FILE_EXTENSIONS = { ".bam",
                                                              ".maf_annotated.vcf",
                                                              ".cleaned.annotated.handfiltered.vcf",
                                                              ".filtered.annotated.vcf"
                                                            };

    // Files whose names start with one of these prefixes will be ignored during indexing:
    public static final String[] FILENAME_BLACKLIST = { BAM_LINKS_SUBDIRECTORY,
                                                        "queueScatterGather",
                                                        "temp-",
                                                        "Scatter",
                                                        "Intermediate"
                                                      };

    // Comparator used to sort Element nodes in the DOM tree as it's constructed. Category nodes
    // (representing directories on the filesystem) always sort before Resource nodes. Nodes
    // of the same type are sorted lexicographically by their "name" attribute (all nodes added
    // to the DOM tree are guaranteed to have a "name" attribute).
    private Comparator<Element> nodeComparator = new Comparator<Element>() {
        public int compare ( Element first, Element second ) {

            // Category nodes sort before Resource nodes:
            if ( isCategoryNode(first) && ! isCategoryNode(second) ) {
                return -1;
            }
            else if ( ! isCategoryNode(first) && isCategoryNode(second) ) {
                return 1;
            }

            // Nodes of the same type are sorted by their "name" attributes:
            return first.getAttribute(NAME_ATTRIBUTE).compareTo(second.getAttribute(NAME_ATTRIBUTE));
        }
    };

    protected int execute() throws Exception {
        Document xmlDocument = createXMLIndex();
        writeXMLFile(xmlDocument);

        return 0;
    }

    private Document createXMLIndex() {
        Document xmlDocument = newXmlDocument();
        Element rootNode = createXmlRootNode(xmlDocument);

        indexDirectoryTree(rootDirectory, xmlDocument, rootNode, 1);

        return xmlDocument;
    }

    private void writeXMLFile ( Document xmlDocument ) {
        try {
            TransformerFactory transformerFactory = TransformerFactory.newInstance();
            Transformer transformer = transformerFactory.newTransformer();

            DOMSource xml = new DOMSource(xmlDocument);
            StreamResult output = new StreamResult(outputFile);

            transformer.setOutputProperty(OutputKeys.INDENT, "yes");
            transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2");

            transformer.transform(xml, output);
        }
        catch ( TransformerConfigurationException e ) {
            throw new ReviewedStingException("Unable to initialize the XML output transformer: " + e.getMessage());
        }
        catch ( TransformerException e ) {
            throw new ReviewedStingException("An unrecoverable error occurred while writing the XML file: " + e.getMessage());
        }
    }

    private Document newXmlDocument() {
        DocumentBuilderFactory documentBuilderFactory = DocumentBuilderFactory.newInstance();

        DocumentBuilder documentBuilder = null;
        try {
            documentBuilder = documentBuilderFactory.newDocumentBuilder();
        }
        catch ( ParserConfigurationException e ) {
            throw new ReviewedStingException("Unable to initialize the XML Document Builder: " + e.getMessage());
        }

        return documentBuilder.newDocument();
    }

    private Element createXmlRootNode ( Document xmlDocument ) {
        Element rootNode = xmlDocument.createElement(ROOT_ELEMENT_TAG_NAME);

        rootNode.setAttribute(NAME_ATTRIBUTE, ROOT_ELEMENT_NAME_ATTRIBUTE);
        rootNode.setAttribute(HYPERLINK_ATTRIBUTE, IGV_HYPERLINK);
        rootNode.setAttribute(VERSION_ATTRIBUTE, ROOT_ELEMENT_VERSION_ATTRIBUTE);

        xmlDocument.appendChild(rootNode);

        return rootNode;
    }

    private void indexDirectoryTree ( File directory, Document xmlDocument, Element parentNode, int depth ) {
        File[] directoryContents = directory.listFiles();

        if ( directoryContents == null || directoryContents.length == 0 ) {
            return;
        }

        ArrayList<Element> childNodes = new ArrayList<Element>(directoryContents.length);

        for ( File currentFile : directoryContents ) {
            String currentFileName = currentFile.getName();
            String currentFilePath = currentFile.getPath();

            if ( skipFile(currentFileName, depth) ) {
                continue;
            }

            if ( currentFile.isDirectory() ) {
                Element directoryNode = newCategoryNode(xmlDocument, currentFileName);

                indexDirectoryTree(currentFile, xmlDocument, directoryNode, depth + 1);

                // Only add this node to its parent if it's non-empty (has child nodes) after the
                // recursive call. This prunes the DOM tree of empty elements as we construct it.
                if ( directoryNode.hasChildNodes() ) {
                    childNodes.add(directoryNode);
                }
            }
            else if ( isTSVFile(currentFileName) ) {
                if ( printPaths ) {
                    System.out.println("Processing " + currentFilePath);
                }

                processTSVFile(currentFile, xmlDocument, childNodes);
            }
            else if ( isResource(currentFilePath) ) {
                Element resourceNode = newResourceNode(xmlDocument, currentFileName, pathToUrl(currentFilePath));
                childNodes.add(resourceNode);
            }
        }

        if ( ! noSort ) {
            Collections.sort(childNodes, nodeComparator);
        }

        for ( Element child : childNodes ) {
            parentNode.appendChild(child);
        }
    }

    private Element newCategoryNode ( Document xmlDocument, String nameAttribute ) {
        Element newCategoryNode = xmlDocument.createElement(CATEGORY_ELEMENT_TAG_NAME);

        newCategoryNode.setAttribute(NAME_ATTRIBUTE, nameAttribute);

        return newCategoryNode;
    }

    private boolean isCategoryNode ( Element node ) {
        return node != null && node.getTagName().equals(CATEGORY_ELEMENT_TAG_NAME);
    }

    private Element newResourceNode ( Document xmlDocument, String nameAttribute, String pathAttribute ) {
        Element newResourceNode = xmlDocument.createElement(RESOURCE_ELEMENT_TAG_NAME);

        newResourceNode.setAttribute(NAME_ATTRIBUTE, nameAttribute);
        newResourceNode.setAttribute(PATH_ATTRIBUTE, pathAttribute);

        return newResourceNode;
    }

    private void processTSVFile ( File tsvFile, Document xmlDocument, ArrayList<Element> childNodes ) {
        try {
            for ( String tsvLine : new XReadLines(tsvFile, true) ) {
                String[] tsvLineTokens = tsvLine.split("\\t");

                if ( tsvLineTokens.length != 2 ) {
                    logger.warn("Malformed line in TSV file " + tsvFile.getAbsolutePath() + ": " + tsvLine);
                    continue;
                }

                String picardProject = tsvLineTokens[0].trim();
                String picardSample = IoUtil.makeFileNameSafe(tsvLineTokens[1].trim());

                String picardPathPrefix = String.format(PICARD_ROOT_DIRECTORY + "%s/%s/current/%s",
                                                       picardProject, picardSample, picardSample);
                File picardBam = new File(picardPathPrefix + ".bam");
                File picardBai = new File(picardPathPrefix + ".bai");

                if ( ! picardBam.isFile() ) {
                    logger.warn("File " + picardBam.getAbsolutePath() + " does not exist (listed in TSV file " + tsvFile.getAbsolutePath() + ")");
                }
                else if ( ! picardBai.isFile() ) {
                    logger.warn("Could not find a bai file for bam " + picardBam.getAbsolutePath());
                }
                else {
                    String bamUrl = createPicardBamSymlinks(picardBam, picardBai, debugSymlinkDir == null ? tsvFile.getParent() : debugSymlinkDir.getPath());
                    if ( bamUrl != null ) {
                        childNodes.add(newResourceNode(xmlDocument, picardBam.getName(), bamUrl));
                    }
                }
            }
        }
        catch ( FileNotFoundException e ) {
            logger.warn("Could not open TSV file " + tsvFile.getAbsolutePath() + " Reason: " + e.getMessage());
        }
    }

    private boolean skipFile ( String fileName, int depth ) {
        if ( useDebugSamples && depth == 1 ) {
            return ! fileName.equals(DEBUG_DIRECTORY_NAME);
        }

        return fileIsBlacklisted(fileName) ||
               (depth == 1 && ! fileName.matches("^P....$"));
    }

    private boolean fileIsBlacklisted ( String fileName ) {
        for ( String blacklistedPrefix : FILENAME_BLACKLIST ) {
            if ( fileName.startsWith(blacklistedPrefix) ) {
                return true;
            }
        }

        return false;
    }

    private boolean isTSVFile ( String fileName ) {
        return fileName.endsWith(".tsv");
    }

    private boolean isResource ( String filePath ) {
        return isProjectVCF(filePath) || hasResourceFileExtension(filePath);
    }

    private boolean isProjectVCF ( String filePath ) {
        File correspondingTSVFile = new File(filePath.replaceAll("\\.vcf$", ".tsv"));

        return filePath.endsWith(".vcf") && correspondingTSVFile.isFile();
    }

    private boolean hasResourceFileExtension ( String filePath ) {
        for ( String resourceFileExtension : RESOURCE_FILE_EXTENSIONS ) {
            if ( filePath.endsWith(resourceFileExtension) ) {
                return true;
            }
        }

        return false;
    }

    private String pathToUrl ( String path ) {
        String rootDirRemoved = path.replaceFirst(rootDirectory.getPath(), "");

        if ( rootDirRemoved.startsWith("/") ) {
            return HOST_PREFIX + webRoot + rootDirRemoved;
        }
        else {
            return HOST_PREFIX + webRoot + "/" + rootDirRemoved;
        }
    }

    private String createPicardBamSymlinks ( File picardBam, File picardBai, String parentDirectory ) {
        String linkDirectoryPath = parentDirectory + "/" + BAM_LINKS_SUBDIRECTORY;
        File linkDirectory = new File(linkDirectoryPath);
        String bamLinkPath = linkDirectoryPath + "/" + picardBam.getName();
        String baiLinkPath = linkDirectoryPath + "/" + picardBai.getName();

        if ( ! noSymlinks && ! linkDirectory.isDirectory() && ! linkDirectory.mkdir() ) {
            logger.warn("Could not create symlink directory " + linkDirectory);
            return null;
        }

        boolean bamLinkCreated = noSymlinks || createSymlink(picardBam.getAbsolutePath(), bamLinkPath);
        boolean baiLinkCreated = noSymlinks || createSymlink(picardBai.getAbsolutePath(), baiLinkPath);

        if ( bamLinkCreated && baiLinkCreated ) {
            return pathToUrl(bamLinkPath);
        }

        return null;
    }

    private boolean createSymlink ( String target, String link ) {
        try {

            // Don't look, my friends -- it's not pretty. Until we move to Java 7 our options for creating
            // symlinks in Java are a choice among evils, however.
            Runtime runtime = Runtime.getRuntime();
            Process lnProcess = runtime.exec("ln -fs " + target + " " + link);

            int lnReturnValue = lnProcess.waitFor();

            if ( lnReturnValue != 0 ) {
                logger.warn("Could not create symlink " + link + " to " + target + ": ln returned " + lnReturnValue);
                return false;
            }
        }
        catch ( InterruptedException e ) {
            logger.warn("Thread interrupted while creating symlink " + link + " to " + target + ": " + e.getMessage());
            return false;
        }
        catch ( IOException e ) {
            logger.warn("I/O error while creating symlink " + link + " to " + target + ": " + e.getMessage());
            return false;
        }

        return true;
    }

    public static void main ( String[] args ) {
        try {
            GSAPipelineIndexer instance = new GSAPipelineIndexer();
            start(instance, args);
            System.exit(CommandLineProgram.result);
        } catch ( UserException e ) {
            exitSystemWithUserError(e);
        } catch ( Exception e ) {
            exitSystemWithError(e);
        }
    }
}
