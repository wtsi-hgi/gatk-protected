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

package org.broadinstitute.sting.gatk.walkers.annotator.genomicannotator;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * This is a container that holds all data corresponding to a single join table as specified by one -J arg (ex: -J bindingName1,/path/to/file,bindingName1.columnName=bindingName2.columnName2).
 * Some terminology:
 * 'bindingName' is an arbitrary label for a given table that is specified on the command line with either the -B or -J arg.
 * In the example above, bindingName1 is the 'local' binding name because it is attached to the join table file provided with this -J arg. bindingName2 is the 'external' binding name because
 * it corresponds to some other table specified previously with another -B or -J arg.
 *
 * The JoinTable object stores a map entry for each record in the join table. The entry's key is the value of the join column in a given record (eg. bindingName1.columnName in the above example),
 * and the entry value is an ArrayList representing the entire join table record.
 * The JoinTable object also stores some other join table parameters such as the column names that were parsed out of the file header, and the bindingNames and columnNames from the -J arg.
 *
 * The join operation is performed by looking up the value of the join column in the external table (the one that this table is being joined to), and then using this value to do a lookup
 * on the map - if there's a hit, it will provide the record from the join table that is to be joined with the record in the external table.
 *
 * More information can be found here: http://www.broadinstitute.org/gsa/wiki/index.php/GenomicAnnotator
 */
public class JoinTable
{
    //the list of join table column names parsed out of the file header.
    private List<String> columnNames; //not fully-qualified

    private String localBindingName;
    private String externalBindingName;
    private String externalColumnName;

    //stores a map entry for each record in the join table. The entry's key is the value of the join column in a given record (eg. bindingName.columnName in the above example),
    //and the entry value is an ArrayList representing the entire join table record.
    private HashMap<String, ArrayList<String>> joinColumnValueToRecords = new HashMap<String, ArrayList<String>>();

    private int maxSize;
    private boolean parsedFromFile = false;

    public JoinTable(int maxSize) {
        this.maxSize = maxSize;
    }

    /**
     * Parses the table from the given file using the JoinTableParser.
     *
     * @param filename The file containing the table.
     * @param localBindingName The binding name within the given file to join on.
     * @param localColumnName The column name within the given file to join on.
     * @param externalBindingName The binding name of another file (previously specified with either -B or -J).
     * @param externalColumnName The column name in this other file to join on.
     */
    public void parseFromFile(String filename, String localBindingName, String localColumnName, String externalBindingName, String externalColumnName)  {
        if(parsedFromFile) {
            throw new ReviewedStingException("parseFromFile(" + filename +", ..) called more than once");
        }
        parsedFromFile = true;

        setLocalBindingName(localBindingName);
        setExternalBindingName(externalBindingName);
        setExternalColumnName(externalColumnName);

        BufferedReader br = null;
        try
        {
            br = new BufferedReader(new FileReader(filename));
            final JoinTableParser parser = new JoinTableParser();

            //read in the header
            columnNames = parser.readHeader(br);

            //get the index of the localJoinColumnName
            int localColumnNameIdx = -1;
            for(int i = 0; i < columnNames.size(); i++) {
                final String columnName = columnNames.get(i);
                if(columnName.equals(localColumnName)) {
                    localColumnNameIdx = i;
                    break;
                }
            }

            if(localColumnNameIdx == -1) {
                throw new UserException.BadArgumentValue("-J", "The -J arg specifies an unknown column name: \"" + localColumnName + "\". It's not one of the column names in the header " + columnNames + " of the file: " + filename);
            }

            //read in all records and create a map entry for each
            String line;
            while((line = br.readLine()) != null) {
                final ArrayList<String> columnValues = parser.parseLine(line);
                if ( columnValues.size() < columnNames.size() )
                    throw new UserException.BadInput("the file: " + filename + " is malformed as there are not a sufficient number of columns for this line: " + line);
                final String joinColumnValue = columnValues.get(localColumnNameIdx);
                put(joinColumnValue, columnValues, filename);
            }
        }
        catch(IOException e)
        {
            throw new UserException.CouldNotReadInputFile(new File(filename), "Unable to parse file", e);
        }
        finally
        {
            try {
                if(br != null) {
                    br.close();
                }
            } catch(IOException e) {
                throw new ReviewedStingException("Unable to close file: " + filename, e);
            }
        }
    }

    /**
     * If the -J arg was:  -J bindingName1,/path/to/file,bindingName1.columnName=bindingName2.columnName2,
     * this returns bindingName1.
     * @return local binding name
     */
    public String getLocalBindingName() {
        return localBindingName;
    }

    public void setLocalBindingName(String localBindingName) {
        this.localBindingName = localBindingName;
    }

    /**
     * @return the list of join table column names parsed out of the file header.
     */
    public List<String> getColumnNames() {
        return columnNames; //not fully-qualified
    }

    protected void setColumnNames(List<String> columnNames) {
        this.columnNames = columnNames;
    }

    /**
     * If the -J arg was:  -J bindingName1,/path/to/file,bindingName1.columnName=bindingName2.columnName2,
     * this returns columnName2.
     * @return external column name
     */
    public String getExternalColumnName() {
        return externalColumnName;
    }

    protected void setExternalColumnName(
            String externalColumnName) {
        this.externalColumnName = externalColumnName;
    }

    /**
     * If the -J arg was:  -J bindingName1,/path/to/file,bindingName1.columnName=bindingName2.columnName2,
     * this returns bindingName2.
     * @return external binding name
     */
    public String getExternalBindingName() {
        return externalBindingName;
    }

    protected void setExternalBindingName(
            String externalBindingName) {
        this.externalBindingName = externalBindingName;
    }

    /**
     * Whether any join table records have the given value in the join column.
     * @param joinColumnValue value
     * @return true if the given name value exists in the file
     */
    public boolean containsJoinColumnValue(String joinColumnValue) {
        return joinColumnValueToRecords.containsKey(joinColumnValue);
    }

    /**
     * Returns all records in the table where the join column has the given value.
     * @param joinColumnValue column value
     * @return row
     */
    public ArrayList<String> get(String joinColumnValue) {
        return joinColumnValueToRecords.get(joinColumnValue);
    }

    /**
     * Adds the given record to the map.
     * @param joinColumnValue value
     * @param record row
     * @param filename the source file name
     */
    protected void put(String joinColumnValue, ArrayList<String> record, String filename) {
        if ( joinColumnValueToRecords.containsKey(joinColumnValue) )
            throw new UserException.BadInput("the file " + filename + " contains non-unique entries for the requested column, which isn't allowed.");
        joinColumnValueToRecords.put(joinColumnValue, record);
        if ( joinColumnValueToRecords.size() > maxSize )
            throw new UserException.BadInput("the file " + filename + " contains more than the maximum number (" + maxSize + ") of allowed rows (see the --maxJoinTableSize argument).");
    }
}
