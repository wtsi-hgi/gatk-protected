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

package org.broadinstitute.sting.utils;

import java.util.*;

/**
 * NamedTable is a utility class for maintaining a table and accessing rows and columns
 * with named identifiers, rather than indicies that must be remembered.  It also grows
 * dynamically; you needn't specify the rows and columns before accessing them, so you can
 * continuously expand the table in situations where you don't necessarily know how many
 * rows or columns you'll have in the end.
 */
public class NamedTable {
    // If in the future, this class gets templatized, the counter variable should really become a CountedObject<T>.
    private HashMap<String, HashMap<String, Double>> table;
    private HashSet<String> rowNames;
    private HashSet<String> colNames;

    public NamedTable() {
        table = new HashMap<String, HashMap<String, Double>>();
        rowNames = new HashSet<String>();
        colNames = new HashSet<String>();
    }

    /**
     * Copy another table into this new table
     * @param ct  the table to copy
     */
    public NamedTable(NamedTable ct) {
        table = new HashMap<String, HashMap<String, Double>>();
        rowNames = new HashSet<String>();
        colNames = new HashSet<String>();

        for (String rowName : ct.getRowNames()) {
            for (String colName : ct.getColumnNames()) {
                this.set(rowName, colName, ct.get(rowName, colName));
            }
        }
    }

    /**
     * If the entry we're trying to access doesn't exist, create it.
     *
     * @param rowName  the name of the row
     * @param colName  the name of the column
     */
    private void verifyEntry(String rowName, String colName) {
        rowNames.add(rowName);
        colNames.add(colName);

        if (!table.containsKey(rowName)) {
            table.put(rowName, new HashMap<String, Double>());
        }

        if (!table.get(rowName).containsKey(colName)) {
            table.get(rowName).put(colName, 0.0);
        }
    }

    /**
     * Set an entry in the table
     *
     * @param rowName  the name of the row
     * @param colName  the name of the column
     * @param value    the value to set for the (row,column)-th entry
     */
    public void set(String rowName, String colName, double value) {
        verifyEntry(rowName, colName);

        table.get(rowName).put(colName, value);
    }

    /**
     * Get an entry in the table
     *
     * @param rowName  the name of the row
     * @param colName  the name of the column
     * @return the value of the (row,column)-th entry
     */
    public double get(String rowName, String colName) {
        verifyEntry(rowName, colName);

        return table.get(rowName).get(colName);
    }

    /**
     * For convenience, increment the (row,column)-th entry in the table so that the
     * user doesn't need to extract, increment, and then reassign the value.  One day,
     * this should probably be rewritten to use Andrey's CountedObject class.
     *
     * @param rowName  the name of the row
     * @param colName  the name of the column
     */
    public void increment(String rowName, String colName) {
        double value = get(rowName, colName);
        
        table.get(rowName).put(colName, value + 1.0);
    }

    /**
     * For convenience, decrement the (row,column)-th entry in the table so that the
     * user doesn't need to extract, increment, and then reassign the value.  One day,
     * this should probably be rewritten to use Andrey's CountedObject class.
     *
     * @param rowName  the name of the row
     * @param colName  the name of the column
     */
    public void decrement(String rowName, String colName) {
        double value = get(rowName, colName);

        table.get(rowName).put(colName, value - 1.0);
    }

    /**
     * Get a sorted list of all the rows in the table
     *
     * @return a sorted list of all the rows in the table
     */
    public ArrayList<String> getRowNames() {
        ArrayList<String> rows = new ArrayList<String>();
        Iterator<String> rowit = rowNames.iterator();
        while (rowit.hasNext()) {
            rows.add(rowit.next());
        }
        Collections.sort(rows);

        return rows;
    }

    /**
     * Get a sorted list of all the columns in the table
     *
     * @return a sorted list of all the columns in the table
     */
    public ArrayList<String> getColumnNames() {
        ArrayList<String> cols = new ArrayList<String>();
        Iterator<String> colit = colNames.iterator();
        while (colit.hasNext()) {
            cols.add(colit.next());
        }
        Collections.sort(cols);

        return cols;
    }

    /**
     * Get a new table representing a subset of the current table
     * @param rowNames  a list of rows to extract
     * @param colNames  a list of columns to extract
     * @return  the subsetted table
     */
    public NamedTable getSubset(ArrayList<String> rowNames, ArrayList<String> colNames) {
        NamedTable ct = new NamedTable();

        for (String rowName : rowNames) {
            for (String colName : colNames) {
                ct.set(rowName, colName, get(rowName, colName));
            }
        }

        return ct;
    }

    /* This stuff doesn't belong in this class, but I don't want
       to delete the code until it's moved somewhere appropriate */
    /*
    public boolean twoTailedFisherExactTest(double pValueLimit) {
        return (twoTailedFisherExactTest() < pValueLimit);
    }

    public double oneTailedFisherExactTestRight() {
        NamedTable ct = new NamedTable(this);

        double pCutoff = pValue();
        double pValue = pCutoff;

        while (ct.rotateRight()) {
            double pValuePiece = ct.pValue();

            if (pValuePiece <= pCutoff) {
                pValue += pValuePiece;
            }
        }

        return pValue;
    }

    public double oneTailedFisherExactTestLeft() {
        NamedTable ct = new NamedTable(this);

        double pCutoff = pValue();
        double pValue = pCutoff;

        while (ct.rotateLeft()) {
            double pValuePiece = ct.pValue();

            if (pValuePiece <= pCutoff) {
                pValue += pValuePiece;
            }
        }

        return pValue;
    }

    public double twoTailedFisherExactTest() {
        return oneTailedFisherExactTestLeft() + oneTailedFisherExactTestRight();
    }

    public double pValue() {
        double p = 0.0;
        
        if (rowNames.size() == 2 && colNames.size() == 2) {
            String[] rows = rowNames.toArray(new String[1]);
            String[] columns = colNames.toArray(new String[1]);

            double a = get(rows[0], columns[0]);
            double b = get(rows[0], columns[1]);
            double c = get(rows[1], columns[0]);
            double d = get(rows[1], columns[1]);
            double n = a + b + c + d;

            double p1 = Arithmetic.binomial(a + b, (long) a);
            double p2 = Arithmetic.binomial(c + d, (long) c);
            double pn = Arithmetic.binomial(n, (long) (a + c));

            p = p1*p2/pn;
        }

        return p;
    }

    public boolean rotateRight() {
        String[] rows = rowNames.toArray(new String[1]);
        String[] columns = colNames.toArray(new String[1]);

        decrement(rows[0], columns[0]);
        increment(rows[1], columns[0]);

        increment(rows[0], columns[1]);
        decrement(rows[1], columns[1]);

        return (get(rows[0], columns[0]) >= 0 && get(rows[1], columns[1]) >= 0);
    }

    public boolean rotateLeft() {
        String[] rows = rowNames.toArray(new String[1]);
        String[] columns = colNames.toArray(new String[1]);

        increment(rows[0], columns[0]);
        decrement(rows[1], columns[0]);

        decrement(rows[0], columns[1]);
        increment(rows[1], columns[1]);

        return (get(rows[0], columns[1]) >= 0 && get(rows[1], columns[0]) >= 0);
    }
    */

    /**
     * Get a nicely-formatted version of the contents of the table.
     *
     * @return  a String representing the contents of the table
     */
    public String toString() {
        String tableString = "";
        boolean headerPrinted = false;

        ArrayList<String> rows = getRowNames();
        ArrayList<String> cols = getColumnNames();

        for (String rowName : rows) {
            if (!headerPrinted) {
                tableString += "rowName  ";
                for (String colName : cols) {
                    tableString += "\t" + colName;
                }
                tableString += "\n";

                headerPrinted = true;
            }

            tableString += rowName;

            for (String colName : cols) {
                tableString += String.format("\t%7.7f", get(rowName, colName));
            }

            tableString += "\n";
        }

        return tableString;
    }
}
