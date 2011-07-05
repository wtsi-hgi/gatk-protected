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

package org.broadinstitute.sting.gatk.walkers.diffengine;

import junit.framework.Test;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: 7/4/11
 * Time: 12:51 PM
 * A generic engine for comparing tree-structured objects
 */
public class DiffEngine {
    final protected static Logger logger = Logger.getLogger(DiffEngine.class);

    private final Map<String, DiffableReader> readers = new HashMap<String, DiffableReader>();
    private final int maxItems;

    public DiffEngine(int maxItems) {
        this.maxItems = maxItems;
        loadDiffableReaders();
    }

    // --------------------------------------------------------------------------------
    //
    // difference calculation
    //
    // --------------------------------------------------------------------------------

    public List<Difference> diff(DiffElement master, DiffElement test) {
        boolean masterIsNode = master instanceof DiffNode;
        boolean testIsNode = test instanceof DiffNode;

        if ( masterIsNode && testIsNode ) {
            return diff((DiffNode)master, (DiffNode)test);
        } else if ( ! masterIsNode && ! testIsNode ) {
            return diff((DiffLeaf) master, (DiffLeaf) test);
        } else {
            // structural difference in types.  one is node, other is leaf
            return Arrays.asList(new Difference(master, test));
        }
    }

    public List<Difference> diff(DiffNode master, DiffNode test) {
        Set<String> allNames = new HashSet<String>(master.getElementNames());
        allNames.addAll(test.getElementNames());
        List<Difference> diffs = new ArrayList<Difference>();

        for ( String name : allNames ) {
            DiffElement masterElt = master.getElement(name);
            DiffElement testElt = test.getElement(name);
            if ( masterElt == null && testElt == null ) {
                throw new ReviewedStingException("BUG: unexceptedly got two null elements for field: " + name);
            } else if ( masterElt == null || testElt == null ) { // if either is null, we are missing a value
                // todo -- should one of these be a special MISSING item?
                diffs.add(new Difference(masterElt, testElt));
            } else {
                diffs.addAll(diff(masterElt, testElt));
            }
        }

        return diffs;
    }

    public List<Difference> diff(DiffLeaf master, DiffLeaf test) {
        if ( master.getValue().equals(test.getValue()) ) {
            return Collections.emptyList();
        } else {
            return Arrays.asList(new Difference(master, test));
        }
    }

    // --------------------------------------------------------------------------------
    //
    // Reporting Structure
    //
    // --------------------------------------------------------------------------------

    public String reportItemizedDifference(List<Difference> diffs) {
        return "";
    }

    public String reportSummarizedDifferences(List<Difference> diff) {
        return "";
    }

    // --------------------------------------------------------------------------------
    //
    // plugin manager
    //
    // --------------------------------------------------------------------------------

    public void loadDiffableReaders() {
        List<Class<? extends DiffableReader>> drClasses = new PluginManager<DiffableReader>( DiffableReader.class ).getPlugins();

        logger.info("Loading diffable modules:");
        for (Class<? extends DiffableReader> drClass : drClasses ) {
            logger.info("\t" + drClass.getSimpleName());

            try {
                DiffableReader dr = drClass.newInstance();
                readers.put(dr.getName(), dr);
            } catch (InstantiationException e) {
                throw new ReviewedStingException("Unable to instantiate module '" + drClass.getSimpleName() + "'");
            } catch (IllegalAccessException e) {
                throw new ReviewedStingException("Illegal access error when trying to instantiate '" + drClass.getSimpleName() + "'");
            }
        }
    }

    protected Map<String, DiffableReader> getReaders() {
        return readers;
    }

    protected DiffableReader getReader(String name) {
        return readers.get(name);
    }

    public DiffableReader findReaderForFile(File file) {
        for ( DiffableReader reader : readers.values() )
            if (reader.canRead(file) )
                return reader;

        return null;
    }

    public DiffNode createDiffableFromFile(File file) {
        DiffableReader reader = findReaderForFile(file);
        if ( reader == null )
            throw new UserException("Unsupported file type: " + file);
        else
            return reader.readFromFile(file);
    }
}
