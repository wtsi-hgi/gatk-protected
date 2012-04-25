/*
 * Copyright (c) 2012, The Broad Institute
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

package org.broadinstitute.sting.utils.bcf2;

import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

public class BCF2Reader {
    // TODO -- add BGZF support

    private InputStream bcf2InputStream;
    private BCF2Codec bcf2Codec = new BCF2Codec();
    private Object header;

    public BCF2Reader( File bcf2File ) {
        openBCF2File(bcf2File);
        readHeader();
    }

    public BCF2Reader( InputStream in ) {
        this.bcf2InputStream = in;
        readHeader();
    }

    public VariantContext readNextRecord() {
        // TODO: Tribble will eventually need to detect end-of-stream in its FeatureSource classes.
        //       Perhaps the codec should detect this for it...
        return bcf2Codec.decode(bcf2InputStream);
    }

    private void openBCF2File( File bcf2File ) {
        try {
            bcf2InputStream = new FileInputStream(bcf2File);
        }
        catch ( IOException e ) {
            throw new UserException.CouldNotReadInputFile(String.format("Failed to open BCF2 file %s for reading",
                                                          bcf2File.getAbsolutePath()));
        }
    }

    private void readHeader() {
        header = bcf2Codec.readHeader(bcf2InputStream);
    }
}
