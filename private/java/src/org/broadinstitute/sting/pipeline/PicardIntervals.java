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

public class PicardIntervals {
    private final String reference;
    private final String targets;

    public PicardIntervals(String reference, String targets) {
        if (reference == null)
            throw new IllegalArgumentException("reference is null");
        this.reference = reference;
        this.targets = targets;
    }

    public String getReference() {
        return reference;
    }

    public String getTargets() {
        return targets;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        PicardIntervals that = (PicardIntervals) o;

        return reference.equals(that.reference) && (targets == null ? that.targets == null : targets.equals(that.targets));
    }

    @Override
    public int hashCode() {
        int result = reference.hashCode();
        result = 31 * result + (targets != null ? targets.hashCode() : 0);
        return result;
    }

    @Override
    public String toString() {
        return String.format("PicardIntervals[reference='%s',targets='%s']", reference, targets);
    }
}
