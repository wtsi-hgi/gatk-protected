package org.broadinstitute.sting.utils.bcf2;

import org.broad.tribble.Feature;

import java.io.InputStream;

public interface GeneralizedFeatureCodec<T extends Feature> {

    public Feature decodeLoc( InputStream inputStream );

    public T decode( InputStream inputStream );

    public Class<T> getFeatureType();

    public Object readHeader( InputStream inputStream );
}
