package org.geojson;

import java.util.List;

public class MultiLineString extends Geometry<List<LngLatAlt>> {
    private static final long serialVersionUID = 1L;

    public MultiLineString() {
        type = "MultiLineString";
    }

    public MultiLineString(List<LngLatAlt> line) {
        type = "MultiLineString";
        add(line);
    }
}