package org.geojson;

public class LineString extends MultiPoint {
    private static final long serialVersionUID = 1L;

    public LineString() {
        type = "LineString";
    }

    public LineString(LngLatAlt... points) {
        super(points);
        type = "LineString";
    }
}