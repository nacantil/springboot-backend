package org.geojson;

public class MultiPoint extends Geometry<LngLatAlt> {
    private static final long serialVersionUID = 1L;

    public MultiPoint() {
        type = "MultiPoint";
    }

    public MultiPoint(LngLatAlt... points) {
        super(points);
        type = "MultiPoint";
    }
}