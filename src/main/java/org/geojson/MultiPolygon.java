package org.geojson;

import java.util.List;

public class MultiPolygon extends Geometry<List<List<LngLatAlt>>> {
    private static final long serialVersionUID = 1L;

    public MultiPolygon() {
        type = "MultiPolygon";
    }

    public MultiPolygon(Polygon polygon) {
        type = "MultiPolygon";
        add(polygon);
    }

    public MultiPolygon add(Polygon polygon) {
        coordinates.add(polygon.getCoordinates());
        return this;
    }
}