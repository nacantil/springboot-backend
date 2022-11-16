package org.geojson;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class GeometryCollection extends GeoJsonObject implements Iterable<GeoJsonObject> {
    private static final long serialVersionUID = 1L;

    private List<GeoJsonObject> geometries = new ArrayList<GeoJsonObject>();

    public GeometryCollection() {
        type = "GeometryCollection";
    }

    public List<GeoJsonObject> getGeometries() {
        return geometries;
    }

    public void setGeometries(List<GeoJsonObject> geometries) {
        this.geometries = geometries;
    }

    @Override
    public Iterator<GeoJsonObject> iterator() {
        return geometries.iterator();
    }

    public GeometryCollection add(GeoJsonObject geometry) {
        geometries.add(geometry);
        return this;
    }
}
