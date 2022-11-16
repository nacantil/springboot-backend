package org.geojson;

import java.util.Arrays;
import java.util.List;

import com.fasterxml.jackson.annotation.JsonIgnore;

public class Polygon extends Geometry<List<LngLatAlt>> {
    private static final long serialVersionUID = 1L;

    public Polygon() {
        type = "Polygon";
    }

    public Polygon(List<LngLatAlt> polygon) {
        type = "Polygon";
        add(polygon);
    }

    public Polygon(LngLatAlt... polygon) {
        type = "Polygon";
        add(Arrays.asList(polygon));
    }

    public void setExteriorRing(List<LngLatAlt> points) {
        coordinates.add(0, points);
    }

    @JsonIgnore
    public List<LngLatAlt> getExteriorRing() {
        assertExteriorRing();
        return coordinates.get(0);
    }

    @JsonIgnore
    public List<List<LngLatAlt>> getInteriorRings() {
        assertExteriorRing();
        return coordinates.subList(1, coordinates.size());
    }

    public List<LngLatAlt> getInteriorRing(int index) {
        assertExteriorRing();
        return coordinates.get(1 + index);
    }

    public void addInteriorRing(List<LngLatAlt> points) {
        assertExteriorRing();
        coordinates.add(points);
    }

    public void addInteriorRing(LngLatAlt... points) {
        assertExteriorRing();
        coordinates.add(Arrays.asList(points));
    }

    private void assertExteriorRing() {
        if (coordinates.isEmpty()) {
          throw new RuntimeException("No exterior ring definied");
        }
    }
}