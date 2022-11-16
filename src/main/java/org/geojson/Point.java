package org.geojson;

import java.util.Objects;


public class Point extends GeoJsonObject {
    private static final long serialVersionUID = 1L;

    private LngLatAlt coordinates;

    public Point() {
        type = "Point";
    }

    public Point(LngLatAlt coordinates) {
        type = "Point";
        this.coordinates = coordinates;
    }

    public Point(double longitude, double latitude) {
        type = "Point";
        coordinates = new LngLatAlt(longitude, latitude);
    }

    public Point(double longitude, double latitude, double altitude) {
        type = "Point";
        coordinates = new LngLatAlt(longitude, latitude, altitude);
    }

    public LngLatAlt getCoordinates() {
        return coordinates;
    }

    public void setCoordinates(LngLatAlt coordinates) {
        this.coordinates = coordinates;
    }

    @Override
    public String toString() {
        return "Point [coordinates=" + coordinates + "]";
    }

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 47 * hash + Objects.hashCode(this.coordinates);
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Point other = (Point) obj;
        if (!Objects.equals(this.coordinates, other.coordinates)) {
            return false;
        }
        return true;
    }

}