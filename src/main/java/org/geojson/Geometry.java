package org.geojson;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

public abstract class Geometry<T> extends GeoJsonObject {

    protected List<T> coordinates = new ArrayList<T>();

    public Geometry() {
    }

    @SafeVarargs
    public Geometry(T... elements) {
	for (T coordinate : elements) {
	    coordinates.add(coordinate);
	}
    }

    public Geometry<T> add(T elements) {
	coordinates.add(elements);
	return this;
    }

    public List<T> getCoordinates() {
	return coordinates;
    }

    public void setCoordinates(List<T> coordinates) {
	this.coordinates = coordinates;
    }

    @Override
    public String toString() {
        return "Geometry{" + "coordinates=" + coordinates + '}';
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 79 * hash + Objects.hashCode(this.coordinates);
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
        final Geometry<?> other = (Geometry<?>) obj;
        if (!Objects.equals(this.coordinates, other.coordinates)) {
            return false;
        }
        return true;
    }
}