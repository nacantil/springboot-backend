package org.geojson;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Objects;

public class FeatureCollection extends GeoJsonObject implements Iterable<Feature> {
    private static final long serialVersionUID = 1L;

    private List<Feature> features = new ArrayList<>();

    public FeatureCollection() {
        type = "FeatureCollection";
    }

    public List<Feature> getFeatures() {
        return features;
    }

    public void setFeatures(List<Feature> features) {
        this.features = features;
    }

    public FeatureCollection add(Feature feature) {
        features.add(feature);
        return this;
    }

    public void addAll(Collection<Feature> features) {
        this.features.addAll(features);
    }

    @Override
    public Iterator<Feature> iterator() {
        return features.iterator();
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 29 * hash + Objects.hashCode(this.features);
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
        final FeatureCollection other = (FeatureCollection) obj;
        if (!Objects.equals(this.features, other.features)) {
            return false;
        }
        return true;
    }
}