package org.geojson;

import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonInclude.Include;
import com.fasterxml.jackson.annotation.JsonSubTypes;
import com.fasterxml.jackson.annotation.JsonSubTypes.Type;
import com.fasterxml.jackson.annotation.JsonTypeInfo;
import com.fasterxml.jackson.annotation.JsonTypeInfo.As;
import com.fasterxml.jackson.annotation.JsonTypeInfo.Id;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;

@JsonTypeInfo(property = "type", use = Id.NAME, include = As.EXISTING_PROPERTY)
@JsonSubTypes({
    @Type(GeometryCollection.class),
    @Type(Geometry.class),
    @Type(Feature.class),
    @Type(Polygon.class),
    @Type(MultiPolygon.class),
    @Type(FeatureCollection.class),
    @Type(Point.class),
    @Type(MultiPoint.class),
    @Type(MultiLineString.class),
    @Type(LineString.class)})
@JsonInclude(Include.NON_NULL)
public abstract class GeoJsonObject implements Serializable {
    private static final long serialVersionUID = 1L;

    private Crs crs;
    private double[] bbox;
    @JsonInclude(Include.NON_EMPTY)
    private Map<String, Object> properties = new TreeMap<>();
    protected String type;

    public String getType() {
        return type;
    }

    public void setType(String type) {
        this.type = type;
    }

    public Crs getCrs() {
        return crs;
    }

    public void setCrs(Crs crs) {
        this.crs = crs;
    }

    public double[] getBbox() {
        return bbox;
    }

    public void setBbox(double[] bbox) {
        this.bbox = bbox;
    }

    public void setProperty(String key, Object value) {
        properties.put(key, value);
    }

    @SuppressWarnings("unchecked")
    public <T> T getProperty(String key) {
        return (T) properties.get(key);
    }

    public Map<String, Object> getProperties() {
        return properties;
    }

    public void setProperties(Map<String, Object> properties) {
        this.properties = properties;
    }

    @Override
    public String toString() {
        return "GeoJsonObject [crs=" + crs + ", bbox=" + Arrays.toString(bbox)
                + ", properties=" + properties + "]";
    }

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 11 * hash + Objects.hashCode(this.crs);
        hash = 11 * hash + Arrays.hashCode(this.bbox);
        hash = 11 * hash + Objects.hashCode(this.properties);
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
        final GeoJsonObject other = (GeoJsonObject) obj;
        if (!Objects.equals(this.crs, other.crs)) {
            return false;
        }
        if (!Arrays.equals(this.bbox, other.bbox)) {
            return false;
        }
        if (!Objects.equals(this.properties, other.properties)) {
            return false;
        }
        return true;
    }

}
