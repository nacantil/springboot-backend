package org.geojson;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;
import java.util.Objects;

public class Crs implements Serializable {
    private static final long serialVersionUID = 1L;

    private Map<String, Object> properties = new HashMap<String, Object>();
    private String type = "name";

    public Map<String, Object> getProperties() {
	return properties;
    }

    public String getType() {
	return type;
    }

    public void setProperties(Map<String, Object> properties) {
	this.properties = properties;
    }

    @Override
    public String toString() {
	return "Crs [type=" + type + ", properties=" + properties + "]";
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 97 * hash + Objects.hashCode(this.properties);
        hash = 97 * hash + Objects.hashCode(this.type);
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
        final Crs other = (Crs) obj;
        if (!Objects.equals(this.properties, other.properties)) {
            return false;
        }
        if (!Objects.equals(this.type, other.type)) {
            return false;
        }
        return true;
    }
}