package org.geojson;

import java.io.Serializable;

import org.geojson.jackson.LngLatAltDeserializer;
import org.geojson.jackson.LngLatAltSerializer;

import com.fasterxml.jackson.databind.annotation.JsonDeserialize;
import com.fasterxml.jackson.databind.annotation.JsonSerialize;

@JsonDeserialize(using = LngLatAltDeserializer.class)
@JsonSerialize(using = LngLatAltSerializer.class)
public class LngLatAlt implements Serializable {
    private static final long serialVersionUID = 1L;

    private double longitude;
    private double latitude;
    private double altitude = Double.NaN;

    public LngLatAlt() {
    }

    public LngLatAlt(double longitude, double latitude) {
	this.longitude = longitude;
	this.latitude = latitude;
    }

    public LngLatAlt(double longitude, double latitude, double altitude) {
	this.longitude = longitude;
	this.latitude = latitude;
	this.altitude = altitude;
    }

    public boolean hasAltitude() {
	return !Double.isNaN(altitude);
    }

    public double getLongitude() {
	return longitude;
    }

    public void setLongitude(double longitude) {
	this.longitude = longitude;
    }

    public double getLatitude() {
	return latitude;
    }

    public void setLatitude(double latitude) {
	this.latitude = latitude;
    }

    public double getAltitude() {
	return altitude;
    }

    public void setAltitude(double altitude) {
	this.altitude = altitude;
    }

    @Override
    public String toString() {
	return "LngLatAlt [longitude=" + longitude + ", latitude=" + latitude + ", altitude=" + altitude + "]";
    }

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 79 * hash +
                (int) (Double.doubleToLongBits(this.longitude) ^
                (Double.doubleToLongBits(this.longitude) >>> 32));
        hash = 79 * hash +
                (int) (Double.doubleToLongBits(this.latitude) ^
                (Double.doubleToLongBits(this.latitude) >>> 32));
        hash = 79 * hash +
                (int) (Double.doubleToLongBits(this.altitude) ^
                (Double.doubleToLongBits(this.altitude) >>> 32));
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
        final LngLatAlt other = (LngLatAlt) obj;
        if (Double.doubleToLongBits(this.longitude) !=
                Double.doubleToLongBits(other.longitude)) {
            return false;
        }
        if (Double.doubleToLongBits(this.latitude) !=
                Double.doubleToLongBits(other.latitude)) {
            return false;
        }
        if (Double.doubleToLongBits(this.altitude) !=
                Double.doubleToLongBits(other.altitude)) {
            return false;
        }
        return true;
    }

}