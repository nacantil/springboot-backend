package mil.navy.spawar.opt.common.utilities;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Date;
import java.util.List;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.JavaType;
import com.fasterxml.jackson.databind.ObjectMapper;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import mil.navy.spawar.opt.common.exception.OptException;
import mil.navy.spawar.opt.common.exception.ValidationException;

/**
 * Geometry Utilities
 * 
 */
public final class GeometryUtilities {
  public static final double EARTH_RADIUS_KILOMETERS = 6371.0;
  public static final double KILOMETERSPERNAUTICALMILE = 1.852;
  public static final double YARDSPERNAUTICALMILE = 2025.3718;
  private static final Logger LOGGER = LogManager.getLogger(GeometryUtilities.class);
  private static final com.esri.core.geometry.SpatialReference m_srsWGS84 =
      com.esri.core.geometry.SpatialReference.create(4326);
  private static final ObjectMapper objectMapper = new ObjectMapper();
  private static final String GEO_JSON_STRING_TEMPLATE = "{\"type\":{type},\"coordinates\":[[{coords}]]}";
  public static final double EARTH_RADIUS_IN_NAUTICAL_MILES = EARTH_RADIUS_KILOMETERS / KILOMETERSPERNAUTICALMILE; // Earth radius 6371 km
  private static final double DEGREESTORADIANS = Math.PI / 180.0;

  /**
   * private constructor - Make sure this class cannot be instantiated from outside the class.
   */
  private GeometryUtilities() {

  }

  /**
   * Convert an ESRI geometry to a geojson feature.
   *
   * @param geometry ESRI geometry to convert
   * @return the converted geojson feature
   */
  private static org.geojson.Feature convertToFeature(com.esri.core.geometry.Geometry geometry) {
    org.geojson.Feature feature = new org.geojson.Feature();
    try {
      org.geojson.GeoJsonObject geoJsonGeometry =
          objectMapper.readValue(GeometryUtilities.geometryToGeoJson(geometry), org.geojson.GeoJsonObject.class);
      feature = new org.geojson.Feature();
      feature.setGeometry(geoJsonGeometry);
    } catch (IOException e) {
      // If we get here, just return an empty feature
      LOGGER.warn("Could not convert ESRI geometry to geoJson feature", e);
    }
    return feature;
  }

  /**
   * Convert a geojson feature into an ESRI geometry.
   *
   * @param geojsonObject object to convert
   * @return the converted ESRI geometry
   */
  private static com.esri.core.geometry.Geometry convertToESRIGeometry(org.geojson.GeoJsonObject geojsonObject) {
    com.esri.core.geometry.Geometry esriGeometry = null;
    com.esri.core.geometry.MapGeometry esriMapGeometry;
    if (geojsonObject != null) {
      if (geojsonObject instanceof org.geojson.Feature) {
        geojsonObject = ((org.geojson.Feature) geojsonObject).getGeometry();
      }
      if (geojsonObject != null) {
        String geojson;
        try {
          geojson = objectMapper.writeValueAsString(geojsonObject);
          if ((geojsonObject instanceof org.geojson.Polygon) || (geojsonObject instanceof org.geojson.MultiPolygon)) {
            esriMapGeometry = com.esri.core.geometry.GeometryEngine.geoJsonToGeometry(geojson, 0,
                com.esri.core.geometry.Geometry.Type.Polygon);
            esriGeometry = esriMapGeometry.getGeometry();
          } else if (geojsonObject instanceof org.geojson.LineString) {
            esriMapGeometry = com.esri.core.geometry.GeometryEngine.geoJsonToGeometry(geojson, 0,
                com.esri.core.geometry.Geometry.Type.Polyline);
            esriGeometry = esriMapGeometry.getGeometry();
          } else if (geojsonObject instanceof org.geojson.MultiLineString) {
            esriMapGeometry = com.esri.core.geometry.GeometryEngine.geoJsonToGeometry(geojson, 0,
                com.esri.core.geometry.Geometry.Type.Polyline);
            esriGeometry = esriMapGeometry.getGeometry();
          } else if (geojsonObject instanceof org.geojson.Point) {
            esriMapGeometry = com.esri.core.geometry.GeometryEngine.geoJsonToGeometry(geojson, 0,
                com.esri.core.geometry.Geometry.Type.Point);
            esriGeometry = esriMapGeometry.getGeometry();
          } else if (geojsonObject instanceof org.geojson.MultiPoint) {
            esriMapGeometry = com.esri.core.geometry.GeometryEngine.geoJsonToGeometry(geojson, 0,
                com.esri.core.geometry.Geometry.Type.MultiPoint);
            esriGeometry = esriMapGeometry.getGeometry();
          }
        } catch (JsonProcessingException e) {
          // Error parsing geojson - just return a null
          LOGGER.warn("Could not parse geoJson string to ESRI geometry", e);
        }
      }
    }
    return esriGeometry;
  }
  
  /**
   * Compute the area of a polygon feature.  If feature does not have a polygon geometry, an exception is thrown
   * 
   * @param polygonFeature feature to compute the area of
   * @return the feature's area.  Since we are using WGS84 coordinates, the results are in degrees squared
   */
  public static double calculateArea(org.geojson.Feature polygonFeature) {
    return GeometryUtilities.calculateArea(polygonFeature, true);
  }
  
  /**
   * Compute the area of a polygon feature.  If feature does not have a polygon geometry, an exception is thrown
   * 
   * @param polygonFeature feature to compute the area of
   * @param checkForCircle true if a check for circle geometry is required
   * @return the feature's area.  Since we are using WGS84 coordinates, the results are in degrees squared
   */
  private static double calculateArea(org.geojson.Feature polygonFeature, boolean checkForCircle) {
    double area = 0.0;
    if (polygonFeature != null) {
      org.geojson.Feature newFeature = polygonFeature;
      if (checkForCircle) {
        // Check to see if the geometry is a circle
        newFeature = GeometryUtilities.checkForCircle(polygonFeature);
      }
      org.geojson.GeoJsonObject geometry = newFeature.getGeometry();
      if (geometry != null) {
        if ((geometry instanceof org.geojson.Polygon) || (geometry instanceof org.geojson.MultiPolygon)) {
          com.esri.core.geometry.Polygon esriPolygon = (com.esri.core.geometry.Polygon)convertToESRIGeometry(geometry);
          area = esriPolygon.calculateArea2D();
          // Positive Area for clockwise external rings, Negative Area for counterclockwise
          if (Math.abs(area) <= 0.0) {
            throw new ValidationException("Area of a Polygon should not be 0.0 - feature: " +
              polygonFeature.getProperty(Constants.OPT_SYS_NAME));
          }
        }
        else {
          throw new ValidationException("Area calculation is allowed only for Polygon or MultiPolygon - feature: "+
             polygonFeature.getProperty(Constants.OPT_SYS_NAME));
        }
      }
      else {
        throw new ValidationException("Geometry is null for feature: " +
            polygonFeature.getProperty(Constants.OPT_SYS_NAME));
      }
    }
    return area;
  }
  
  /**
   * Compute the 2D length of a line feature
   * 
   * @param lineFeature the feature to compute length of
   * @return the length of the feature
   */
  private static double calculateLength2D(org.geojson.Feature lineFeature) {
    double length = 0.0;
    if (lineFeature != null) {
      org.geojson.GeoJsonObject geometry = lineFeature.getGeometry();
      if (geometry != null) {
        if ((geometry instanceof org.geojson.LineString) || (geometry instanceof org.geojson.MultiLineString)) {
          com.esri.core.geometry.Polyline esriPolyline = (com.esri.core.geometry.Polyline)convertToESRIGeometry(geometry);
          length = esriPolyline.calculateLength2D();
          if (length <= 0.0) {
            throw new ValidationException("Length of a Line should not be 0.0 or less - feature: " +
              lineFeature.getProperty(Constants.OPT_SYS_NAME));
          }
        }
        else {
          throw new ValidationException("Length calculation is allowed only for LineString or MultiLineString - feature: "+
              lineFeature.getProperty(Constants.OPT_SYS_NAME));
        }
      }
      else {
        throw new ValidationException("Geometry is null for feature: " + lineFeature.getProperty(Constants.OPT_SYS_NAME));
      }
    }
    return length;
  }

  /**
   * Create a feature with a geometry of the envelope of the passed in geometry.
   *
   * @param geojsonObject object to compute the envelope for
   * @return the feature with the envelope as its geometry
   */
  public static org.geojson.Feature getEnvelope(org.geojson.Feature geojsonObject) {
    // Make sure a circle is converted to a polygon
    org.geojson.Feature newGeojsonObject = GeometryUtilities.checkForCircle(geojsonObject);
    org.geojson.Feature tmpGeojsonObject = GeometryUtilities.checkCrossesAntimeridian(newGeojsonObject);
    com.esri.core.geometry.Geometry geometry = convertToESRIGeometry(tmpGeojsonObject);
    com.esri.core.geometry.Envelope envelope = new com.esri.core.geometry.Envelope();
    geometry.queryEnvelope(envelope);
    org.geojson.LngLatAlt lowerLeft =
        new org.geojson.LngLatAlt(envelope.getLowerLeft().getX(), envelope.getLowerLeft().getY());
    org.geojson.LngLatAlt lowerRight =
        new org.geojson.LngLatAlt(envelope.getLowerRight().getX(), envelope.getLowerRight().getY());
    org.geojson.LngLatAlt upperLeft =
        new org.geojson.LngLatAlt(envelope.getUpperLeft().getX(), envelope.getUpperLeft().getY());
    org.geojson.LngLatAlt upperRight =
        new org.geojson.LngLatAlt(envelope.getUpperRight().getX(), envelope.getUpperRight().getY());
    org.geojson.Polygon geoJsonEnvelopeGeometry =
        new org.geojson.Polygon(lowerLeft, lowerRight, upperRight, upperLeft, lowerLeft);
    org.geojson.Feature geoJsonEnvelope = GeometryUtilities.createFeature();
    geoJsonEnvelope.setGeometry(geoJsonEnvelopeGeometry);
    geoJsonEnvelope.setProperty("lowerLeft", lowerLeft);
    geoJsonEnvelope.setProperty("lowerRight", lowerRight);
    geoJsonEnvelope.setProperty("upperLeft", upperLeft);
    geoJsonEnvelope.setProperty("upperRight", upperRight);
    return geoJsonEnvelope;
  }

  /**
   * Given a geojson object, figure out what type the geometry is.
   *
   * @param geojsonObject geojson object to determine the geometry type
   * @return the geometry type of the geojson object
   */
  public static String getGeometryType(org.geojson.GeoJsonObject geojsonObject) {
    String geometryType = "unknown";
    if (geojsonObject != null) {
      if (geojsonObject instanceof org.geojson.Feature) {
        geojsonObject = ((org.geojson.Feature) geojsonObject).getGeometry();
      }
      if (geojsonObject != null) {
        if (geojsonObject instanceof org.geojson.MultiPolygon) {
          geometryType = Constants.MULTIPOLYGON_TYPE;
        } else if (geojsonObject instanceof org.geojson.Polygon) {
          geometryType = Constants.POLYGON_TYPE;
        } else if (geojsonObject instanceof org.geojson.LineString) {
          geometryType = Constants.LINESTRING_TYPE;
        } else if (geojsonObject instanceof org.geojson.MultiLineString) {
          geometryType = Constants.MULTILINESTRING_TYPE;
        } else if (geojsonObject instanceof org.geojson.Point) {
          geometryType = Constants.POINT_TYPE;
        } else if (geojsonObject instanceof org.geojson.MultiPoint) {
          geometryType = Constants.MULTIPOINT_TYPE;
        }
      }
    }
    return geometryType;
  }

  /**
   * Intersect two features and return the resulting geometry.
   *
   * @param feature1 the first feature to intersect
   * @param feature2 the second feature to intersect
   * @return the returned feature that is the intersection of the input features
   */
  public static org.geojson.Feature intersect(org.geojson.Feature feature1, org.geojson.Feature feature2) {
    org.geojson.Feature intersectionFeature = new org.geojson.Feature();
    // If one or both of the geometries are circles, convert to polygons first
    org.geojson.Feature newFeature1 = GeometryUtilities.checkForCircle(feature1);
    org.geojson.Feature newFeature2 = GeometryUtilities.checkForCircle(feature2);
    if ((newFeature1 != null) && (newFeature2 != null)) {
      org.geojson.GeoJsonObject geoJsonGeometry1 = newFeature1.getGeometry();
      org.geojson.GeoJsonObject geoJsonGeometry2 = newFeature2.getGeometry();
      if (locationsCrossDateline(geoJsonGeometry1, geoJsonGeometry2)) {
        org.geojson.Feature tmpFeature1 = GeometryUtilities.makeLongitudesPositive(feature1);
        org.geojson.Feature tmpFeature2 = GeometryUtilities.makeLongitudesPositive(feature2);
        geoJsonGeometry1 = tmpFeature1.getGeometry();
        geoJsonGeometry2 = tmpFeature2.getGeometry();
      }
      com.esri.core.geometry.Geometry geometry1 = GeometryUtilities.convertToESRIGeometry(geoJsonGeometry1);
      com.esri.core.geometry.Geometry geometry2 = GeometryUtilities.convertToESRIGeometry(geoJsonGeometry2);
      com.esri.core.geometry.Geometry intersection =
          com.esri.core.geometry.GeometryEngine.intersect(geometry1, geometry2, m_srsWGS84);
      intersectionFeature = null;
      if ((intersection != null) && !intersection.isEmpty()) {
        intersectionFeature = GeometryUtilities.normalizeGeometry(GeometryUtilities.convertToFeature(intersection));
      }
      intersectionFeature = (intersectionFeature != null) ? intersectionFeature : new org.geojson.Feature();
    }
    return intersectionFeature;
  }

  /**
   * Union a list of geometries and return the resulting geometry.
   *
   * @param geometries the list of features containing the geometries to union
   * @return the feature containing the union of the geometries
   */
  public static org.geojson.Feature union(List<org.geojson.Feature> geometries) {
    org.geojson.Feature union = null;
    if ((geometries != null) && (geometries.size() > 0)) {
      List<org.geojson.Feature> newGeometries = GeometryUtilities.checkForCircles(geometries);
      List<com.esri.core.geometry.Geometry> checkedGeometries = new ArrayList<com.esri.core.geometry.Geometry>();
      boolean crossesDateline = false;
      for (org.geojson.Feature feature : newGeometries) {
        if (GeometryUtilities.isCrossDateLine(feature)) {
          crossesDateline = true;
          break;
        }
      }
      for (org.geojson.Feature feature : newGeometries) {
        org.geojson.Feature newFeature = (crossesDateline) ? GeometryUtilities.makeLongitudesPositive(feature) : feature;
        checkedGeometries.add(GeometryUtilities.convertToESRIGeometry(newFeature));
      }
      if (checkedGeometries.size() > 0) {
        com.esri.core.geometry.Geometry esriUnion =
            GeometryUtilities.normalizeGeometry(com.esri.core.geometry.GeometryEngine
                .union(checkedGeometries.toArray(new com.esri.core.geometry.Geometry[0]), m_srsWGS84));
        union = GeometryUtilities.convertToFeature(esriUnion);
        // When GeometryCollection items are unioned, a Crs value is inserted but not all geometries have one.
        // Hibernate seems to have a problem with using the default (no Crs) and having one specified (with the default)
        union.getGeometry().setCrs(null);
      }
    }
    return union;
  }

  /**
   * Buffer a geometry.
   *
   * @param feature the feature containing the geometry to buffer
   * @param bufferSizeNM the size of the buffer in nautical miles
   * @return the feature containing the buffered geometry
   */
  public static org.geojson.Feature bufferGeometry(org.geojson.Feature feature, double bufferSizeNM) {
    return GeometryUtilities.bufferGeometry(feature, bufferSizeNM, true);
  }

  /**
   * Buffer a geometry.
   *
   * @param feature the feature containing the geometry to buffer
   * @param bufferSizeNM the size of the buffer in nautical miles
   * @param chekcForCircle true if the circle geometry check is required
   * @return the feature containing the buffered geometry
   */
  private static org.geojson.Feature bufferGeometry(org.geojson.Feature feature, double bufferSizeNM, boolean checkForCircle) {
    org.geojson.Feature newFeature = feature;
    // Encountered some exceptions within this method (for whatever reasons), 
    // so I'll put a try-catch for safety reasons ...
    try {
      if ((feature != null) && (bufferSizeNM > 0.0)) {
        if (checkForCircle) {
          newFeature = GeometryUtilities.checkForCircle(feature);
        }
        newFeature = GeometryUtilities.checkCrossesAntimeridian(newFeature);
        if (newFeature != null) {
          com.esri.core.geometry.Geometry geometry = GeometryUtilities.convertToESRIGeometry(newFeature);
          if (geometry != null) {
            double bufferSizeRadians = bufferSizeNM * DEGREESTORADIANS;
            com.esri.core.geometry.Geometry bufferedGeometry =
                com.esri.core.geometry.GeometryEngine.buffer(geometry, m_srsWGS84, bufferSizeRadians);
            if (bufferedGeometry != null) {
              newFeature = GeometryUtilities.convertToFeature(bufferedGeometry);
              newFeature = GeometryUtilities.normalizeGeometry(newFeature);
            }
          }
        }
      }
     } catch (Exception e) {
	    System.out.println("EXCEPTION CAUGHT in the mil.navy.spawar.opt.common.utilities.GeometryUtilities.bufferGeometry method: " 
	       + e.getMessage());
	     newFeature = feature;
     }
    return newFeature;
  }

  /**
   * Buffer multiple geometries by multiple distances.
   *
   * @param features the array of features containing geometries to buffer
   * @param bufferSizesNM the array of buffer sizes in nautical miles to use for buffering
   * @return the array of features with buffered geometries
   */
  public static org.geojson.Feature[] bufferGeometries(org.geojson.Feature[] features, Double[] bufferSizesNM,
      boolean unionResult) {
    List<org.geojson.Feature> newBufferedFeatures = new ArrayList<org.geojson.Feature>();
    if ((features != null) && (bufferSizesNM != null) && (features.length == bufferSizesNM.length)) {
      List<org.geojson.Feature> newFeatures = GeometryUtilities.checkForCircles(Arrays.asList(features));
      List<Double> bufferSizesRadians = new ArrayList<Double>();
      List<com.esri.core.geometry.Geometry> geometries = new ArrayList<com.esri.core.geometry.Geometry>();
      for (int i = 0; i < newFeatures.size(); i++) {
        if (bufferSizesNM[i] > 0.0) {
          bufferSizesRadians.add(bufferSizesNM[i] * DEGREESTORADIANS);
          org.geojson.Feature newFeature = GeometryUtilities.checkCrossesAntimeridian(newFeatures.get(i));
          geometries.add(GeometryUtilities.convertToESRIGeometry(newFeature));
        }
      }
      if (geometries.size() > 0) {
        double[] bufferSizes = new double[geometries.size()];
        int ii = 0;
        for (Double bufferSizeRadians : bufferSizesRadians) {
          bufferSizes[ii++] = bufferSizeRadians;
        }
        com.esri.core.geometry.Geometry[] bufferedGeometries =
            com.esri.core.geometry.GeometryEngine.buffer(geometries.toArray(new com.esri.core.geometry.Geometry[0]),
                m_srsWGS84, bufferSizes, unionResult);
        if ((bufferedGeometries != null) && (bufferedGeometries.length > 0)) {
          for (int i = 0; i < bufferedGeometries.length; i++) {
            org.geojson.Feature newFeature = GeometryUtilities.convertToFeature(bufferedGeometries[i]);
            newFeature = GeometryUtilities.normalizeGeometry(newFeature);
            newBufferedFeatures.add(newFeature);
          }
        }
      }
    }
    return newBufferedFeatures.toArray(new org.geojson.Feature[0]);
  }
  
  /**
   * Build a geodesic buffer around a point
   * 
   * @param pointFeature center point of the buffered area
   * @param radiusNM radius of the buffer
   * @return the buffered feature
   */
  public static org.geojson.Feature bufferPointGeodesic(org.geojson.Feature pointFeature, double radiusNM) {
    org.geojson.Feature bufferedFeature = pointFeature;  // return the point if buffer size is zero
    if (radiusNM > 0.0) {
      bufferedFeature = new org.geojson.Feature();
      org.geojson.Point center = (org.geojson.Point)pointFeature.getGeometry();
      double radiusKM = radiusNM * GeometryUtilities.KILOMETERSPERNAUTICALMILE;
      double bearing = 0.0;
      double deltaBearing = 360.0 / 96;
      List<org.geojson.LngLatAlt> coordinates = new ArrayList<org.geojson.LngLatAlt>();
      while (bearing <= 360.0) {
        org.geojson.Point point = GeometryUtilities.pointFromRangeBearing(center, radiusKM, bearing);
        org.geojson.LngLatAlt coordinate = point.getCoordinates();
        coordinates.add(coordinate);
        bearing += deltaBearing;
      }
      org.geojson.Polygon polygon = new org.geojson.Polygon(coordinates);
      bufferedFeature.setGeometry(polygon);
    }
    return bufferedFeature;
  }
  
  /**
   * Buffer an array of geometries in a geodesic manner (instead of Euclidean)
   * 
   * @param features the array of features to buffer
   * @param bufferSizesNM the array of sizes to apply to the corresponding features
   * @param unionResult true if the results are to be unioned together
   * @return
   */
  public static org.geojson.Feature[] bufferGeometriesGeodesic(org.geojson.Feature[] features, Double[] bufferSizesNM,
      boolean unionResult) {
    int count = 0;
    List<org.geojson.Feature> returnedFeatures = new ArrayList<>();
    List<org.geojson.Feature> bufferedFeatures = new ArrayList<>();
    if ((features != null) && (bufferSizesNM != null) && (features.length == bufferSizesNM.length)) {
      List<org.geojson.Feature> newFeatures = GeometryUtilities.checkForCircles(Arrays.asList(features));
      for (org.geojson.Feature feature : newFeatures) {
        if (bufferSizesNM[count] > 0.0) {
          if (feature.getGeometry() instanceof org.geojson.Point) { // buffer a point in a geodesic manner
            bufferedFeatures.add(GeometryUtilities.bufferPointGeodesic(feature, bufferSizesNM[count]));
          }
          else { // don't have routines to buffer anything other than point geodesically.  Must use Euclidean buffer
            bufferedFeatures.add(GeometryUtilities.bufferGeometry(feature, bufferSizesNM[count], false));
          }
        }
        count++;
      }
    }
    if (bufferedFeatures.size() > 0) {
      if (unionResult) {
        org.geojson.Feature unionFeature = GeometryUtilities.union(bufferedFeatures);
        returnedFeatures.add(unionFeature);
      }
      else {
        returnedFeatures.addAll(bufferedFeatures);
      }
    }
    return returnedFeatures.toArray(new org.geojson.Feature[0]);
  }

  /**
   * Determine if two geometries are disjoint.
   *
   * @param feature1 a feature containing the first geometry to check
   * @param feature2 a feature containing the second geometry to check
   * @return true if the geometries are disjoint
   */
  public static boolean disjoint(org.geojson.Feature feature1, org.geojson.Feature feature2) {
    boolean isDisjoint = true;
    if ((feature1 != null) && (feature2 != null)) {
      org.geojson.Feature newFeature1 = GeometryUtilities.checkForCircle(feature1);
      org.geojson.Feature newFeature2 = GeometryUtilities.checkForCircle(feature2);
      org.geojson.GeoJsonObject geoJsonGeometry1 = newFeature1.getGeometry();
      org.geojson.GeoJsonObject geoJsonGeometry2 = newFeature2.getGeometry();
      if (locationsCrossDateline(geoJsonGeometry1, geoJsonGeometry2)) {
        org.geojson.Feature tmpFeature1 = GeometryUtilities.makeLongitudesPositive(newFeature1);
        org.geojson.Feature tmpFeature2 = GeometryUtilities.makeLongitudesPositive(newFeature2);
        geoJsonGeometry1 = tmpFeature1.getGeometry();
        geoJsonGeometry2 = tmpFeature2.getGeometry();
      }
      com.esri.core.geometry.Geometry geometry1 = GeometryUtilities.convertToESRIGeometry(geoJsonGeometry1);
      com.esri.core.geometry.Geometry geometry2 = GeometryUtilities.convertToESRIGeometry(geoJsonGeometry2);
      isDisjoint = com.esri.core.geometry.GeometryEngine.disjoint(geometry1, geometry2, m_srsWGS84);
    }
    return isDisjoint;
  }

  /**
   * Determine if a geometry is contained within another.
   *
   * @param container the feature with the geometry to check as the container
   * @param isContainedIn the feature with the geometry to check if it is contained in the container
   * @return true if isContained is contained in the container
   */
  public static boolean contains(org.geojson.Feature container, org.geojson.Feature isContainedIn) {
    boolean isContains = false;
    if ((container != null) && (container.getGeometry() != null) && (isContainedIn != null) && (isContainedIn.getGeometry() != null)) {
      org.geojson.Feature newContainer = GeometryUtilities.checkForCircle(container);
      org.geojson.Feature newIsContainedIn = GeometryUtilities.checkForCircle(isContainedIn);
      org.geojson.Feature newContainerFeature = GeometryUtilities.checkCrossesAntimeridian(newContainer);
      org.geojson.Feature newIsContainedFeature = newIsContainedIn;
      if (GeometryUtilities.isCrossDateLine(newContainer)) {
        // Since the container crosses, need to make sure the isContainedIn feature is handled properly
        newIsContainedFeature = GeometryUtilities.makeLongitudesPositive(newIsContainedIn);
      }
      com.esri.core.geometry.Geometry containerGeometry = GeometryUtilities.convertToESRIGeometry(newContainerFeature);
      com.esri.core.geometry.Geometry isContainedGeometry =
          GeometryUtilities.convertToESRIGeometry(newIsContainedFeature);
      isContains = com.esri.core.geometry.GeometryEngine.contains(containerGeometry, isContainedGeometry, m_srsWGS84);
    }
    return isContains;
  }

  /**
   * Determines whether two geometries are equal.
   *
   * @param feature1 the first feature to check
   * @param feature2 the second feature to check
   * @return true, if the two geometries are equal
   */
  public static boolean equals(org.geojson.Feature feature1, org.geojson.Feature feature2) {
    boolean equals = false;
    if ((feature1 != null) && (feature2 != null)) {
      org.geojson.Feature newFeature1 = GeometryUtilities.checkCrossesAntimeridian(GeometryUtilities.checkForCircle(feature1));
      org.geojson.Feature newFeature2 = GeometryUtilities.checkCrossesAntimeridian(GeometryUtilities.checkForCircle(feature2));
      com.esri.core.geometry.Geometry geometry1 = GeometryUtilities.convertToESRIGeometry(newFeature1);
      com.esri.core.geometry.Geometry geometry2 = GeometryUtilities.convertToESRIGeometry(newFeature2);
      equals = com.esri.core.geometry.GeometryEngine.equals(geometry1, geometry2, m_srsWGS84);
    }
    return equals;
  }
  
  /**
   * Simplify a feature geometry.
   * 
   * @param feature feature with geometry to simplify
   * @return a new feature object with a geometry that is simplified
   */
  public static org.geojson.Feature simplify(org.geojson.Feature feature) {
    org.geojson.Feature newFeature = GeometryUtilities.checkCrossesAntimeridian(GeometryUtilities.checkForCircle(feature));
    com.esri.core.geometry.Geometry geometry = GeometryUtilities.convertToESRIGeometry(newFeature);
    com.esri.core.geometry.Geometry simplifiedGeometry = com.esri.core.geometry.GeometryEngine.simplify(geometry, m_srsWGS84);
    org.geojson.Feature simplifiedFeature = GeometryUtilities.convertToFeature(simplifiedGeometry);
    return simplifiedFeature;
  }
  
  /**
   * Determine if a point is in a list of points
   * 
   * @param point the point to check
   * @param pointList the list of points to see if it contains the point to check
   * @return true if the point is in the list
   */
  public static boolean listContainsPoint(org.geojson.Point point, List<org.geojson.Point> pointList) {
    boolean isInList = false;
    org.geojson.Feature pointFeature = new org.geojson.Feature();
    pointFeature.setGeometry(point);
    org.geojson.Feature listPointFeature = new org.geojson.Feature();
    for (org.geojson.Point listPoint : pointList) {
      listPointFeature.setGeometry(listPoint);
      if (GeometryUtilities.equals(pointFeature, listPointFeature)) {
        isInList = true;
        break;
      }
    }
    return isInList;
  }

  /**
   * Get the distance to the nearest coordinate.
   *
   * @param geoJsonGeometry the geometry to check against
   * @param geoJsonPoint the point to determine the distance from
   * @return distance to geometry
   */
  public static double getDistanceToNearestCoordinate(org.geojson.GeoJsonObject geoJsonGeometry, org.geojson.Point geoJsonPoint) {
    org.geojson.Point nearestPoint = GeometryUtilities.getNearestCoordinate(geoJsonGeometry, geoJsonPoint);
    double distance = GeometryUtilities.rangeFromPoints(geoJsonPoint, nearestPoint);
    return distance;
  }

  /**
   * Get the distance to the nearest coordinate.
   *
   * @param Feature the feature whose geometry to check against
   * @param geoJsonPoint the point to determine the distance from
   * @return distance to geometry
   */
  public static double getDistanceToNearestCoordinate(org.geojson.Feature feature, org.geojson.Point geoJsonPoint) {
    double distance = Double.MAX_VALUE;
    if (feature != null) {
      org.geojson.Feature newFeature = GeometryUtilities.checkForCircle(feature);
      org.geojson.Point nearestPoint = GeometryUtilities.getNearestCoordinate(newFeature.getGeometry(), geoJsonPoint);
      distance = GeometryUtilities.rangeFromPoints(geoJsonPoint, nearestPoint);
    }
    return distance;
  }

  /**
   * Get the nearest coordinate.
   *
   * @param geoJsonGeometry the geometry to check against
   * @param geoJsonPoint the point to determine the nearest coordinate to
   * @return the nearest point
   */
  public static org.geojson.Point getNearestCoordinate(org.geojson.GeoJsonObject geoJsonGeometry, org.geojson.Point geoJsonPoint) {
    org.geojson.Point nearestPoint = null;
    com.esri.core.geometry.Proximity2DResult result = GeometryUtilities.getProximity2DResult(geoJsonGeometry, geoJsonPoint);
    com.esri.core.geometry.Point esriNearestPoint =
        (result != null) ? (com.esri.core.geometry.Point) result.getCoordinate() : null;
    if (esriNearestPoint != null) {
      org.geojson.Feature nearestFeaturePoint = GeometryUtilities.convertToFeature(esriNearestPoint);
      nearestFeaturePoint = GeometryUtilities.normalizeGeometry(nearestFeaturePoint);
      nearestPoint = (org.geojson.Point) nearestFeaturePoint.getGeometry();
    }
    return nearestPoint;
  }

  /**
   * Get the Proximity2DResult structure for determining the nearest coordinate.
   *
   * @param geoJsonGeometry geoJson object to determine the nearest point to
   * @param geoJsonPoint the point to find the nearest point from
   * @return the Proximity2DResult structure
   */
  private static com.esri.core.geometry.Proximity2DResult getProximity2DResult(
      org.geojson.GeoJsonObject geoJsonGeometry, org.geojson.Point geoJsonPoint) {
    com.esri.core.geometry.Proximity2DResult result = null;
    org.geojson.GeoJsonObject newGeoJsonGeometry = geoJsonGeometry;
    org.geojson.GeoJsonObject newGeoJsonPoint = geoJsonPoint;
    // if either feature or point crosses the anti-meridian, convert the longitudes to be positive for computational
    // consistency
    if (GeometryUtilities.locationsCrossDateline(geoJsonPoint, geoJsonGeometry)) {
      newGeoJsonGeometry = GeometryUtilities
          .createFeature(GeometryUtilities.makeDatelineCrossingLocationsPositive(extractLatLons(geoJsonGeometry)));
      newGeoJsonPoint = GeometryUtilities
          .createFeature(GeometryUtilities.makeDatelineCrossingLocationsPositive(extractLatLons(geoJsonPoint)));
    }
    com.esri.core.geometry.Geometry geometry = GeometryUtilities.convertToESRIGeometry(newGeoJsonGeometry);
    com.esri.core.geometry.Point point = (com.esri.core.geometry.Point) GeometryUtilities.convertToESRIGeometry(newGeoJsonPoint);
    result = com.esri.core.geometry.GeometryEngine.getNearestCoordinate(geometry, point, false);
    return result;
  }

  /**
   * Determine if a feature crosses the antimeridian.
   *
   * @param feature the feature to check
   * @return a new feature with geometry adjustments made if it crosses the antimeridian
   */
  private static org.geojson.Feature checkCrossesAntimeridian(org.geojson.Feature feature) {
    org.geojson.Feature newFeature = feature;
    if (feature != null) {
      // Assume that the feature has been checked for Circle geometry
      if (GeometryUtilities.isCrossDateLine(feature)) {
        newFeature = GeometryUtilities.makeLongitudesPositive(feature);
      }
    }
    return newFeature;
  }

  /**
   * Set the shape type - Circle, Relative, Oparea or a shape (Point, LineString, Polygon, etc.)
   * 
   * @param feature the feature to set the shapetype of
   */
  public static void setShapeType(org.geojson.Feature feature) {
    Map<String, Object> properties = feature.getProperties();
    boolean radiusType = properties.containsKey(Constants.OPT_SYS_RADIUS);
    boolean geomRelativeType = properties.containsKey(Constants.GEOJSON_GEOMRELATIVE);
    boolean opareasType = properties.containsKey(Constants.OPT_SYS_OPAREA);
    String shapeType = GeometryUtilities.getGeometryType(feature);
    shapeType = ("unknown".equals(shapeType)) ? null : shapeType;
    if (radiusType) {
      shapeType = Constants.CIRCLE_TYPE;
    } else if (geomRelativeType) {
      shapeType = Constants.RELATIVE_TYPE;
    } else if (opareasType) {
      shapeType = Constants.OPAREA_TYPE;
    }
    if (shapeType != null) {
      feature.setProperty(Constants.OPT_SHAPETYPE, shapeType);
    }
  }

  /**
   * Compute a point given a start point, speed, and bearing.
   *
   * @param startPoint starting point
   * @param speedKts speed in knots
   * @param bearing bearing in degrees
   * @param timeStart starting time
   * @param timeEnd ending time
   * @return the new point computed
   */
  public static org.geojson.Point pointFromSpeedBearing(org.geojson.Point startPoint, double speedKts, double bearing,
      Date timeStart, Date timeEnd) {
    double speedKPH = speedKts * GeometryUtilities.KILOMETERSPERNAUTICALMILE; // 1.852 kph per knot
    double distanceKM = speedKPH * ((timeEnd.getTime() - timeStart.getTime()) / 3600000.0); // 60*60*1000 = 3600000
                                                                                            // milliseconds per hour
    return GeometryUtilities.pointFromRangeBearing(startPoint, distanceKM, bearing);
  }

  /**
   * Compute a point given a start point, range, and bearing along a rhumb line.
   *
   * @param startPoint starting point
   * @param distanceKM distance in kilometers
   * @param direction bearing in degrees
   * @return the new point computed
   */
  public static org.geojson.Point pointFromRangeBearing(org.geojson.Point startPoint, double distanceKM,
      double direction) {
    double latStart = GeometryUtilities.toRadians(startPoint.getCoordinates().getLatitude());
    double lonStart = GeometryUtilities.toRadians(startPoint.getCoordinates().getLongitude());
    double distance = distanceKM / EARTH_RADIUS_KILOMETERS;
    double bearing = GeometryUtilities.toRadians(direction);
    // Calculation along the rhumb line See http://www.movable-type.co.uk/scripts/latlong.html and
    // look for the code under Rhumb lines - Distance
    double lat2 = latStart + distance * Math.cos(bearing);
    double delpsi = Math.log(Math.tan(lat2 / 2.0 + Math.PI / 4.0) / Math.tan(latStart / 2.0 + Math.PI / 4.0));
    double q = Math.abs(delpsi) > 10e-12 ? ((lat2 - latStart) / delpsi) : Math.cos(latStart); // E-W course becomes
                                                                                              // ill-conditioned with
                                                                                              // 0/0
    double lon2 = lonStart + distance * Math.sin(bearing) / q;
    // check for going past the pole, normalize latitude if so
    if (Math.abs(lat2) > Math.PI / 2) {
      lat2 = (lat2 > 0) ? (Math.PI - lat2) : (-Math.PI - lat2);
    }

    // Great circle distance See http://www.movable-type.co.uk/scripts/latlong.html and
    // look for the code under Destination point given distance and bearing from start point
    // var lat2 = Math.asin(Math.sin(latStart) * Math.cos(distance) + Math.cos(latStart) * Math.sin(distance) *
    // Math.cos(bearing));
    // var a = Math.atan2(Math.sin(bearing) * Math.sin(distance) * Math.cos(latStart),
    // Math.cos(distance) - Math.sin(latStart) * Math.sin(lat2));
    // var lon2 = lonStart + a;

    // normalize longitude to -180 to 180.
    lon2 = (lon2 + 3 * Math.PI) % (2 * Math.PI) - Math.PI;

    return new org.geojson.Point(GeometryUtilities.toDegrees(lon2), GeometryUtilities.toDegrees(lat2));
  }

  /**
   * Determine the distance between two points.
   *
   * @param point1 first point
   * @param point2 second point
   * @return distance computed
   */
  public static double rangeFromPoints(org.geojson.Point point1, org.geojson.Point point2) {
    double lat1 = GeometryUtilities.toRadians(point1.getCoordinates().getLatitude());
    double lon1 = GeometryUtilities.toRadians(point1.getCoordinates().getLongitude());
    double lat2 = GeometryUtilities.toRadians(point2.getCoordinates().getLatitude());
    double lon2 = GeometryUtilities.toRadians(point2.getCoordinates().getLongitude());
    double R = EARTH_RADIUS_IN_NAUTICAL_MILES;
    double dellat = lat2 - lat1;
    double delpsi = Math.log(Math.tan(Math.PI / 4.0 + lat2 / 2.0) / Math.tan(Math.PI / 4.0 + lat1 / 2.0));
    double q = Math.abs(delpsi) > 10e-12 ? ((dellat) / delpsi) : Math.cos(lat1); // E-W course becomes ill-conditioned
                                                                               // with 0/0
    // if dellon over 180 degrees, take shorter rhumb line across the anti-meridian:
    double dellon = lon2 - lon1;
    if (Math.abs(dellon) > Math.PI) {
      dellon = dellon > 0.0 ? -(2.0 * Math.PI - dellon) : (2.0 * Math.PI + dellon);
    }
    double dist = Math.sqrt(dellat * dellat + q * q * dellon * dellon) * R;
    return dist;
  }

  /**
   * Determine the bearing between two points along a rhumb line.
   *
   * @param point1 first point
   * @param point2 second point
   * @return bearing computed
   */
  public static double bearingFromPoints(org.geojson.Point point1, org.geojson.Point point2) {
    double lat1 = GeometryUtilities.toRadians(point1.getCoordinates().getLatitude());
    double lon1 = GeometryUtilities.toRadians(point1.getCoordinates().getLongitude());
    double lat2 = GeometryUtilities.toRadians(point2.getCoordinates().getLatitude());
    double lon2 = GeometryUtilities.toRadians(point2.getCoordinates().getLongitude());

    double delpsi = Math.log(Math.tan(Math.PI / 4.0 + lat2 / 2.0) / Math.tan(Math.PI / 4.0 + lat1 / 2.0));

    // if dLon over 180 degrees take shorter rhumb line across the anti-meridian:
    double dellon = lon2 - lon1;
    if (Math.abs(dellon) > Math.PI) {
      dellon = dellon > 0.0 ? -(2.0 * Math.PI - dellon) : (2.0 * Math.PI + dellon);
    }

    double brng = GeometryUtilities.toDegrees(Math.atan2(dellon, delpsi));
    return brng;
  }

  /**
   * Determine the initial bearing in radians between two points along a great circle.
   * 
   * @param point1 the first point.
   * @param point2 the second point.
   * @return the initial bearing in radians.
   */
  public static double initialBearingInRadiansFromPointsAlongGreatCircle(org.geojson.Point point1,
      org.geojson.Point point2) {
    double course;

    // Based on Ed Williams Aviation Formulary which assumes North latitudes and
    // West longitudes are positive and South latitudes and East longitudes negative.
    // Hence we change the signs of the longitudes in the calculations below.
    // https://edwilliams.org/avform.htm
    double lat1_radians = GeometryUtilities.toRadians(point1.getCoordinates().getLatitude());
    double lon1_radians = -GeometryUtilities.toRadians(point1.getCoordinates().getLongitude());
    double lat2_radians = GeometryUtilities.toRadians(point2.getCoordinates().getLatitude());
    double lon2_radians = -GeometryUtilities.toRadians(point2.getCoordinates().getLongitude());

    // Check if initial point is a pole
    if (Math.cos(lat1_radians) < 0.000001) {
      if (lat1_radians > 0) {
        course = Math.PI;
      } else {
        course = 2.0 * Math.PI;
      }
    } else {
      course = Math.atan2(Math.sin(lon1_radians - lon2_radians) * Math.cos(lat2_radians),
          Math.cos(lat1_radians) * Math.sin(lat2_radians)
              - Math.sin(lat1_radians) * Math.cos(lat2_radians) * Math.cos(lon1_radians - lon2_radians));
      // if (course < 0) {
      // course += 2.0 * Math.PI;
      // }
    }
    return course;
  }

  /**
   * Determine the initial bearing in degrees between two points along a great circle.
   * 
   * @param point1 the first point.
   * @param point2 the second point.
   * @return the initial bearing in degrees.
   */
  public static double initialBearingInDegreesFromPointsAlongGreatCircle(org.geojson.Point point1, org.geojson.Point point2) {
    return GeometryUtilities.toDegrees(GeometryUtilities.initialBearingInRadiansFromPointsAlongGreatCircle(point1, point2));
  }

  /**
   * Determines the angular cross track distance in radians from a point D to a great circle route defined by points A
   * and B.
   * 
   * @param pointA the start point of the great circle route; must not be a pole.
   * @param pointB the end point of the great circle route.
   * @param pointD the point to measure from.
   * @return the angular cross track distance in radians.
   */
  public static double angularCrossTrackDistanceInRadiansFromGreatCircle(org.geojson.Point pointA, org.geojson.Point pointB,
      org.geojson.Point pointD) {
    com.esri.core.geometry.Point esriPointA =
        (com.esri.core.geometry.Point) GeometryUtilities.convertToESRIGeometry(pointA);
    com.esri.core.geometry.Point esriPointD =
        (com.esri.core.geometry.Point) GeometryUtilities.convertToESRIGeometry(pointD);

    double distAD_radians = GeometryUtilities.calculateGreatCircleDistanceInRadians(esriPointA, esriPointD);
    double courseAD = GeometryUtilities.initialBearingInRadiansFromPointsAlongGreatCircle(pointA, pointD);
    double courseAB = GeometryUtilities.initialBearingInRadiansFromPointsAlongGreatCircle(pointA, pointB);
    return Math.asin(Math.sin(distAD_radians) * Math.sin(courseAD - courseAB));
  }

  /**
   * Determines the cross track distance in nautical miles from a point D to a great circle route defined by points A
   * and B.
   * 
   * @param pointA the start point of the great circle route; must not be a pole.
   * @param pointB the end point of the great circle route.
   * @param pointD the point to measure from.
   * @return the angular cross track distance in nautical miles.
   */
  public static double crossTrackDistanceInNauticalMilesFromGreatCircle(org.geojson.Point pointA, org.geojson.Point pointB,
      org.geojson.Point pointD) {
    double crossTrackDistance_nm =
        GeometryUtilities.angularCrossTrackDistanceInRadiansFromGreatCircle(pointA, pointB, pointD) * EARTH_RADIUS_IN_NAUTICAL_MILES;
    return crossTrackDistance_nm;
  }

  /**
   * Determines the angular along track distance in radians, the distance from point A along the great circle route
   * towards point B to the point closest to point D.
   * 
   * @param pointA the start point of the great circle route.
   * @param pointB the end point of the great circle route.
   * @param pointD the reference point.
   * @return the angular along track distance in radians.
   */
  public static double angularAlongTrackDistanceInRadiansForGreatCircle(org.geojson.Point pointA, org.geojson.Point pointB,
      org.geojson.Point pointD) {
    com.esri.core.geometry.Point esriPointA =
        (com.esri.core.geometry.Point) GeometryUtilities.convertToESRIGeometry(pointA);
    com.esri.core.geometry.Point esriPointD =
        (com.esri.core.geometry.Point) GeometryUtilities.convertToESRIGeometry(pointD);
    double distAD_radians = GeometryUtilities.calculateGreatCircleDistanceInRadians(esriPointA, esriPointD);
    double crossTrackDistance_radians = GeometryUtilities.angularCrossTrackDistanceInRadiansFromGreatCircle(pointA, pointB, pointD);
    double alongTrackDistance_radians = Math.acos(Math.cos(distAD_radians) / Math.cos(crossTrackDistance_radians));
    return alongTrackDistance_radians;
  }

  /**
   * Determines the along track distance in nautical miles, the distance from point A along the great circle route
   * towards point B to the point closest to point D.
   * 
   * @param pointA the start point of the great circle route.
   * @param pointB the end point of the great circle route.
   * @param pointD the reference point.
   * @return the along track distance in nautical miles.
   */
  public static double alongTrackDistanceInNauticalMilesForGreatCircle(org.geojson.Point pointA, org.geojson.Point pointB,
      org.geojson.Point pointD) {
    double alongTrackDistance_nm =
        GeometryUtilities.angularAlongTrackDistanceInRadiansForGreatCircle(pointA, pointB, pointD) * EARTH_RADIUS_IN_NAUTICAL_MILES;
    return alongTrackDistance_nm;
  }

  /**
   * Compute a point given a start point, angular distance in radians, and initial bearing in radians along a great
   * circle.
   * 
   * @param point1 the start point.
   * @param angularDistanceInRadians the angular range in radians.
   * @param initialBearingInRadians the initial bearing in radians.
   * @return the new point computed.
   */
  public static org.geojson.Point pointFromAngularRangeBearingAlongGreatCircle(org.geojson.Point point1,
      double angularDistanceInRadians, double initialBearingInRadians) {
    double lat1_radians = GeometryUtilities.toRadians(point1.getCoordinates().getLatitude());
    double lon1_radians = -GeometryUtilities.toRadians(point1.getCoordinates().getLongitude());
    double lat_radians = Math.asin(Math.sin(lat1_radians) * Math.cos(angularDistanceInRadians)
        + Math.cos(lat1_radians) * Math.sin(angularDistanceInRadians) * Math.cos(initialBearingInRadians));
    double dlon_radians =
        Math.atan2(Math.sin(initialBearingInRadians) * Math.sin(angularDistanceInRadians) * Math.cos(lat1_radians),
            Math.cos(angularDistanceInRadians) - Math.sin(lat1_radians) * Math.sin(lat_radians));
    double lon_radians = -(((lon1_radians - dlon_radians + Math.PI) % (2.0 * Math.PI)) - Math.PI);
    double lat_degrees = Math.toDegrees(lat_radians);
    double lon_degrees = Math.toDegrees(lon_radians);
    return new org.geojson.Point(lon_degrees, lat_degrees);
  }

  /**
   * Compute a point given a start point, range in nautical miles, and initial bearing in radians along a great circle.
   * 
   * @param point1 the start point.
   * @param rangeInNauticalMiles the range in nautical miles.
   * @param initialBearingInRadians the initial bearing in radians.
   * @return the new point computed.
   */
  public static org.geojson.Point pointFromRangeBearingAlongGreatCircle(org.geojson.Point point1,
      double rangeInNauticalMiles, double initialBearingInRadians) {
    return GeometryUtilities.pointFromAngularRangeBearingAlongGreatCircle(point1,
        rangeInNauticalMiles / EARTH_RADIUS_IN_NAUTICAL_MILES, initialBearingInRadians);
  }

  /**
   * Computes the point on the great circle route determined by point A and point B that is closest to point D.
   * 
   * @param pointA the start point of the great circle route.
   * @param pointB the end point of the great circle route.
   * @param pointD the reference point.
   * @return the closest point.
   */
  public static org.geojson.Point nearestPointOnGreatCircle(org.geojson.Point pointA, org.geojson.Point pointB,
      org.geojson.Point pointD) {
    double alongTrackDistance_radians = GeometryUtilities.angularAlongTrackDistanceInRadiansForGreatCircle(pointA, pointB, pointD);
    double initialBearingInRadians = GeometryUtilities.initialBearingInRadiansFromPointsAlongGreatCircle(pointA, pointB);
    org.geojson.Point nearest =
        GeometryUtilities.pointFromAngularRangeBearingAlongGreatCircle(pointA, alongTrackDistance_radians, initialBearingInRadians);
    return nearest;
  }

  /**
   * Gets an object property and converts it based on the clazz arguments.
   *
   * @param <T> the generic type
   * @param feature the feature
   * @param property the property
   * @param clazz the clazz
   * @return the object property
   */
  @SuppressWarnings("unchecked")
  public static <T> T getObjectProperty(org.geojson.Feature feature, String property, Class<?> clazz) {
    T propObject = feature.getProperty(property);
    if (propObject != null) {
      // If property object isn't correct class type, then convert it to correct type and save it back to the map.
      if (propObject.getClass().isInstance(new LinkedHashMap())) {
        propObject = (T) objectMapper.convertValue(propObject, clazz);
      } else if (propObject.getClass().isInstance(new ArrayList()) && ((ArrayList<?>) propObject).size() > 0 &&
          ((ArrayList<?>) propObject).get(0).getClass().isInstance(new LinkedHashMap())) {
        JavaType type = objectMapper.getTypeFactory().constructCollectionType(List.class, clazz);
        propObject = objectMapper.convertValue(propObject, type);
      }
    }
    return propObject;
  }

  /**
   * Convert degrees to radians.
   *
   * @param angleInDegrees degrees to convert
   * @return angle in radians
   */
  private static double toRadians(double angleInDegrees) {
    return Math.toRadians(angleInDegrees);
  }

  /**
   * Convert radians t degrees.
   *
   * @param angleInRadians angle in radians to convert
   * @return angle in degrees
   */
  private static double toDegrees(double angleInRadians) {
    return Math.toDegrees(angleInRadians);
  }

  /**
   * Convert between -180 to 180 to 0 to 360.
   *
   * @param feature the feature containing the geometry to convert
   * @return a feature with the geometry between 0 and 360
   */
  private static org.geojson.Feature makeLongitudesPositive(org.geojson.Feature feature) {
    // Assume that feature has already been checked for Circle geometry
    org.geojson.Feature newFeature = feature;
    if (feature != null) {
      newFeature = GeometryUtilities.createFeature();
      newFeature.setBbox(feature.getBbox());
      newFeature.setCrs(feature.getCrs());
      newFeature.setId(feature.getId());
      newFeature.setProperties(feature.getProperties());
      if (feature.getGeometry() instanceof org.geojson.MultiPolygon) {
        org.geojson.MultiPolygon multiPolygon = (org.geojson.MultiPolygon) feature.getGeometry();
        List<List<List<org.geojson.LngLatAlt>>> lngLatAltList = multiPolygon.getCoordinates();
        org.geojson.MultiPolygon newMultiPolygon = new org.geojson.MultiPolygon();
        for (List<List<org.geojson.LngLatAlt>> lngLatAlts : lngLatAltList) {
          org.geojson.Polygon polygon = new org.geojson.Polygon();
          for (List<org.geojson.LngLatAlt> lngLatAlt : lngLatAlts) {
            List<List<org.geojson.LngLatAlt>> tmpLngLatAlts = new ArrayList<>();
            org.geojson.Polygon tmpPolygon = new org.geojson.Polygon();
            tmpLngLatAlts.add(lngLatAlt);
            tmpPolygon.setCoordinates(tmpLngLatAlts);
            List<gov.nasa.worldwind.geom.LatLon> latLons =
                GeometryUtilities.makeDatelineCrossingLocationsPositive(GeometryUtilities.extractLatLons(tmpPolygon));
            org.geojson.Polygon newPolygon =
                (org.geojson.Polygon) ((org.geojson.Feature) GeometryUtilities.createFeature(latLons)).getGeometry();
            polygon.add(newPolygon.getCoordinates().get(0));
          }
          newMultiPolygon.add(polygon);
        }
        newFeature.setGeometry(newMultiPolygon);
      } else {
        newFeature.setGeometry((org.geojson.GeoJsonObject) GeometryUtilities.createFeature(
           GeometryUtilities.makeDatelineCrossingLocationsPositive(GeometryUtilities.extractLatLons(feature))).getGeometry());
      }
    }
    return newFeature;
  }

  /**
   * Check to see if the geometry crosses the antimeridian.
   *
   * @param feature the feature containing the geometry to check
   * @return true if the geometry crosses the antimeridian
   */
  private static boolean isCrossDateLine(org.geojson.Feature feature) {
    boolean crossesDateLine = false;
    if (feature != null) {
      // Assume the feature has already been converted from a circle
      crossesDateLine = GeometryUtilities.isCrossDateLine(feature.getGeometry());
    }
    return crossesDateLine;
  }


  /**
   * Check to see if the geometry crosses the antimeridian.
   *
   * @param geoJsonObject the geometry to check
   * @return true if the geometry crosses the antimeridian
   */
  private static boolean isCrossDateLine(org.geojson.GeoJsonObject geoJsonObject) {
    boolean crosses = false;
    if (geoJsonObject != null) {
      crosses = GeometryUtilities.isCrossDateLine(GeometryUtilities.extractLatLons(geoJsonObject, false));
    }
    return crosses;
  }

  /**
   * Check to see if the LatLon list crosses the antimeridian.
   *
   * @param latLons the lat lons
   * @return true if the geometry crosses the antimeridian
   */
  private static boolean isCrossDateLine(List<gov.nasa.worldwind.geom.LatLon> latLons) {
    return gov.nasa.worldwind.geom.LatLon.locationsCrossDateLine(latLons);
  }

  /**
   * Make the list of LatLon locations positive by transforming between 0 and 360 degrees.
   *
   * @param locations the list of LatLons to transform
   * @return a list of transformed LatLons
   */
  private static List<gov.nasa.worldwind.geom.LatLon> makeDatelineCrossingLocationsPositive(
      Iterable<? extends gov.nasa.worldwind.geom.LatLon> locations) {
    return gov.nasa.worldwind.geom.LatLon.makeDatelineCrossingLocationsPositive(locations);
  }

  /**
   * Create a geojson feature.
   *
   * @return the created feature
   */
  private static org.geojson.Feature createFeature() {
    org.geojson.Feature feature = new org.geojson.Feature();
    // feature.setCrs(DEFAULT_CRS);
    // feature.setProperty("name", "EPSG:4326");
    return feature;
  }

  /**
   * Create a geojson feature with geometry from a list of LatLons.
   *
   * @param latLons the list of LatLons to use for creating the geometry
   * @return the created feature
   */
  private static org.geojson.Feature createFeature(List<gov.nasa.worldwind.geom.LatLon> latLons) {
    org.geojson.Feature feature = createFeature();
    org.geojson.GeoJsonObject geoJsonObject = GeometryUtilities.constructGeoJsonObject(latLons);
    feature.setGeometry(geoJsonObject);
    return feature;
  }

  /**
   * Get the list of LatLons from the geojson object geometry.
   *
   * @param geoJsonObject geojson object containing the geometry to get the LatLon list from
   * @return the list of LatLns extracted
   */
  public static List<gov.nasa.worldwind.geom.LatLon> extractLatLons(org.geojson.GeoJsonObject geoJsonObject) {
    return GeometryUtilities.extractLatLons(geoJsonObject, true);
  }

  /**
   * Get the list of LatLons from the geojson object geometry.
   *
   * @param geoJsonObject geojson object containing the geometry to get the LatLon list from
   * @param checkForCircle true if a feature should be checked for circle geometry
   * @return the list of LatLns extracted
   */
  public static List<gov.nasa.worldwind.geom.LatLon> extractLatLons(org.geojson.GeoJsonObject geoJsonObject, boolean checkForCircle) {
    List<gov.nasa.worldwind.geom.LatLon> latLons = new ArrayList<gov.nasa.worldwind.geom.LatLon>();
    if (geoJsonObject instanceof org.geojson.Feature) {
      if (checkForCircle) {
        geoJsonObject = GeometryUtilities.checkForCircle((org.geojson.Feature) geoJsonObject);
      }
      geoJsonObject = ((org.geojson.Feature) geoJsonObject).getGeometry();
    }
    if (geoJsonObject instanceof org.geojson.Point) {
      org.geojson.Point point = (org.geojson.Point) geoJsonObject;
      org.geojson.LngLatAlt coordinates = point.getCoordinates();
      gov.nasa.worldwind.geom.LatLon latlon = gov.nasa.worldwind.geom.LatLon.fromDegrees(coordinates.getLatitude(),
          coordinates.getLongitude());
      latLons.add(latlon);
    } else if (geoJsonObject instanceof org.geojson.LineString) {
      org.geojson.LineString line = (org.geojson.LineString) geoJsonObject;
      List<org.geojson.LngLatAlt> coordinates = line.getCoordinates();
      latLons.addAll(coordinates.stream()
          .map(lla -> gov.nasa.worldwind.geom.LatLon.fromDegrees(lla.getLatitude(), lla.getLongitude()))
          .collect(Collectors.toList()));
    } else if (geoJsonObject instanceof org.geojson.MultiLineString) {
      org.geojson.MultiLineString multiLine = (org.geojson.MultiLineString) geoJsonObject;
      latLons.addAll(GeometryUtilities.collectLatLons(multiLine.getCoordinates()));
    } else if (geoJsonObject instanceof org.geojson.Polygon) {
      org.geojson.Polygon polygon = (org.geojson.Polygon) geoJsonObject;
      latLons.addAll(GeometryUtilities.collectLatLons(polygon.getCoordinates()));
    } else if (geoJsonObject instanceof org.geojson.MultiPolygon) {
      org.geojson.MultiPolygon multiPolygon = (org.geojson.MultiPolygon) geoJsonObject;
      List<List<List<org.geojson.LngLatAlt>>> coordinates = multiPolygon.getCoordinates();
      for (List<List<org.geojson.LngLatAlt>> polygonCoords : coordinates) {
        latLons.addAll(GeometryUtilities.collectLatLons(polygonCoords));
      }
    }
     return latLons;
    }

  public static List<gov.nasa.worldwind.geom.LatLon> collectLatLons(List<List<org.geojson.LngLatAlt>> coordinates) {
    return 
        coordinates
        .stream()
        .flatMap(Collection::stream)
        .map(lla -> gov.nasa.worldwind.geom.LatLon.fromDegrees(lla.getLatitude(), lla.getLongitude()))
        .collect(Collectors.toList());
  }

  /**
   * Get the list of LatLons from the geojson object geometry.
   *
   * @param multiPolygon MultiPolygon object containing the geometry to get the LatLon list from
   * @return the list of LatLns extracted
   */
  private static List<gov.nasa.worldwind.geom.LatLon> extractLatLons(org.geojson.MultiPolygon multiPolygon) {
    List<gov.nasa.worldwind.geom.LatLon> latLons = new ArrayList<gov.nasa.worldwind.geom.LatLon>();
    if (multiPolygon != null) {
      List<List<List<org.geojson.LngLatAlt>>> coordinatesList = multiPolygon.getCoordinates();
      for (List<List<org.geojson.LngLatAlt>> coordinates : coordinatesList) {
        org.geojson.Polygon polygon = new org.geojson.Polygon();
        polygon.setCoordinates(coordinates);
        latLons.addAll(GeometryUtilities.extractLatLons(polygon));
      }
    }
    return latLons;
  }

  /**
   * Build a geojson object from a LatLon object.
   *
   * @param latLon the LatLon object used to create the geojson object
   * @return the created geojson object
   */
  private static org.geojson.GeoJsonObject constructGeoJsonObject(gov.nasa.worldwind.geom.LatLon latLon) {
    org.geojson.GeoJsonObject geoJsonObject = null;
    String geoJsonString = GEO_JSON_STRING_TEMPLATE;
    geoJsonString = geoJsonString.replace("{type}", "\"Point\"");
    StringBuilder coordinates = new StringBuilder();
    coordinates.append("[").append(latLon.longitude.degrees).append(",").append(latLon.latitude.degrees).append("]");
    geoJsonString = geoJsonString.replace("[[{coords}]]", coordinates.toString());
    geoJsonObject = GeometryUtilities.createGeoJsonObject(geoJsonString);
    return geoJsonObject;
  }

  /**
   * Construct a geojson object from a list of LatLons.
   *
   * @param latLons the list of LatLons to use to construct the geojson object
   * @return the constructed geojson object
   */
  private static org.geojson.GeoJsonObject constructGeoJsonObject(List<gov.nasa.worldwind.geom.LatLon> latLons) {
    org.geojson.GeoJsonObject geoJsonObject = null;
    if (latLons != null && latLons.size() > 0) {
      boolean lineStringType = false;
      String geoJsonString = GEO_JSON_STRING_TEMPLATE;
      int size = latLons.size();
      if (size == 1) {
        geoJsonString = geoJsonString.replace("{type}", "\"Point\"");
      } else {
        gov.nasa.worldwind.geom.LatLon firstLatLon = latLons.get(0);
        gov.nasa.worldwind.geom.LatLon lastLatLon = latLons.get(size - 1);
        if (gov.nasa.worldwind.geom.LatLon.equals(firstLatLon, lastLatLon)) {
          geoJsonString = geoJsonString.replace("{type}", "\"Polygon\"");
        } else {
          lineStringType = true;
          geoJsonString = geoJsonString.replace("{type}", "\"LineString\"");
        }
      }
      StringBuilder coordinates = new StringBuilder();
      for (int i = 0; i < size; i++) {
        gov.nasa.worldwind.geom.LatLon latLon = latLons.get(i);
        coordinates.append("[").append(latLon.longitude.degrees).append(",").append(latLon.latitude.degrees)
            .append("]");
        if (i < (size - 1)) {
          coordinates.append(",");
        }
      }
      if (size == 1) {
        geoJsonString = geoJsonString.replace("[[{coords}]]", coordinates.toString());
      } else if (lineStringType) {
        geoJsonString = geoJsonString.replace("[{coords}]", coordinates.toString());
      } else {
        geoJsonString = geoJsonString.replace("{coords}", coordinates.toString());
      }
      geoJsonObject = GeometryUtilities.createGeoJsonObject(geoJsonString);
    }
    return geoJsonObject;
  }

  /**
   * Create a geojson object from a string.
   *
   * @param geoJsonString a string in geojson format
   * @return the constructed geojson object
   * @throws IllegalArgumentException the illegal argument exception
   */
  private static org.geojson.GeoJsonObject createGeoJsonObject(String geoJsonString) throws IllegalArgumentException {
    org.geojson.GeoJsonObject geoJsonObject = null;
    try {
      geoJsonObject = objectMapper.readValue(geoJsonString, org.geojson.GeoJsonObject.class);
    } catch (IOException e) {
      throw new IllegalArgumentException(e);
    }
    return geoJsonObject;
  }

  /**
   * Convert an ESRI geometry to a geojson string.
   *
   * @param p_geom the ESRI geometry to convert
   * @return the geojson string
   */
  private static String geometryToGeoJson(com.esri.core.geometry.Geometry p_geom) {
    String geoJson = com.esri.core.geometry.GeometryEngine.geometryToGeoJson(m_srsWGS84, p_geom);
    if (p_geom instanceof com.esri.core.geometry.Polygon) {
      try {
        // GeometryEngine.geometryToGeoJson is self-closing in that it always adds the first point as the last point
        // whether or not the initial last point already matches the first point. This results in multiple points at the
        // same location.  When there are two there isn't a problem on the javascript client, but if there are 3 or
        // more, they are turned into "NaN" (not a number) and this plays havoc with the client.
        //
        // Need to check if the last two points of the exterior ring are the same. If so, remove the last point and get
        // the GeoJson string again
        org.geojson.GeoJsonObject geometry = (org.geojson.GeoJsonObject) GeometryUtilities.createGeoJsonObject(geoJson);
        org.geojson.GeoJsonObject newGeometry = geometry;
        if (geometry instanceof org.geojson.Polygon) {
          org.geojson.Polygon polygon = (org.geojson.Polygon) geometry;
          newGeometry = GeometryUtilities.removeLastPoint(polygon);
        } else if (geometry instanceof org.geojson.MultiPolygon) {
          List<List<List<org.geojson.LngLatAlt>>> multiPolygonCoordinates =
              ((org.geojson.MultiPolygon) geometry).getCoordinates();
          org.geojson.MultiPolygon newMultiPolygon = new org.geojson.MultiPolygon();
          for (List<List<org.geojson.LngLatAlt>> polygonCoordinates : multiPolygonCoordinates) {
            org.geojson.Polygon polygon = new org.geojson.Polygon();
            polygon.setCoordinates(polygonCoordinates);
            org.geojson.Polygon newPolygon = GeometryUtilities.removeLastPoint(polygon);
            newMultiPolygon.add(newPolygon);
          }
          newGeometry = newMultiPolygon;
        }
        geoJson = objectMapper.writeValueAsString(newGeometry);
      } catch (JsonProcessingException e) {
        // do nothing if an error occurs
        LOGGER.error("GeometryUtilities.geometryToGeoJson failed with JsonProcessingException: ", e);
      }
    }
    return geoJson;
  }

  /**
   * Remove the last point from a polygon is necessary.
   *
   * @param polygon the polygon to check
   * @return the possibly modified polygon
   */
  private static org.geojson.Polygon removeLastPoint(org.geojson.Polygon polygon) {
    org.geojson.Polygon newPolygon = polygon;
    List<org.geojson.LngLatAlt> exteriorRing = polygon.getExteriorRing();
    if ((exteriorRing != null) && (exteriorRing.size() > 1)) {
      if (exteriorRing.get(exteriorRing.size() - 2).equals(exteriorRing.get(exteriorRing.size() - 1))) {
        exteriorRing.remove(exteriorRing.size() - 1);
        newPolygon = new org.geojson.Polygon();
        newPolygon.setExteriorRing(exteriorRing);
        List<List<org.geojson.LngLatAlt>> interiorRings = polygon.getInteriorRings();
        if (interiorRings != null) {
          for (List<org.geojson.LngLatAlt> interiorRing : interiorRings) {
            newPolygon.addInteriorRing(interiorRing);
          }
        }
      }
    }
    return newPolygon;
  }

  /**
   * Convert geometry from 0 to 360 to -180 to 180.
   *
   * @param feature the feature to normalize
   * @return the normalized feature
   */
  public static org.geojson.Feature normalizeGeometry(org.geojson.Feature feature) {
    org.geojson.Feature newFeature = null;
    if (feature != null) {
      // Assume feature has been checked for Circle geometry
      if (feature.getGeometry() instanceof org.geojson.MultiPolygon) {
        newFeature = GeometryUtilities.normalizeGeometry((org.geojson.MultiPolygon)feature.getGeometry());
      } else {
        List<gov.nasa.worldwind.geom.LatLon> latlons =
            GeometryUtilities.normalizeGeometry(GeometryUtilities.extractLatLons(feature, false));
        newFeature = GeometryUtilities.createFeature(latlons);
      }
      if (newFeature != null) {
        newFeature.setBbox(feature.getBbox());
        newFeature.setCrs(feature.getCrs());
        newFeature.setId(feature.getId());
        newFeature.setProperties(feature.getProperties());
      }
    }
    return newFeature;
  }

  /**
   * Convert a multipolygon geometry from 0 to 360 to -180 to 180.
   *
   * @param multiPolygon the multi polygon
   * @return the normalized feature
   */
  private static org.geojson.Feature normalizeGeometry(org.geojson.MultiPolygon multiPolygon) {
    org.geojson.Feature newFeature = null;
    if (multiPolygon != null) {
      List<List<List<org.geojson.LngLatAlt>>> lngLatAltListListList = multiPolygon.getCoordinates();
      org.geojson.MultiPolygon newMultiPolygon = new org.geojson.MultiPolygon();
      List<List<List<org.geojson.LngLatAlt>>> lngLatAltListListListNorm =
          new ArrayList<List<List<org.geojson.LngLatAlt>>>();
      for (List<List<org.geojson.LngLatAlt>> lngLatAltListList : lngLatAltListListList) {
        List<List<org.geojson.LngLatAlt>> lngLatAltListListNorm = new ArrayList<List<org.geojson.LngLatAlt>>();
        for (List<org.geojson.LngLatAlt> lngLatAltList : lngLatAltListList) {
          List<org.geojson.LngLatAlt> lngLatAltListNorm = new ArrayList<org.geojson.LngLatAlt>();
          for (org.geojson.LngLatAlt lngLatAlt : lngLatAltList) {
            gov.nasa.worldwind.geom.LatLon latLon =
                gov.nasa.worldwind.geom.LatLon.fromDegrees(lngLatAlt.getLatitude(), lngLatAlt.getLongitude());
            gov.nasa.worldwind.geom.LatLon latLonNorm = GeometryUtilities.normalizeGeometry(latLon);
            org.geojson.LngLatAlt lngLatAltNorm = new org.geojson.LngLatAlt(latLonNorm.getLongitude().getDegrees(),
                latLonNorm.getLatitude().getDegrees());
            lngLatAltListNorm.add(lngLatAltNorm);
          }
          lngLatAltListListNorm.add(lngLatAltListNorm);
        }
        lngLatAltListListListNorm.add(lngLatAltListListNorm);
      }
      newMultiPolygon.setCoordinates(lngLatAltListListListNorm);
      newFeature = GeometryUtilities.createFeature();
      newFeature.setGeometry(newMultiPolygon);
    }
    return newFeature;
  }

  /**
   * Normalize a list of LatLons to -180 to 180.
   *
   * @param latLons the list of LatLons to normalize
   * @return the list of normalized LatLons
   */
  private static List<gov.nasa.worldwind.geom.LatLon> normalizeGeometry(List<gov.nasa.worldwind.geom.LatLon> latLons) {
    List<gov.nasa.worldwind.geom.LatLon> newLatLons = new ArrayList<gov.nasa.worldwind.geom.LatLon>();
    for (gov.nasa.worldwind.geom.LatLon latLon : latLons) {
      newLatLons.add(GeometryUtilities.normalizeGeometry(latLon));
    }
    return newLatLons;
  }

  /**
   * Normalize a LatLono entry.
   *
   * @param latLon the LatLon entry to normalize
   * @return the normalized LatLon entry
   */
  private static gov.nasa.worldwind.geom.LatLon normalizeGeometry(gov.nasa.worldwind.geom.LatLon latLon) {
    gov.nasa.worldwind.geom.Angle latitude = latLon.getLatitude().normalizedLatitude();
    gov.nasa.worldwind.geom.Angle longitude = latLon.getLongitude().normalizedLongitude();
    return new gov.nasa.worldwind.geom.LatLon(latitude, longitude);
  }

  /**
   * Normalize an ESRI geometry.
   *
   * @param geometry the ESRI geometry to normalize
   * @return the normalized ESRI geometry
   */
  private static com.esri.core.geometry.Geometry normalizeGeometry(com.esri.core.geometry.Geometry geometry) {
    com.esri.core.geometry.Geometry newGeometry = geometry;
    if ((geometry != null) && !geometry.isEmpty()) {
      org.geojson.Feature feature = GeometryUtilities.convertToFeature(geometry);
      org.geojson.Feature newFeature = GeometryUtilities.normalizeGeometry(feature);
      newGeometry = GeometryUtilities.convertToESRIGeometry(newFeature);
    }
    return newGeometry;
  }

  /**
   * Get the centroid for a feature.
   *
   * @param feature the feature to compute the centroid for
   * @return a new feature with the centroid as the geometry
   */
  public static org.geojson.Feature getCentroidCoordinate(org.geojson.Feature feature) {
    org.geojson.Feature newFeature = new org.geojson.Feature();
    if (feature != null) {
      org.geojson.Feature fTmp = GeometryUtilities.checkForCircle(feature);
      if (fTmp.getGeometry() != null) {
        newFeature.setGeometry(GeometryUtilities.checkCrossesAntimeridian(fTmp).getGeometry());
        gov.nasa.worldwind.geom.Sector sector = gov.nasa.worldwind.geom.Sector.boundingSector(extractLatLons(newFeature, false));
        gov.nasa.worldwind.geom.LatLon latlon = sector.getCentroid();
        if (latlon != null) {
          latlon = GeometryUtilities.normalizeGeometry(latlon);
          newFeature.setGeometry((org.geojson.Point) GeometryUtilities.constructGeoJsonObject(latlon));
        }
      }
    }
    return newFeature;
  }

  /**
   * Check if two geometries cross or are on either side of the antimeridian.
   *
   * @param geom1 the first geometry to check
   * @param geom2 the second geometry to check
   * @return true if the geometries cross or are on either side of the antimeridian
   */
  private static boolean locationsCrossDateline(org.geojson.GeoJsonObject geom1, org.geojson.GeoJsonObject geom2) {
    boolean crosses = false;
    if ((geom1 != null) && (geom2 != null)) {
      if ((geom1 instanceof org.geojson.Point) && (geom2 instanceof org.geojson.Point)) {
        gov.nasa.worldwind.geom.LatLon latLon1 =
            GeometryUtilities.normalizeGeometry(GeometryUtilities.extractLatLonFromPoint(geom1));
        gov.nasa.worldwind.geom.LatLon latLon2 =
            GeometryUtilities.normalizeGeometry(GeometryUtilities.extractLatLonFromPoint(geom2));
        crosses = gov.nasa.worldwind.geom.LatLon.locationsCrossDateline(latLon1, latLon2);
      } else if ((geom1 instanceof org.geojson.Point) && !(geom2 instanceof org.geojson.Point)) {
        crosses = (geom2 instanceof org.geojson.MultiPolygon)
            ? GeometryUtilities.pointToNonPointCrosses((org.geojson.Point) geom1, (org.geojson.MultiPolygon) geom2)
            : GeometryUtilities.pointToNonPointCrosses((org.geojson.Point) geom1, geom2);
      } else if (!(geom1 instanceof org.geojson.Point) && (geom2 instanceof org.geojson.Point)) {
        crosses = (geom1 instanceof org.geojson.MultiPolygon)
            ? GeometryUtilities.pointToNonPointCrosses((org.geojson.Point) geom2, (org.geojson.MultiPolygon) geom1)
            : GeometryUtilities.pointToNonPointCrosses((org.geojson.Point) geom2, geom1);
      } else {
        crosses = GeometryUtilities.geomToGeomCrosses(geom1, geom2);
      }
    }
    return crosses;
  }

  /**
   * Get the LaLon coordinate from a Point.
   *
   * @param geoJsonObjectGeometry the Point to get the LatLon coordinate from
   * @return the LatLon coordinate
   */
  private static gov.nasa.worldwind.geom.LatLon extractLatLonFromPoint(org.geojson.GeoJsonObject geoJsonObjectGeometry) {
    String geometry = geoJsonObjectGeometry.toString();
    Pattern latPattern = Pattern.compile("latitude=[-+]?\\d*\\.\\d*");
    Pattern lonPattern = Pattern.compile("longitude=[-+]?\\d*\\.\\d*");
    Matcher latMatcher = latPattern.matcher(geometry);
    Matcher lonMatcher = lonPattern.matcher(geometry);
    latMatcher.find();
    lonMatcher.find();
    String latitude = latMatcher.group().split("=")[1];
    String longitude = lonMatcher.group().split("=")[1];
    return new gov.nasa.worldwind.geom.LatLon(gov.nasa.worldwind.geom.Angle.fromDegrees(Double.parseDouble(latitude)),
        gov.nasa.worldwind.geom.Angle.fromDegrees(Double.parseDouble(longitude)));
  }

  /**
   * Determine if a line between a point and multipolygon will cross the antimeridian.
   *
   * @param point the point to consider
   * @param multiPolygon the multipolygon to consider
   * @return true if the line crosses the antimeridian
   */
  private static boolean pointToNonPointCrosses(org.geojson.Point point, org.geojson.MultiPolygon multiPolygon) {
    boolean crosses = false;
    gov.nasa.worldwind.geom.LatLon pointLatLon =
        GeometryUtilities.normalizeGeometry(GeometryUtilities.extractLatLonFromPoint(point));
    List<List<List<org.geojson.LngLatAlt>>> lngLatAltList = multiPolygon.getCoordinates();
    for (List<List<org.geojson.LngLatAlt>> lngLatAlts : lngLatAltList) {
      org.geojson.Polygon polygon = new org.geojson.Polygon();
      polygon.setCoordinates(lngLatAlts);
      List<gov.nasa.worldwind.geom.LatLon> latLons =
          GeometryUtilities.normalizeGeometry(GeometryUtilities.extractLatLons(polygon, false));
      for (gov.nasa.worldwind.geom.LatLon latLon : latLons) {
        if (gov.nasa.worldwind.geom.LatLon.locationsCrossDateline(pointLatLon, latLon)) {
          crosses = true;
          break;
        }
      }
      if (crosses) {
        break;
      }
    }
    return crosses;
  }

  /**
   * Determine if a line between a point and geojson object will cross the antimeridian.
   *
   * @param point the point to consider
   * @param geom the geojson object to consider
   * @return true if the line crosses the antimeridian
   */
  private static boolean pointToNonPointCrosses(org.geojson.Point point, org.geojson.GeoJsonObject geom) {
    boolean crosses = false;
    gov.nasa.worldwind.geom.LatLon pointLatLon =
        GeometryUtilities.normalizeGeometry(GeometryUtilities.extractLatLonFromPoint(point));
    List<gov.nasa.worldwind.geom.LatLon> geomLatLons =
        GeometryUtilities.normalizeGeometry(GeometryUtilities.extractLatLons(geom, false));
    for (gov.nasa.worldwind.geom.LatLon latLon : geomLatLons) {
      if (gov.nasa.worldwind.geom.LatLon.locationsCrossDateline(pointLatLon, latLon)) {
        crosses = true;
        break;
      }
    }
    return crosses;
  }

  /**
   * Determine if a line between two geojson objects will cross the antimeridian.
   *
   * @param geom1 the first geojson object to consider
   * @param geom2 the second geojson object to consider
   * @return true if the line crosses the antimeridian
   */
  private static boolean geomToGeomCrosses(org.geojson.GeoJsonObject geom1, org.geojson.GeoJsonObject geom2) {
    boolean crosses = false;
    List<gov.nasa.worldwind.geom.LatLon> geom1LatLons = (geom1 instanceof org.geojson.MultiPolygon)
        ? GeometryUtilities.normalizeGeometry(GeometryUtilities.extractLatLons((org.geojson.MultiPolygon) geom1))
        : GeometryUtilities.normalizeGeometry(GeometryUtilities.extractLatLons(geom1));
    List<gov.nasa.worldwind.geom.LatLon> geom2LatLons = (geom2 instanceof org.geojson.MultiPolygon)
        ? GeometryUtilities.normalizeGeometry(GeometryUtilities.extractLatLons((org.geojson.MultiPolygon) geom2))
        : GeometryUtilities.normalizeGeometry(GeometryUtilities.extractLatLons(geom2));
    for (gov.nasa.worldwind.geom.LatLon latLon1 : geom1LatLons) {
      for (gov.nasa.worldwind.geom.LatLon latLon2 : geom2LatLons) {
        if (gov.nasa.worldwind.geom.LatLon.locationsCrossDateline(latLon1, latLon2)) {
          crosses = true;
          break;
        }
      }
      if (crosses) {
        break;
      }
    }
    return crosses;
  }

  /**
   * Compute intermediate points on a transit.
   *
   * @param point1 starting point
   * @param point2 ending point
   * @param spacingInNauticalMiles spacing between points
   * @return the intermediate points along the transit
   */
  public static List<org.geojson.Point> calculateIntermediatePointsOnGreatCircle(org.geojson.Point point1,
      org.geojson.Point point2, double spacingInNauticalMiles) {
    com.esri.core.geometry.Point esriPoint1 =
        (com.esri.core.geometry.Point) GeometryUtilities.convertToESRIGeometry(point1);
    com.esri.core.geometry.Point esriPoint2 =
        (com.esri.core.geometry.Point) GeometryUtilities.convertToESRIGeometry(point2);
    List<com.esri.core.geometry.Point> esriPoints =
        GeometryUtilities.calculateIntermediatePointsOnGreatCircle(esriPoint1, esriPoint2, spacingInNauticalMiles);
    List<org.geojson.Point> points = new ArrayList<org.geojson.Point>();
    for (com.esri.core.geometry.Point esriPoint : esriPoints) {
      org.geojson.Feature feature = GeometryUtilities.convertToFeature(esriPoint);
      points.add((org.geojson.Point) feature.getGeometry());
    }
    return points;
  }

  /**
   * Returns a list of intermediate points along the great circle connecting the two specified points. The starting
   * point and ending point are included in the list.
   *
   * @param p_point1 starting ESRI point
   * @param p_point2 ending ESRI point
   * @param spacingInNauticalMiles spacing between points
   * @return the list
   */
  private static List<com.esri.core.geometry.Point> calculateIntermediatePointsOnGreatCircle(
      com.esri.core.geometry.Point p_point1, com.esri.core.geometry.Point p_point2, double spacingInNauticalMiles) {
    if (spacingInNauticalMiles <= 0) {
      throw new OptException("spacingInNauticalMiles must be positive");
    }
    if ((EqualsUtilities.areEqual(p_point1.getY(), 90.0) && EqualsUtilities.areEqual(p_point2.getY(), -90.0))
        || (EqualsUtilities.areEqual(p_point1.getY(), -90.0) && EqualsUtilities.areEqual(p_point2.getY(), 90.0))) {
      throw new OptException("The two point latitudes cannot be antipodal");
    }
    if (EqualsUtilities.areEqual(p_point1.getY() + p_point2.getY(), 0.0)
        && EqualsUtilities.areEqual(Math.abs(p_point1.getX() - p_point2.getX()), 180.0)) {
      throw new OptException("The two point longitudes cannot be antipodal");
    }
    double lat1 = Math.toRadians(p_point1.getY());
    double lon1 = Math.toRadians(p_point1.getX());
    double lat2 = Math.toRadians(p_point2.getY());
    double lon2 = Math.toRadians(p_point2.getX());
    double greatCircleDistanceInRadians = GeometryUtilities.calculateGreatCircleDistanceInRadians(p_point1, p_point2);
    double greatCircleDistanceInNauticalMiles = greatCircleDistanceInRadians * EARTH_RADIUS_IN_NAUTICAL_MILES;
    double deltaF = spacingInNauticalMiles / greatCircleDistanceInNauticalMiles;
    List<com.esri.core.geometry.Point> geoPoints = new ArrayList<>();
    geoPoints.add(p_point1);
    for (double f = deltaF; f < 1.0; f += deltaF) {
      double A = Math.sin((1.0 - f) * greatCircleDistanceInRadians) / Math.sin(greatCircleDistanceInRadians);
      double B = Math.sin(f * greatCircleDistanceInRadians) / Math.sin(greatCircleDistanceInRadians);
      double x = A * Math.cos(lat1) * Math.cos(lon1) + B * Math.cos(lat2) * Math.cos(lon2);
      double y = A * Math.cos(lat1) * Math.sin(lon1) + B * Math.cos(lat2) * Math.sin(lon2);
      double z = A * Math.sin(lat1) + B * Math.sin(lat2);
      double lat = Math.atan2(z, Math.sqrt(x * x + y * y));
      double lon = Math.atan2(y, x);
      geoPoints.add(new com.esri.core.geometry.Point(Math.toDegrees(lon), Math.toDegrees(lat)));
    }
    geoPoints.add(p_point2);
    return geoPoints;
  }

  /**
   * Calculate the great circle distance.
   *
   * @param p_point1 ESRI starting point
   * @param p_point2 ESRI ending point
   * @return the great circle distance
   */
  private static double calculateGreatCircleDistanceInRadians(com.esri.core.geometry.Point p_point1,
      com.esri.core.geometry.Point p_point2) {
    double lat1 = Math.toRadians(p_point1.getY());
    double lon1 = Math.toRadians(p_point1.getX());
    double lat2 = Math.toRadians(p_point2.getY());
    double lon2 = Math.toRadians(p_point2.getX());
    double greatCircleDistanceInRadians =
        Math.acos(Math.sin(lat1) * Math.sin(lat2) + Math.cos(lat1) * Math.cos(lat2) * Math.cos(lon1 - lon2));
    return greatCircleDistanceInRadians;
  }

  /**
   * Calculate the great circle angular distance between two points.
   *
   * @param point1 GeoJSON starting point.
   * @param point2 GeoJSON ending point.
   * @return the great circle angular distance in radians.
   */
  public static double calculateGreatCircleAngularDistanceInRadians(org.geojson.Point point1,
      org.geojson.Point point2) {
    double lat1 = Math.toRadians(point1.getCoordinates().getLatitude());
    double lon1 = Math.toRadians(point1.getCoordinates().getLongitude());
    double lat2 = Math.toRadians(point2.getCoordinates().getLatitude());
    double lon2 = Math.toRadians(point2.getCoordinates().getLongitude());
    double greatCircleDistanceInRadians =
        Math.acos(Math.sin(lat1) * Math.sin(lat2) + Math.cos(lat1) * Math.cos(lat2) * Math.cos(lon1 - lon2));
    return greatCircleDistanceInRadians;
  }

  /**
   * Extract the Polygons out of a MultiPolygon
   * 
   * @param multiPolygon the MultiPolygon to extract the Polygons from
   * @return the list of Polygons in the MultiPolygon
   */
  public static List<org.geojson.Polygon> multiPolygonToPolygons(org.geojson.MultiPolygon multiPolygon) {
    List<org.geojson.Polygon> polygons = new ArrayList<>();
    List<List<List<org.geojson.LngLatAlt>>> coordinatesListArray = multiPolygon.getCoordinates();
    for (List<List<org.geojson.LngLatAlt>> coordinatesList : coordinatesListArray) {
      int count = 0;
      org.geojson.Polygon polygon = new org.geojson.Polygon();
      for (List<org.geojson.LngLatAlt> coordinates : coordinatesList) {
        if (count++ == 0) {
          polygon.setExteriorRing(coordinates);
        } else {
          polygon.addInteriorRing(coordinates);
        }
      }
      polygons.add(polygon);
    }
    return polygons;
  }
  
  /**
   * Replace the first point of a line in a feature with the passed in point.  An empty feature is returned if the
   * geometry is null or if the geometry is not a LineString or a MultiLineString
   * 
   * @param feature the feature containing the line
   * @param point the replacement point
   * @return a new feature with a new geometry having the first point replaced
   */
  public static org.geojson.Feature replaceFirstPointInLine(org.geojson.Feature feature, org.geojson.Point point) {
    org.geojson.Feature newFeature = new org.geojson.Feature();
    org.geojson.GeoJsonObject geometry = feature.getGeometry();
    if (geometry != null) {
      org.geojson.LngLatAlt newPointCoordinates = new org.geojson.LngLatAlt();
      org.geojson.LngLatAlt pointCoordinates = point.getCoordinates();
      newPointCoordinates.setLatitude(pointCoordinates.getLatitude());
      newPointCoordinates.setLongitude(pointCoordinates.getLongitude());
      if (geometry instanceof org.geojson.LineString) {
        org.geojson.LineString lineString = (org.geojson.LineString)geometry;
        List<org.geojson.LngLatAlt> coordinates = lineString.getCoordinates();
        org.geojson.LineString newLineString = new org.geojson.LineString();
        if ((coordinates != null) && (coordinates.size() > 0)) {
          newLineString.add(newPointCoordinates);
          for (int index=1; index < coordinates.size(); index++) {
            newLineString.add(coordinates.get(index));
          }
          newFeature.setGeometry(newLineString);
        }
      }
      else if (geometry instanceof org.geojson.MultiLineString) {
        // Not having better information, replace the first point of the MultiLineString
        org.geojson.MultiLineString multiLineString = (org.geojson.MultiLineString)geometry;
        List<List<org.geojson.LngLatAlt>> coordinatesList = multiLineString.getCoordinates();
        org.geojson.MultiLineString newMultiLineString = new org.geojson.MultiLineString();
        if ((coordinatesList != null) && (coordinatesList.size() > 0)) {
          List<List<org.geojson.LngLatAlt>> newCoordinatesList = new ArrayList<List<org.geojson.LngLatAlt>>();
          List<org.geojson.LngLatAlt> newCoordinates = new ArrayList<org.geojson.LngLatAlt>();
          newCoordinates.add(newPointCoordinates);
          List<org.geojson.LngLatAlt> coordinates = coordinatesList.get(0);
          for (int index=1; index < coordinates.size(); index++) {
            newCoordinates.add(coordinates.get(index));
          }
          newCoordinatesList.add(newCoordinates);
          for (int index=1; index < coordinatesList.size(); index++) {
            newCoordinatesList.add(coordinatesList.get(index));
          }
          newMultiLineString.setCoordinates(newCoordinatesList);
          newFeature.setGeometry(newMultiLineString);
        }
      }
    }
    return newFeature;
  }
  
  /**
   * Replace the last point of a line in a feature with the passed in point.  An empty feature is returned if the
   * geometry is null or if the geometry is not a LineString or a MultiLineString
   * 
   * @param feature the feature containing the line
   * @param point the replacement point
   * @return a new feature with a new geometry having the last point replaced
   */
  public static org.geojson.Feature replaceLastPointInLine(org.geojson.Feature feature, org.geojson.Point point) {
    org.geojson.Feature newFeature = new org.geojson.Feature();
    org.geojson.GeoJsonObject geometry = feature.getGeometry();
    if (geometry != null) {
      org.geojson.LngLatAlt newPointCoordinates = new org.geojson.LngLatAlt();
      org.geojson.LngLatAlt pointCoordinates = point.getCoordinates();
      newPointCoordinates.setLatitude(pointCoordinates.getLatitude());
      newPointCoordinates.setLongitude(pointCoordinates.getLongitude());
      if (geometry instanceof org.geojson.LineString) {
        org.geojson.LineString lineString = (org.geojson.LineString)geometry;
        List<org.geojson.LngLatAlt> coordinates = lineString.getCoordinates();
        org.geojson.LineString newLineString = new org.geojson.LineString();
        if ((coordinates != null) && (coordinates.size() > 0)) {
          for (org.geojson.LngLatAlt lngLatAlt : coordinates) {
            org.geojson.LngLatAlt newLngLatAlt = new org.geojson.LngLatAlt();
            newLngLatAlt.setLatitude(lngLatAlt.getLatitude());
            newLngLatAlt.setLongitude(lngLatAlt.getLongitude());
            newLineString.add(newLngLatAlt);
          }
          newLineString.getCoordinates().set(coordinates.size() - 1, newPointCoordinates);
          newFeature.setGeometry(newLineString);
        }
      }
      else if (geometry instanceof org.geojson.MultiLineString) {
        // Not having better information, replace the last point of the MultiLineString
        org.geojson.MultiLineString multiLineString = (org.geojson.MultiLineString)geometry;
        List<List<org.geojson.LngLatAlt>> coordinatesList = multiLineString.getCoordinates();
        org.geojson.MultiLineString newMultiLineString = new org.geojson.MultiLineString();
        if ((coordinatesList != null) && (coordinatesList.size() > 0)) {
          List<List<org.geojson.LngLatAlt>> newCoordinatesList = new ArrayList<List<org.geojson.LngLatAlt>>();
          List<org.geojson.LngLatAlt> newCoordinates = new ArrayList<org.geojson.LngLatAlt>();
          List<org.geojson.LngLatAlt> coordinates = coordinatesList.get(coordinatesList.size() - 1);
          for (org.geojson.LngLatAlt coordinate : coordinates) {
            org.geojson.LngLatAlt newLngLatAlt = new org.geojson.LngLatAlt();
            newLngLatAlt.setLatitude(coordinate.getLatitude());
            newLngLatAlt.setLongitude(coordinate.getLongitude());
            newCoordinates.add(newLngLatAlt);
          }
          newCoordinates.set(coordinates.size() - 1, newPointCoordinates);
          for (List<org.geojson.LngLatAlt> coords : coordinatesList) {
            List<org.geojson.LngLatAlt> newCoords = new ArrayList<org.geojson.LngLatAlt>();
            for (org.geojson.LngLatAlt coord : coords) {
              org.geojson.LngLatAlt newLngLatAlt = new org.geojson.LngLatAlt();
              newLngLatAlt.setLatitude(coord.getLatitude());
              newLngLatAlt.setLongitude(coord.getLongitude());
              newCoords.add(newLngLatAlt);
            }
            newCoordinatesList.add(newCoords);
          }
          newCoordinatesList.set(coordinatesList.size() - 1, newCoordinates);
          newMultiLineString.setCoordinates(newCoordinatesList);
          newFeature.setGeometry(newMultiLineString);
        }
      }
    }
    return newFeature;
  }
  
  /**
   * Make sure valid geometry is presented
   * 
   * @param feature feature with geometry to validate
   * @exception ValidationException
   */
  public static void validateGeometry(org.geojson.Feature feature) {
    if (feature != null) {
      org.geojson.GeoJsonObject geometry = feature.getGeometry();
      if (geometry != null) {  // geometry is allowed to be null for some features
        if ((geometry instanceof org.geojson.Polygon) || (geometry instanceof org.geojson.MultiPolygon)) {
          // make sure that the area can be calculated properly
          GeometryUtilities.calculateArea(feature, false);
        } else if ((geometry instanceof org.geojson.LineString) || (geometry instanceof org.geojson.MultiLineString)) {
          // make sure that the length can be calculated properly
          GeometryUtilities.calculateLength2D(feature);
        } else if (geometry instanceof org.geojson.MultiPoint) {
          // Nothing to validate in this case
        } else if (geometry instanceof org.geojson.Point) {
          // See if this is a circle and if a radius is provided
          if (Constants.CIRCLE_TYPE.equalsIgnoreCase(feature.getProperty(Constants.OPT_SHAPETYPE))) {
            // Get the radius and use that to build a circle
            double radius = 0.0;
            Object radiusObj = feature.getProperty(Constants.OPT_SYS_RADIUS);
            if (radiusObj instanceof Double) {
              radius = (double)radiusObj;
            } else if (radiusObj instanceof Integer) {
              radius = ((Integer)radiusObj).intValue();
            }
            if (radius <= 0.0) {
              throw new ValidationException("Invalid geometry for feature: " + feature.getProperty(Constants.OPT_SYS_NAME) +
                  " - does not have a positive value as radius");
            }
          }
        } else {
          throw new ValidationException("Invalid geometry type for feature: " + feature.getProperty(Constants.OPT_SYS_NAME));
        }
      }
    }
    else {
      throw new ValidationException("Feature to Validate is null");
    }
  }

  /**
   * Build the polygon representing a sector.  innerradius of zero and begin/end arc is 0/360 is a circle
   *
   * @param center the point from which the sector is offset
   * @param innerradius inner radius in NM
   * @param outerradius outer radius in NM
   * @param beginarc starting bearing
   * @param endarc ending bearing
   * @param circle true if this is supposed to be circle
   * @return sector polygon
   */
  public static org.geojson.Polygon buildSectorGeometry(org.geojson.Point center, double innerradius, double outerradius,
      double beginarc, double endarc, boolean circle) {
    final double DIVS = 36.0;
    final double AZINC = 360.0 / DIVS;
    double az = 0.0;
    org.geojson.Polygon geometry = null;
    double innerRadius = innerradius * GeometryUtilities.KILOMETERSPERNAUTICALMILE;
    double outerRadius = outerradius * GeometryUtilities.KILOMETERSPERNAUTICALMILE;
    List<org.geojson.LngLatAlt> points = new ArrayList<org.geojson.LngLatAlt>();
    while (endarc < beginarc) {
      endarc += 360.0;
    }
    // Get the starting point - if inner radius is 0.0, the center is the starting point.  But only if this isn't a
    // circle
    if (!circle) {
      if (EqualsUtilities.areEqual(innerRadius, 0.0)) {
        points.add(center.getCoordinates());
      } else {
        org.geojson.Point point = GeometryUtilities.pointFromRangeBearing(center, innerRadius, beginarc);
        points.add(point.getCoordinates());
      }
    }
    az = beginarc;
    while (az < endarc) {
      org.geojson.Point point = GeometryUtilities.pointFromRangeBearing(center, outerRadius, az);
      points.add(point.getCoordinates());
      az += AZINC;
    }
    org.geojson.Point point = GeometryUtilities.pointFromRangeBearing(center, outerRadius, endarc);
    points.add(point.getCoordinates());
    if (innerRadius > 0.0) {
      az = endarc;
      while (az > beginarc) {
        point = GeometryUtilities.pointFromRangeBearing(center, innerRadius, az);
        points.add(point.getCoordinates());
        az -= AZINC;
      }
      point = GeometryUtilities.pointFromRangeBearing(center, innerRadius, beginarc);
      points.add(point.getCoordinates());
    }
    geometry = new org.geojson.Polygon(points);
    return geometry;
  }
  
  /**
   * Check a list of features to see if any are circles and return them with polygon geometries representing the circle
   * 
   * @param features the list of features to check
   * @return a list of features with all features and those that are circles have the geometry replaced with a polygon
   */
  private static List<org.geojson.Feature> checkForCircles(List<org.geojson.Feature> features) {
    List<org.geojson.Feature> newFeatures = null;
    if (features != null) {
      newFeatures = new ArrayList<>();
      for (org.geojson.Feature feature : features) {
        newFeatures.add(GeometryUtilities.checkForCircle(feature));
      }
    }
    return newFeatures;
  }
  
  /**
   * Check a feature to see if it has a circle geometry - return it with a polygon geometry representing the circle
   * 
   * @param feature the feature to check
   * @return a new feature with any circle geometry replaced with a polygon representing the circle
   */
  private static org.geojson.Feature checkForCircle(org.geojson.Feature feature) {
    org.geojson.Feature newFeature = null;
    if (feature != null) {
      newFeature = new org.geojson.Feature();
      newFeature.setProperties(feature.getProperties());
      org.geojson.GeoJsonObject featureGeometry = feature.getGeometry();
      org.geojson.GeoJsonObject polygon = featureGeometry;
      if (featureGeometry != null) {
        // If this is a circle, build a polygon that represents it.
        if ((featureGeometry instanceof org.geojson.Point) && Constants.CIRCLE_TYPE.equals(feature.getProperty(Constants.OPT_SHAPETYPE))) {
          // Get the radius and use that to build a circle
          double radius = 0.0;
          Object radiusObj = feature.getProperty(Constants.OPT_SYS_RADIUS);
          if (radiusObj instanceof Double) {
            radius = (double)radiusObj;
          } else if (radiusObj instanceof Integer) {
            radius = ((Integer)radiusObj).intValue();
          }
          polygon = GeometryUtilities.buildSectorGeometry((org.geojson.Point)featureGeometry, 0.0, radius, 0.0, 360.0, true);
        }
      }
      newFeature.setGeometry(polygon);
    }
    return newFeature;
  }
  
  /**
   * Convert the decimal latitude/longitude to DDMMSS or DDDMMSS based on the value of isLon
   * 
   * @param lat latitude or longitude decimal value
   * @param isLon true if this is a longitude
   * @return string representing the latitude/longitude as DDMMSS or DDDMMSS
   */
  public static String decimalToDDMMSS(double lat, boolean isLon) {
    String direction = "N";
    if (lat < 0.0) {
      direction = "S";
      lat = Math.abs(lat);
    }
    double deg = Math.floor(lat);
    double min = Math.floor((lat - deg) * 60);
    double sec = ((lat - deg - min / 60) * 3600);
    if (sec >= 60.0) {
      sec = 0;
      min++;
    }
    if (min >= 60.0) {
      min = 0.0;
      deg++;
    }
    String degree = String.valueOf(deg);
    if (deg < 10.0) {
      degree = (isLon) ? ("00" + degree) : ("0" + degree);
    } else if (deg < 100.0) {
      degree = (isLon) ? ("0" + degree) : degree;
    }
    String minute = String.valueOf(min);
    if (min < 10.0) {
      minute = "0" + minute;
    }
    String second = String.valueOf(sec);
    if (sec < 10.0) {
      second = "0" + second;
    }
    String ddmmss = degree + "\u00B0" + minute + "'" + second + "\"" + direction;
    return ddmmss;
  }
  
  /**
   * Convert the decimal latitude/longitude to DDMMSS or DDDMMSS based on the value of isLon
   * 
   * @param lat latitude or longitude decimal value
   * @param isLon true if this is a longitude
   * @return string representing the latitude/longitude as 'Direction DD MM SS' or 'Direction DDD MM SS'
   */
  public static String decimalToDirectionDDMMSS(double lat, boolean isLon) {
    String direction = "";
    if (isLon) {
      direction = "E";
      if (lat < 0.0) {
        direction = "W";
        lat = Math.abs(lat);
      }
    }
    else {
      direction = "N";
      if (lat < 0.0) {
        direction = "S";
        lat = Math.abs(lat);
      }      
    }
    
    double deg = Math.floor(lat);
    double min = Math.floor((lat - deg) * 60);
    double sec = ((lat - deg - min / 60) * 3600);
    if (sec >= 60.0) {
      sec = 0;
      min++;
    }
    if (min >= 60.0) {
      min = 0.0;
      deg++;
    }
        
    String degree = String.valueOf(deg);
    degree = degree.substring(0, degree.indexOf("."));
    if (deg < 10.0) {
      degree = (isLon) ? ("00" + degree) : ("0" + degree);
    } else if (deg < 100.0) {
      degree = (isLon) ? ("0" + degree) : degree;
    }
    
    String minute = String.valueOf(min);
    minute = minute.substring(0, minute.indexOf("."));
    if (min < 10.0) {
      minute = "0" + minute;
    }
    String second = String.valueOf(sec);
    second = second.substring(0, second.indexOf("."));
    if (sec < 10.0) {
      second = "0" + second;
    }
    
    String ddmmss = direction + " " + degree + " " + minute + " " + second;    
    return ddmmss;
  }
 
  
  /**
   * Find the radius of a Circle represented by a polygon by getting the envelope of the poiygon and then finding the
   * distance from the centroid to the line forming the bottom of the envelope 
   * 
   * @param feature the feature with a Polygon geometry
   * @return the radius of the polygon or -1.0 if the geometry is not a Polygon
   */
  public static double findRadiusOfPolygon(org.geojson.Feature feature) {
    double radius = -1.0;
    org.geojson.GeoJsonObject geometry = feature.getGeometry();
    if ((geometry != null) && (Constants.POLYGON_TYPE.equalsIgnoreCase(geometry.getType()))) {
      org.geojson.Point point1 = feature.getProperty(Constants.JSON_CENTROID_GEOMETRY);
      // find the envelope of the polygon
      org.geojson.Feature envelope = GeometryUtilities.getEnvelope(feature);
      // the lowerLeft will give us the latitude, the centroid will give us the longitude
      org.geojson.LngLatAlt lowerLeft = envelope.getProperty("lowerLeft");
      org.geojson.Point point2 = new org.geojson.Point( point1.getCoordinates().getLongitude(), lowerLeft.getLatitude());
      radius = GeometryUtilities.rangeFromPoints(point1, point2);
    }
    return radius;
  }
  
  /**
   * Ensure that the value is between -180 and +180
   * 
   * @param degrees the value to normalize
   * @return the normalized value
   */
  public static double normalizeDegrees(double degrees) {
    while (degrees > 180.0) {
      degrees -= 360.0;
    }
    while (degrees < -180.0) {
      degrees += 360.0;
    }
    return degrees;
  }
}
