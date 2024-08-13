///////////////////////////////////////////////////////////////////////////
// Copyright © Esri. All Rights Reserved.
//
// Licensed under the Apache License Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
///////////////////////////////////////////////////////////////////////////

define([
  'esri/SpatialReference',
  'esri/geometry/Polygon',
  'esri/geometry/Polyline',
  'esri/geometry/geometryEngine',
  'esri/geometry/webMercatorUtils'
],

  function (SpatialReference, Polygon, Polyline, GeometryEngine, webMercatorUtils
  ) {
    var mo = {};

    mo.getGeodesicDistance = function (selectedGeometry, intersectedGeomery, distanceUnit) {
      var geodesicDistance = 0;
      //retrun 0 if any of the geometry is not valid
      //or if the geometires intersects
      if (!selectedGeometry || !intersectedGeomery ||
        GeometryEngine.intersects(selectedGeometry, intersectedGeomery)) {
        return 0;
      }
      var outSR = new SpatialReference(4326);
      if (webMercatorUtils.canProject(selectedGeometry, outSR) &&
        webMercatorUtils.canProject(intersectedGeomery, outSR)) {
        //project geometries in 4326 SR
        selectedGeometry = webMercatorUtils.project(selectedGeometry, outSR);
        intersectedGeomery = webMercatorUtils.project(intersectedGeomery, outSR);
        //Based on geometry types call the method to get distance between two geometries
        if (selectedGeometry.type === "point") {
          if (intersectedGeomery.type === 'point') {
            geodesicDistance =
              mo.getGeodesicDistanceBetweenPointToPoint(selectedGeometry, intersectedGeomery);
          } else if (intersectedGeomery.type === 'polyline') {
            geodesicDistance =
              mo.getGeodesicDistanceBetweenPointToLine(selectedGeometry, intersectedGeomery);
          } else if (intersectedGeomery.type === 'polygon') {
            geodesicDistance =
              mo.getGeodesicDistanceBetweenPointToPolygon(selectedGeometry, intersectedGeomery);
          }
        } else if (selectedGeometry.type === "polyline") {
          if (intersectedGeomery.type === 'point') {
            geodesicDistance =
              mo.getGeodesicDistanceBetweenPointToLine(intersectedGeomery, selectedGeometry);
          } else if (intersectedGeomery.type === 'polyline') {
            geodesicDistance =
              mo.getGeodesicDistanceBetweenLineToLine(selectedGeometry, intersectedGeomery);
          } else if (intersectedGeomery.type === 'polygon') {
            geodesicDistance =
              mo.getGeodesicDistanceBetweenLineToPolygon(selectedGeometry, intersectedGeomery);
          }
        } else if (selectedGeometry.type === "polygon") {
          if (intersectedGeomery.type === 'point') {
            geodesicDistance =
              mo.getGeodesicDistanceBetweenPointToPolygon(intersectedGeomery, selectedGeometry);
          } else if (intersectedGeomery.type === 'polyline') {
            geodesicDistance =
              mo.getGeodesicDistanceBetweenLineToPolygon(intersectedGeomery, selectedGeometry);
          } else if (intersectedGeomery.type === 'polygon') {
            geodesicDistance =
              mo.getGeodesicDistanceBetweenPolygonToPolygon(selectedGeometry, intersectedGeomery);
          }
        }
        //By default geodesicDistance will be calculated in meters,
        //If distance unit is defined, convert the the the value from meters to the defined unit
        if (distanceUnit) {
          geodesicDistance = mo.convertSingle(geodesicDistance, distanceUnit);
        }
      }
      return geodesicDistance;
    };

    mo.getGeodesicDistanceBetweenPolygonToPolygon = function (polygon1, polygon2) {
      var newDistance, minimumDistance;
      //if two geometries intersects then the distance will be zero
      if (GeometryEngine.intersects(polygon1, polygon2)) {
        return 0;
      }
      var i, j;
      for (i = 0; i < polygon1.rings.length; i++) {
        for (j = 0; j < polygon1.rings[i].length; j++) {
          var polygon1Point = polygon1.getPoint(i, j);
          newDistance = mo.getGeodesicDistanceBetweenPointToPolygon(polygon1Point, polygon2);
          if (minimumDistance === undefined || newDistance < minimumDistance) {
            minimumDistance = newDistance;
          }
        }
      }
      for (i = 0; i < polygon2.rings.length; i++) {
        for (j = 0; j < polygon2.rings[i].length; j++) {
          var polygon2Point = polygon2.getPoint(i, j);
          newDistance = mo.getGeodesicDistanceBetweenPointToPolygon(polygon2Point, polygon1);
          if (minimumDistance === undefined || newDistance < minimumDistance) {
            minimumDistance = newDistance;
          }
        }
      }
      return minimumDistance;
    };

    mo.getGeodesicDistanceBetweenLineToPolygon = function (line, polygon) {
      var newDistance, minimumDistance, i, j;
      //if two geometries intersects then the distance will be zero
      if (GeometryEngine.intersects(line, polygon)) {
        return 0;
      }
      for (i = 0; i < line.paths.length; i++) {
        for (j = 0; j < line.paths[i].length; j++) {
          var line1Point = line.getPoint(i, j);
          newDistance = mo.getGeodesicDistanceBetweenPointToPolygon(line1Point, polygon);
          if (minimumDistance === undefined || newDistance < minimumDistance) {
            minimumDistance = newDistance;
          }
        }
      }
      for (i = 0; i < polygon.rings.length; i++) {
        for (j = 0; j < polygon.rings[i].length; j++) {
          var polygon1Point = polygon.getPoint(i, j);
          newDistance = mo.getGeodesicDistanceBetweenPointToLine(polygon1Point, line);
          if (minimumDistance === undefined || newDistance < minimumDistance) {
            minimumDistance = newDistance;
          }
        }
      }
      return minimumDistance;
    };

    mo.getGeodesicDistanceBetweenLineToLine = function (line1, line2) {
      var i, j, newDistance, minimumDistance;
      //if two geometries intersects then the distance will be zero
      if (GeometryEngine.intersects(line1, line2)) {
        return 0;
      }
      for (i = 0; i < line1.paths.length; i++) {
        for (j = 0; j < line1.paths[i].length; j++) {
          var line1Point = line1.getPoint(i, j);
          newDistance = mo.getGeodesicDistanceBetweenPointToLine(line1Point, line2);
          if (minimumDistance === undefined || newDistance < minimumDistance) {
            minimumDistance = newDistance;
          }
        }
      }

      for (i = 0; i < line2.paths.length; i++) {
        for (j = 0; j < line2.paths[i].length; j++) {
          var line2Point = line2.getPoint(i, j);
          newDistance = mo.getGeodesicDistanceBetweenPointToLine(line2Point, line1);
          if (minimumDistance === undefined || newDistance < minimumDistance) {
            minimumDistance = newDistance;
          }
        }
      }
      return minimumDistance;
    };

    mo.getGeodesicDistanceBetweenPointToPolygon = function (point, polygon) {
      var i, newDistance, minimumDistance;
      if (polygon.rings.length > 0) {
        for (i = 0; i < polygon.rings.length; i++) {
          var polygonJson = {
            "rings": [polygon.rings[i]],
            "spatialReference": polygon.spatialReference
          };
          var newPolygon = new Polygon(polygonJson);
          var nearestCoordinate = GeometryEngine.nearestCoordinate(newPolygon, point).coordinate;
          newDistance = mo.getGeodesicDistanceBetweenPointToPoint(point, nearestCoordinate);
          if (minimumDistance === undefined || newDistance < minimumDistance) {
            minimumDistance = newDistance;
          }
        }
      }
      return minimumDistance;
    };

    mo.getGeodesicDistanceBetweenPointToLine = function (point, line) {
      var i, newDistance, minimumDistance;
      if (line.paths.length > 0) {
        for (i = 0; i < line.paths.length; i++) {
          var polylineJson = {
            "paths": [line.paths[i]],
            "spatialReference": line.spatialReference
          };
          var polyline = new Polyline(polylineJson);
          var nearestCoordinate = GeometryEngine.nearestCoordinate(polyline, point).coordinate;
          newDistance = mo.getGeodesicDistanceBetweenPointToPoint(point, nearestCoordinate);
          if (minimumDistance === undefined || newDistance < minimumDistance) {
            minimumDistance = newDistance;
          }
        }
      }
      return minimumDistance;
    };

    mo.getGeodesicDistanceBetweenPointToPoint = function (point1, point2) {
      //Using Vincenty Direct and Inverse Solution of Geodesics on the Ellipsoid
      var geodesicDist = mo.getDistance(point1, point2);
      return isNaN(geodesicDist) ? 0 : geodesicDist;
    };

    /*---------------------------------------------------------------------------------------------
    * Vincenty Direct and Inverse Solution of Geodesics on the Ellipsoid (c) Chris Veness 2002-2016
    *                                                                                   MIT Licence
    *
    * www.movable-type.co.uk/scripts/latlong-vincenty.html
    * www.movable-type.co.uk/scripts/geodesy/docs/module-latlon-vincenty.html
    *---------------------------------------------------------------------------------------------
    * Returns the Distance using Vincenty inverse solution.
    **/
    mo.getDistance = function (startPoint, endPoint) {
      var φ1 = startPoint.y * Math.PI / 180, λ1 = startPoint.x * Math.PI / 180;
      var φ2 = endPoint.y * Math.PI / 180, λ2 = endPoint.x * Math.PI / 180;

      var a = 6378137; var b = 6356752.314245; var f = (a - b) / a;

      var L = λ2 - λ1;
      var tanU1 = (1 - f) * Math.tan(φ1),
        cosU1 = 1 / Math.sqrt((1 + tanU1 * tanU1)), sinU1 = tanU1 * cosU1;
      var tanU2 = (1 - f) * Math.tan(φ2),
        cosU2 = 1 / Math.sqrt((1 + tanU2 * tanU2)), sinU2 = tanU2 * cosU2;

      var sinλ, cosλ, sinSqσ, sinσ, cosσ, σ, sinα, cosSqα, cos2σM, C;

      var λ = L, λʹ, iterations = 0;
      do {
        sinλ = Math.sin(λ);
        cosλ = Math.cos(λ);
        sinSqσ = (cosU2 * sinλ) * (cosU2 * sinλ) +
          (cosU1 * sinU2 - sinU1 * cosU2 * cosλ) * (cosU1 * sinU2 - sinU1 * cosU2 * cosλ);
        sinσ = Math.sqrt(sinSqσ);
        if (sinσ == 0) { // jshint ignore:line
          return { distance: 0, initialBearing: 0, finalBearing: 0 };  // co-incident points
        }
        cosσ = sinU1 * sinU2 + cosU1 * cosU2 * cosλ;
        σ = Math.atan2(sinσ, cosσ);
        sinα = cosU1 * cosU2 * sinλ / sinσ;
        cosSqα = 1 - sinα * sinα;
        cos2σM = cosσ - 2 * sinU1 * sinU2 / cosSqα;
        if (isNaN(cos2σM)) {
          cos2σM = 0;  // equatorial line: cosSqα=0 (§6)
        }
        C = f / 16 * cosSqα * (4 + f * (4 - 3 * cosSqα));
        λʹ = λ;
        λ = L + (1 - C) * f * sinα *
          (σ + C * sinσ * (cos2σM + C * cosσ * (-1 + 2 * cos2σM * cos2σM)));
      } while (Math.abs(λ - λʹ) > 1e-12 && ++iterations < 200);
      if (iterations >= 200) {
        return null;
      }

      var uSq = cosSqα * (a * a - b * b) / (b * b);
      var A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));
      var B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));
      var Δσ = B * sinσ * (cos2σM + B / 4 * (cosσ * (-1 + 2 * cos2σM * cos2σM) -
        B / 6 * cos2σM * (-3 + 4 * sinσ * sinσ) * (-3 + 4 * cos2σM * cos2σM)));

      var s = b * A * (σ - Δσ);

      //commenting the below code as we dont need the angles
      // var α1 = Math.atan2(cosU2 * sinλ, cosU1 * sinU2 - sinU1 * cosU2 * cosλ);
      // var α2 = Math.atan2(cosU1 * sinλ, -sinU1 * cosU2 + cosU1 * sinU2 * cosλ);

      //α1 = (α1 + 2 * Math.PI) % (2 * Math.PI); // normalize to 0..360
      //α2 = (α2 + 2 * Math.PI) % (2 * Math.PI); // normalize to 0..360

      s = Number(s.toFixed(3)); // round to 1mm precision
      //α1 = α1 * 180 / Math.PI;
      //α2 = α2 * 180 / Math.PI;
      return s;//{ distance: s, initialBearing: α1, finalBearing: α2 };
    };


    mo.convertSingle = function (fromValue, toUnits) {
      var fromUnits = 'meters';
      if (fromUnits === toUnits) {
        return +fromValue;
      } else {
        return (+fromValue / mo.perMeter(fromUnits)) * mo.perMeter(toUnits);
      }
    };

    mo.perMeter = function (units) {
      var conversionFactor = 1.0;
      switch (units) {
        case 'meters':
          conversionFactor = 1.0;
          break;
        case 'kilometers':
          conversionFactor = 1 / 1000;
          break;
        case 'feet':
          conversionFactor = 1 / 0.3048;
          break;
        case 'yards':
          conversionFactor = 1 / 0.9144;
          break;
        case 'miles':
          conversionFactor = 1 / 1609.344;
          break;
      }
      return conversionFactor;
    };

    return mo;
  });
