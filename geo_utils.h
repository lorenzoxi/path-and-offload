#ifndef INC_005_CCNC_DRONE_ENERGY_OPT_GEO_UTILS_H
#define INC_005_CCNC_DRONE_ENERGY_OPT_GEO_UTILS_H

#include "Point.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

inline double distance(const Point p1, const Point &p2) {

    double lat1 = p1.lat;
    double lon1 = p1.lon;
    double lat2 = p2.lat;
    double lon2 = p2.lon;

    const double R = 6371.0; // Earth's radius in kilometers

    // Convert latitude and longitude from degrees to radians
    double lat1_rad = lat1 * M_PI / 180.0;
    double lon1_rad = lon1 * M_PI / 180.0;
    double lat2_rad = lat2 * M_PI / 180.0;
    double lon2_rad = lon2 * M_PI / 180.0;

    // Haversine formula
    double dlat = lat2_rad - lat1_rad;
    double dlon = lon2_rad - lon1_rad;

    double a = sin(dlat / 2.0) * sin(dlat / 2.0) +
               cos(lat1_rad) * cos(lat2_rad) *
               sin(dlon / 2.0) * sin(dlon / 2.0);

    double c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a));

    return R * c > 0 ? R * c : 0.1; // Ensure non-negative distance

}

inline double
distance(const double &lat1, const double &lon1, const double &lat2, const double &lon2) {


    const double R = 6371.0; // Earth's radius in kilometers

    // Convert latitude and longitude from degrees to radians
    double lat1_rad = lat1 * M_PI / 180.0;
    double lon1_rad = lon1 * M_PI / 180.0;
    double lat2_rad = lat2 * M_PI / 180.0;
    double lon2_rad = lon2 * M_PI / 180.0;

    // Haversine formula
    double dlat = lat2_rad - lat1_rad;
    double dlon = lon2_rad - lon1_rad;

    double a = sin(dlat / 2.0) * sin(dlat / 2.0) +
               cos(lat1_rad) * cos(lat2_rad) *
               sin(dlon / 2.0) * sin(dlon / 2.0);

    double c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a));

    return R * c > 0 ? R * c : 0.1; // Ensure non-negative distance

}

/*
inline double
distance(const double &lat1, const double &lon1, const double &alt1,
         const double &lat2, const double &lon2, const double &alt2) {

    const double R = 6371000.0; // Earth's radius in meters

    // Convert latitude and longitude from degrees to radians
    double lat1_rad = lat1 * M_PI / 180.0;
    double lon1_rad = lon1 * M_PI / 180.0;
    double lat2_rad = lat2 * M_PI / 180.0;
    double lon2_rad = lon2 * M_PI / 180.0;

    // Haversine formula for horizontal distance
    double dlat = lat2_rad - lat1_rad;
    double dlon = lon2_rad - lon1_rad;

    double a = sin(dlat / 2.0) * sin(dlat / 2.0) +
               cos(lat1_rad) * cos(lat2_rad) *
               sin(dlon / 2.0) * sin(dlon / 2.0);

    double c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a));

    double horizontal_distance = R * c;

    // Altitude difference
    double dalt = alt2 - alt1;

    // Pythagorean distance considering altitude
    double total_distance = sqrt(horizontal_distance * horizontal_distance + dalt * dalt);

    return total_distance > 0 ? total_distance : 0.1; // Ensure non-zero, non-negative distance
}
*/

#include <cmath>

inline double
distance(const double &lat1, const double &lon1, const double &alt1,
         const double &lat2, const double &lon2, const double &alt2) {

    const double R = 6371000.0; // Earth's radius in meters

    // Convert degrees to radians
    auto to_rad = [](double deg) {
        return deg * M_PI / 180.0;
    };

    double lat1_rad = to_rad(lat1);
    double lon1_rad = to_rad(lon1);
    double lat2_rad = to_rad(lat2);
    double lon2_rad = to_rad(lon2);

    // Convert geodetic coordinates to 3D Cartesian coordinates (ECEF)
    double x1 = (R + alt1) * cos(lat1_rad) * cos(lon1_rad);
    double y1 = (R + alt1) * cos(lat1_rad) * sin(lon1_rad);
    double z1 = (R + alt1) * sin(lat1_rad);

    double x2 = (R + alt2) * cos(lat2_rad) * cos(lon2_rad);
    double y2 = (R + alt2) * cos(lat2_rad) * sin(lon2_rad);
    double z2 = (R + alt2) * sin(lat2_rad);

    // Euclidean distance in 3D space
    double dx = x2 - x1;
    double dy = y2 - y1;
    double dz = z2 - z1;

    double total_distance = std::sqrt(dx * dx + dy * dy + dz * dz);

    return total_distance > 0 ? total_distance : 0.1; // Ensure small non-zero return
}


inline double euc_distance(const Point &p1, const Point &p2) {
    double dx = p2.lat - p1.lat;
    double dy = p2.lon - p1.lon;
    double dz = p2.z - p1.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}


inline double euc_distance(const Sensor &s1, const Point &s2) {
    double dx = s2.lat - s1.lat;
    double dy = s2.lon - s1.lon;
    double dz = s2.z - s1.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

#endif //INC_005_CCNC_DRONE_ENERGY_OPT_GEO_UTILS_H
