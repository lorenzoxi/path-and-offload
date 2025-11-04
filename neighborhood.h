#ifndef INC_005_CCNC_DRONE_ENERGY_OPT_NEIGHBORHOOD_H
#define INC_005_CCNC_DRONE_ENERGY_OPT_NEIGHBORHOOD_H

#include <vector>
#include "Sensor.h"
#include "Point.h"
#include "geo_utils.h"
#include <iostream>
#include <tuple>

inline std::vector<int>
neighborhood(const std::vector<Sensor> &sensors, const Point &p, double radius_min,
             double radius_max) {

    double lat1 = p.lat;
    double lon1 = p.lon;

    std::vector<int> neighbors = {};


    for (const auto &sensor: sensors) {

        double lat2 = sensor.lat;
        double lon2 = sensor.lon;

        double dist_in_km = distance(lat1, lon1, lat2, lon2);

        //std::cout << "Sensor ID: " << sensor._id << ", Distance: " << dist_in_km << std::endl;

        if (dist_in_km >= radius_min && dist_in_km <= radius_max) {
            neighbors.emplace_back(sensor._id);
        }
    }

    return neighbors;
}

inline std::vector<int>
neighborhood(const std::vector<Sensor> &sensors, const Sensor &p, double radius_min,
             double radius_max) {

    double lat1 = p.lat;
    double lon1 = p.lon;

    std::vector<int> neighbors = {};


    for (const auto &sensor: sensors) {

        double lat2 = sensor.lat;
        double lon2 = sensor.lon;

        double dist = distance(lat1, lon1, lat2, lon2);

        // std::cout << "Sensor ID: " << std::get<0>(sensor) << ", Distance: " << dist << std::endl;

        if (dist >= radius_min && dist <= radius_max) {
            neighbors.emplace_back(sensor._id);
        }
    }

    return neighbors;
}

inline std::vector<int>
neighborhood(const std::vector<Point> &points, const Sensor &p, double radius_min,
             double radius_max) {
    double lat1 = p.lat;
    double lon1 = p.lon;
    std::vector<int> neighbors = {};
    for (const auto &point: points) {
        double lat2 = point.lat;
        double lon2 = point.lon;

        double dist = distance(lat1, lon1, lat2, lon2);

        // std::cout << "Sensor ID: " << std::get<0>(sensor) << ", Distance: " << dist << std::endl;

        if (dist >= radius_min && dist <= radius_max) {
            neighbors.emplace_back(point._id);
        }
    }
    return neighbors;
}


#endif