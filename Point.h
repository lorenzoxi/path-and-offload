#ifndef INC_005_CCNC_DRONE_ENERGY_OPT_POINT_H
#define INC_005_CCNC_DRONE_ENERGY_OPT_POINT_H

#include <tuple>
#include "include/json.hpp"

using json = nlohmann::json;

struct Point {
    int x; // x-coordinate
    int y; // y-coordinate
    int z; // z-coordinate (optional, can be used for altitude or other purposes)
    int _id; // Unique _identifier for the point
    double lat; // Latitude (optional, can be used for geographical coordinates)
    double lon; // Longitude (optional, can be used for geographical coordinates)

    Point(int x, int y, int z, int _id, double lat = 0.0, double lon = 0.0)
            : x(x), y(y), z(z), _id(_id), lat(lat), lon(lon) {}

};

inline void from_json(const nlohmann::json &j, Point &p) {
    p.x = j.at("x").get<int>();
    p.y = j.at("y").get<int>();
    p.z = j.value("z", 0);
    p._id = j.at("_id").get<int>();
    p.lat = j.value("lat", 0.0);
    p.lon = j.value("lon", 0.0);
}

inline void to_json(nlohmann::json &j, const Point &p) {
    j = nlohmann::json{{"x",   p.x},
                       {"y",   p.y},
                       {"z",   p.z},
                       {"_id", p._id},
                       {"lat", p.lat},
                       {"lon", p.lon}};
}


#endif //INC_005_CCNC_DRONE_ENERGY_OPT_POINT_H
