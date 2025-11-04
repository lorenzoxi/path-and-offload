#ifndef INC_005_CCNC_DRONE_ENERGY_OPT_SENSOR_H
#define INC_005_CCNC_DRONE_ENERGY_OPT_SENSOR_H

#include "include/json.hpp"

using json = nlohmann::json;

struct Sensor {
    unsigned int _id;    // Unique identifier for the sensor
    double lat;
    double lon;
    double z;
    double payload;   // payload value (e.g., processing offload percentage)

    Sensor(unsigned int _id, double lat, double lon, double z, double payload = 0.0)
            : _id(_id), lat(lat), lon(lon), z(z), payload(payload) {}
};

inline void from_json(const nlohmann::json &j, Sensor &s) {
    s._id = j.at("_id").get<unsigned int>();
    s.lat = j.at("lat").get<double>();
    s.lon = j.at("lon").get<double>();
    s.z = j.at("z").get<double>();
    s.payload = j.value("payload", 0.0);
}

inline void to_json(nlohmann::json &j, const Sensor &s) {
    j = nlohmann::json{
            {"_id",        s._id},
            {"lat",        s.lat},
            {"y",          s.lon},
            {"z",          s.z},
            {"payload", s.payload}
    };
}


#endif
