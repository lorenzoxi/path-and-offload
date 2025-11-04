#ifndef READ_DATA_H // Ensure you have include guards for read_data.h as well
#define READ_DATA_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <sstream>
#include <utility>
#include <iostream>
#include <vector>
#include <string>
#include <filesystem>

namespace fs = std::filesystem;
// 1. Include nlohmann/json.hpp FIRST
#include "include/json.hpp"

// 2. Then, define your 'json' alias
using json = nlohmann::json;

// 3. Now, include your custom types that have from_json/to_json defined
//    (Point.h and Sensor.h should ideally already contain their from_json/to_json functions)
#include "Point.h"
#include "Sensor.h" // Assuming Sensor also has its own from_json/to_json or isn't deserialized directly with get<std::vector<Sensor>>()

// Function to load the ROI and generate a discretization grid as vector<int>
inline std::vector<Point> load_roi_discretization(const std::string &path) {
    // Open and parse the JSON file
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open JSON file: " + path);
    }

    json j;
    file >> j;

    // Extract ROI attributes
    if (!j.contains("roi")) {
        throw std::runtime_error("JSON does not contain 'roi' field.");
    }

    // This duplicate check is unnecessary:
    // if (!j.contains("roi")) {
    //     throw std::runtime_error("JSON does not contain 'roi' field.");
    // }

    auto roi = j["roi"];

    // Check if the 'roi' object contains the 'points' array
    if (!roi.contains("points") || !roi["points"].is_array()) {
        throw std::runtime_error("'roi' field does not contain a 'points' array.");
    }

    std::vector<Point> points;
    try {
        // Iterate over each element in the "points" JSON array
        for (const auto &point_json: roi["points"]) {
            Point p(0, 0, 0, 0); // Create a default Point object
            // Manually call from_json or access members
            // Option 1: Call your existing from_json function if it's visible (it should be)
            from_json(point_json, p);
            points.push_back(p);

        }

    } catch (const json::exception &e) {
        throw std::runtime_error("Failed to parse 'points' array: " + std::string(e.what()));
    }


    if (points.empty()) {
        throw std::runtime_error("No points found in the 'roi' field.");
    }


    return points;
}


std::vector<Sensor> load_sensors(const std::string &path, const double &default_height = 0.0) {
    // Open and parse the JSON file
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open JSON file: " + path);
    }

    json j;
    file >> j;

    // Check if "sensors" field exists
    if (!j.contains("sensors") || !j["sensors"].is_array()) {
        throw std::runtime_error("'sensors' field is missing or not an array.");
    }

    std::vector<Sensor> sensors;

    // Iterate over each sensor entry
    for (const auto &item: j["sensors"]) {

        unsigned int id = item.value("_id", 0);
        double lat = item.value("lat", 0.0);
        double lon = item.value("lon", 0.0);
        double offload = item.value("offload", 0.0);

        sensors.emplace_back(id, lat, lon, default_height,offload);
    }

    return sensors;
}


inline std::map<std::pair<int, int>, double> load_energy_map(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    json j;
    file >> j;

    std::map<std::pair<int, int>, double> energy_map;

    for (auto it = j.begin(); it != j.end(); ++it) {
        std::string key = it.key();  // format "id1_id2"
        double value = it.value();

        std::istringstream iss(key);
        std::string token;
        int id1, id2;

        if (std::getline(iss, token, '_')) {
            id1 = std::stoi(token);
        }
        if (std::getline(iss, token, '_')) {
            id2 = std::stoi(token);
        }

        energy_map[{id1, id2}] = value;
    }

    return energy_map;
}


inline std::vector<Sensor> read_sensors_from_file(const std::string& filename) {
    std::vector<Sensor> sensors;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: could not open file " << filename << std::endl;
        return sensors;
    }

    double lat, lon, height;
    double payload = 30.0;
    unsigned int id = 0;

    while (file >> lat >> lon >> height) {
        sensors.emplace_back(id++, lat, lon, height, payload);
    }

    file.close();
    return sensors;
}


inline std::vector<std::string> get_txt_files_in_folder(const std::string& folderPath) {
    std::vector<std::string> txtFiles;

    try {
        for (const auto& entry : fs::directory_iterator(folderPath)) {
            if (entry.is_regular_file() && entry.path().extension() == ".txt") {
                txtFiles.push_back(entry.path().string());
            }
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
    }

    return txtFiles;
}

inline std::vector<std::string> escape_backslashes(const std::vector<std::string>& input) {
    std::vector<std::string> output;
    output.reserve(input.size());

    for (const auto& s : input) {
        std::string escaped;
        escaped.reserve(s.size()); // avoid reallocations

        for (char c : s) {
            if (c == '\\') {
                escaped += "\\\\"; // double backslash
            } else {
                escaped += c;
            }
        }
        output.push_back(std::move(escaped));
    }
    return output;
}



#endif // READ_DATA_H