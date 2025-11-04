#ifndef INC_2025_JOURNAL_DRONE_GENERATE_ROI_GRID_H
#define INC_2025_JOURNAL_DRONE_GENERATE_ROI_GRID_H

#include "Point.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <limits>

inline std::vector<Point>
generate_roi_grid(const int &grid_side_meters, const int &spacing_meters) {
    std::vector<Point> grid_points;
    int id_counter = 0;

    // Calculate the number of points along each axis
    int num_points_x = grid_side_meters / spacing_meters;
    int num_points_y = grid_side_meters / spacing_meters;

    // Generate grid points
    for (int i = 0; i <= num_points_y; ++i) {
        for (int j = 0; j <= num_points_x; ++j) {
            double lat = 0 + i * spacing_meters;
            double lon = 0 + j * spacing_meters;
            Point p(j * spacing_meters, i * spacing_meters, 0, id_counter, lat, lon);
            grid_points.push_back(p);
            id_counter++;
        }
    }

    // Print the generated points
    std::cout << "Generated " << grid_points.size() << " grid points:" << std::endl;
    std::cout << "num_points_x: " << num_points_x << ", num_points_y: " << num_points_y
              << std::endl;
    for (const auto &point: grid_points) {
        std::cout << "Point ID: " << point._id << ", Lat: " << point.lat << ", Lon: " << point.lon
                  << std::endl;
    }

    return grid_points;
}

inline double points_distance(const Point &a, const Point &b) {
    double dx = static_cast<double>(a.x - b.x);
    double dy = static_cast<double>(a.y - b.y);
    return std::sqrt(dx * dx + dy * dy);
}

inline int find_farthest_point(const std::vector<Point> &grid_points, int start_id) {
    // Find the reference point
    bool found = false;
    Point start_point(0, 0, 0, -1);
    for (const auto &p : grid_points) {
        if (p._id == start_id) {
            start_point = p;
            found = true;
            break;
        }
    }

    if (!found) {
        throw std::runtime_error("Point with given ID not found in grid.");
    }

    // Find the farthest point
    Point farthest_point = start_point;
    double max_dist = -1e6;

    for (const auto &p : grid_points) {
        if (p._id == start_point._id) continue; // skip itself

        double d = points_distance(start_point, p);
        if (d > max_dist) {
            max_dist = d;
            farthest_point = p;
        }
    }

    if (farthest_point._id == start_point._id) {
        throw std::runtime_error("No farthest point found (grid may contain only one point).");
    }
    std::cout << "Farthest point from ID " << start_id << " is ID " << farthest_point._id
              << " at distance " << max_dist << " meters." << std::endl;
    return farthest_point._id;
}

#endif //INC_2025_JOURNAL_DRONE_GENERATE_ROI_GRID_H
