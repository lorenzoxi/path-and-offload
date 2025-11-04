#include "gurobi_c++.h"
#include <iostream>
#include <vector>
#include <tuple>
#include <map>
#include <regex>
#include <filesystem>
#include "read_data.h"
#include "neighborhood.h"
#include "SpreadingFactor.h"
#include "SF.h"
#include "LoRa.h"
#include "generate_roi_grid.h"
#include "PW.h"
#include "energy.h"

using namespace std;
namespace fs = std::filesystem;

int main() {
    //try {
    // ------------------------------------------------------------
    // DATA
    // ------------------------------------------------------------
    auto points = ::generate_roi_grid(1000,
                                      50);

    double tmp = 2.1;
    int tmp1 = (int) 3.9 + 1;
    std::cout << "tmp1: " << tmp1 << std::endl;

    vector<int> C = {1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<SpreadingFactor> F = SF125;

    const std::vector<int> vec_n_sensors = {5,50,100};

    double speed = 0.0;

    std::vector<std::pair<int, int> > start_end_points = {};
    std::vector<int> start_points = {0};

    //make starting point 0,0 => 20, 20 ajnd ending 1000,1000 => 980, 980
    for (auto p : points) {
        if (p.lat == 0 && p.lon == 0) {
            p.lat = 20;
            p.lon = 20;
        }
        if (p.lat == 1000 && p.lon == 1000) {
            p.lat = 980;
            p.lon = 980;
        }
    }
    points[0].lat = 20;
    points[0].lon = 20;
    points[points.size() - 1].lat = 980;
    points[points.size() - 1].lon = 980;

    for (auto sp: start_points) {
        int ep = find_farthest_point(points, sp);
        start_end_points.emplace_back(sp, ep);
        std::cout << "Start point: " << sp << ", End point: " << ep << std::endl;
        std::cout << "Coordinates Start point: (" << points[sp].lat << ", " << points[sp].lon << "), End point: ("
                  << points[ep].lat << ", " << points[ep].lon << ")" << std::endl;
    }

    std::vector<std::string> velocity_modes = {"variable", "fixed"};


    for (const auto& velocity : velocity_modes) {
        for (auto [start_point, end_point]: start_end_points) {
            for (auto N_SENSORS: vec_n_sensors) {

                if (velocity == "variable") {
                    switch (N_SENSORS) {
                        case 5:
                            speed = 29.997723;
                            break;
                        case 50:
                            speed = 29.461136624999998;
                            break;
                        case 100:
                            speed = 27.904765;
                            break;
                        default:
                            speed = 0.0;
                    }
                } else {
                    speed = 69;
                }
                const std::string instances_folder =
                        "C:\\Users\\perin\\CLionProjects\\2025-uav-journal\\instances\\" + std::to_string(N_SENSORS);
                std::vector<string> instances = escape_backslashes(get_txt_files_in_folder(instances_folder));

                // create solutions folder in ../solutions if it doesn't exist
                const std::string solutions_folder =
                        "C:\\\\Users\\\\perin\\\\CLionProjects\\\\2025-uav-journal\\\\solutions\\\\" + std::to_string(
                            N_SENSORS) + "_vel_"+ velocity + "_start_" + std::to_string(start_point) + "_end_" + std::to_string(end_point);

                if (!fs::exists(solutions_folder)) {
                    fs::create_directories(solutions_folder);
                }
                const int T = ((int) N_SENSORS / 8 + 1) * 2; // time horizon

                for (const auto &instance: instances) {
                    const std::vector<Sensor> S = read_sensors_from_file(instance);


                    std::string instance_name = fs::path(instance).stem().string();
                    std::cout << "======================================================================" << std::endl;
                    std::cout << "Instance: " << instance_name << ", num sensors: " << S.size() << ", T: " << T <<
                            std::endl;
                    std::cout << "Start point: " << start_point << ", End point: " << end_point << std::endl;
                    std::cout << "======================================================================" << std::endl;

                    //create a folder for the solutions of this instance if it doesn't exist
                    std::string instance_folder = solutions_folder + "\\\\" + fs::path(instance).stem().string() + "\\\\T_"
                                                  + std::to_string(T);
                    if (!fs::exists(instance_folder)) {
                        fs::create_directories(instance_folder);
                    }

                    std::vector<double> Pw_levels = {
                        0,
                        1.5849, // 2 dBm
                        3.1623, // 5 dBm
                        5.0119, // 7 dBm
                        10.0000, // 10 dBm
                        15.8489, // 12 dBm
                        25.1189, // 14 dBm
                        39.8107, // 16 dBm
                        50.1187 // 17 dBm
                    };

                    double H_cost = hovering_energy_per_sec();
                    std::cout << "Hovering energy per second: " << H_cost << " J" << std::endl;

                    double c_1 = 2.5; // coefficient for the first term in the transmission power constraint
                    double c_2 = 10; // coefficient for the second term in the transmission power constraint
                    double alpha = 27.66; // coefficient for the transmission power constraint

                    bool hasStart = any_of(points.begin(), points.end(),
                                           [&](auto &p) { return p._id == start_point; });
                    bool hasEnd = any_of(points.begin(), points.end(),
                                         [&](auto &p) { return p._id == end_point; });
                    if (!hasStart)
                        throw std::runtime_error("Start point not found: " + std::to_string(start_point));
                    if (!hasEnd) throw std::runtime_error("End point not found: " + std::to_string(end_point));

                    // Time horizon:
                    //   y: t = 0..T
                    //   x: t = 0..T-1 (transitions)
                    double bigM = 1e6;

                    // ------------------------------------------------------------
                    // GUROBI
                    // ------------------------------------------------------------
                    GRBEnv env = GRBEnv(true);
                    env.start();
                    GRBModel model = GRBModel(env);

                    // tuning (optional)
                    model.set(GRB_DoubleParam_TimeLimit, 60 * 120);
                    model.set(GRB_DoubleParam_MIPGap, 1e-3);

                    std::string log_file_name = instance_name + "_log.log";
                    std::string solution_file_name_sol = instance_name + "_sol.sol";
                    std::string solution_file_name_json = instance_name + "_sol.json";
                    std::string model_file_name = instance_name + "_model.lp";
                    model.set(GRB_StringParam_LogFile, instance_folder + "\\\\" + log_file_name);

                    // ------------------------------------------------------------
                    // VARIABLES
                    // ------------------------------------------------------------
                    // x[t,i,j]     {0,1} : move from i at time t to j at time t+1
                    // y[t,i]       {0,1} : uav is at location i at time t
                    // d[t,s,c,f]   {0,1} : at time t sensor s transmit data to uav (at i?) using channel c and spreading factor f
                    // p[i,j]       []    : power level used to reach from i to j
                    // h[t]         []    : hoovering time at time t
                    // w[t,s,c,f]   []    : linearization of p[i,j] * d[t,s,c,f]
                    // z[t,s,k]     {0,1} : select power level k for sensor s at time t

                    map<tuple<int, int, int>, GRBVar> x; // (t,i,j)
                    map<pair<int, int>, GRBVar> y; // (t,i)
                    map<tuple<int, int, int, int>, GRBVar> d; // (t,s,c,f)
                    map<tuple<int>, GRBVar> h; // (t) hoovering time at time t
                    map<tuple<int, int, int, int>, GRBVar> w;
                    map<tuple<int, int>, GRBVar> p; // (i,j) power level from i to j
                    map<tuple<int, int, int>, GRBVar> z; // (t, s, k)

                    // Create y[t,i] for t=0..T
                    for (int t = 0; t <= T; ++t) {
                        for (auto &p: points) {
                            y[{t, p._id}] = model.addVar(0.0, 1.0, 1.0, GRB_BINARY,
                                                         "y_" + to_string(t) + "_" + to_string(p._id));
                        }
                    }

                    // Create x[t,i,j] for t=0..T-1
                    for (int t = 0; t < T; ++t) {
                        for (auto &i: points)
                            for (auto &j: points) {
                                x[{t, i._id, j._id}] = model.addVar(0.0, 1.0, 1.0, GRB_BINARY,
                                                                    "x_" + to_string(t) + "_" +
                                                                    to_string(i._id) + "_" + to_string(j._id));
                            }
                    }

                    // Create d[t,s,i,c,f] variables
                    for (int t = 0; t <= T; ++t) {
                        for (const auto &s: S) {
                            for (const auto &c: C) {
                                for (const auto &f: F) {
                                    string varName = "d_" + to_string(t) + "_" + to_string(s._id) + "_" +
                                                     to_string(c) + "_" + to_string(f._id);
                                    d[{t, s._id, c, f._id}] = model.addVar(0.0, 1.0, 1.0,
                                                                           GRB_BINARY, varName);
                                }
                            }
                        }
                    }

                    // Create p[t,s] variables
                    for (int t = 0; t <= T; ++t) {
                        for (const auto &s: S) {
                            string varName = "p_" + to_string(t) + "_" + to_string(s._id);
                            p[{t, s._id}] = model.addVar(0.0, 20.0, 1.0, GRB_CONTINUOUS, varName);
                        }
                    }

                    // Create h[t] variables
                    for (int t = 0; t <= T; ++t) {
                        h[{t}] = model.addVar(0.0, 36, 1.0, GRB_INTEGER,
                                              "h_" + to_string(t));
                    }

                    // Create w[t,s,c,f] variables
                    for (int t = 0; t <= T; ++t) {
                        for (const auto &s: S) {
                            for (const auto &c: C) {
                                for (const auto &f: F) {
                                    string varName = "w_" + to_string(t) + "_" + to_string(s._id) + "_" +
                                                     to_string(c) + "_" + to_string(f._id);
                                    w[{t, s._id, c, f._id}] = model.addVar(0.0, 50.0, 1.0,
                                                                           GRB_CONTINUOUS, varName);
                                }
                            }
                        }
                    }


                    // ------------------------------------------------------------
                    // CONSTRAINTS
                    // ------------------------------------------------------------

                    // (1) Exactly one location at each time t
                    for (int t = 0; t <= T; ++t) {
                        GRBLinExpr at_t = 0;
                        for (auto &i: points) at_t += y[{t, i._id}];
                        model.addConstr(at_t == 1, "one_loc_t" + to_string(t));
                    }

                    // (2) Leave where you are: sum_j x[t,i,j] = y[t,i]
                    for (int t = 0; t < T; ++t) {
                        for (auto &i: points) {
                            GRBLinExpr out = 0;
                            for (auto &j: points) out += x[{t, i._id, j._id}];
                            model.addConstr(out == y[{t, i._id}],
                                            "leave_t" + to_string(t) + "_i" + to_string(i._id));
                        }
                    }

                    // (3) Arrive where you'll be: sum_i x[t,i,j] = y[t+1,j]
                    for (int t = 0; t < T; ++t) {
                        for (auto &j: points) {
                            GRBLinExpr in = 0;
                            for (auto &i: points) in += x[{t, i._id, j._id}];
                            model.addConstr(in == y[{t + 1, j._id}],
                                            "arrive_t" + to_string(t + 1) + "_j" + to_string(j._id));
                        }
                    }

                    // (4) Start & End
                    model.addConstr(y[{0, start_point}] == 1, "start_at_s");
                    model.addConstr(y[{T, end_point}] == 1, "end_at_e");

                    // (5) Make 'end' absorbing (can arrive early and then stay via self-loops)
                    for (int t = 0; t < T; ++t) {
                        for (auto &j: points)
                            if (j._id != end_point) {
                                model.addConstr(x[{t, end_point, j._id}] == 0, "no_leave_end_t" + to_string(t));
                            }
                        // x[{t, end, end}] allowed
                    }

                    // (6)  forbid waiting (i.e. self-loops) except at 'end'
                    for (int t = 0; t < T; ++t) {
                        for (auto &i: points)
                            if (i._id != end_point) {
                                model.addConstr(x[{t, i._id, i._id}] == 0,
                                                "no_self_loop_t" + to_string(t) + "_i" + to_string(i._id));
                            }
                    }

                    // (7) each sensor must offload once during the whole mission
                    for (const auto &s: S) {
                        GRBLinExpr offload_sum = 0;
                        for (int t = 0; t <= T; ++t) {
                            for (const auto &c: C) {
                                for (const auto &f: F) {
                                    offload_sum += d[{t, s._id, c, f._id}];
                                }
                            }
                        }
                        model.addConstr(offload_sum == 1, "sensor_offload_once_" + to_string(s._id));
                    }

                    // (8) at most 8 sensor at time can offload
                    for (int t = 0; t <= T; ++t) {
                        GRBLinExpr offload_t = 0;
                        for (const auto &s: S) {
                            for (const auto &c: C) {
                                for (const auto &f: F) {
                                    offload_t += d[{t, s._id, c, f._id}];
                                }
                            }
                        }
                        model.addConstr(offload_t <= 8, "max_8_offload_at_t_" + to_string(t));
                    }

                    // (9) at time, same channel and spreading factor can be used at most once
                    for (int t = 0; t <= T; ++t) {
                        for (const auto &c: C) {
                            for (const auto &f: F) {
                                GRBLinExpr cf_sum = 0;
                                for (const auto &s: S) {
                                    cf_sum += d[{t, s._id, c, f._id}];
                                }
                                model.addConstr(cf_sum <= 1, "cf_once_t" + to_string(t) + "_c" +
                                                             to_string(c) + "_f" + to_string(f._id));
                            }
                        }
                    }

                    // (10) toa
                    for (int t = 0; t <= T; ++t) {
                        GRBLinExpr toa_sum = 0;
                        for (const auto &s: S) {
                            for (const auto &c: C) {
                                for (const auto &f: F) {
                                    double toa = get_time_on_air(s.payload, f._id, f.bandwidth * 1000);
                                    model.addConstr(d[{t, s._id, c, f._id}] * toa <= 36,
                                                    "d_toa_lorawan_" + std::to_string(t) + "_" +
                                                    std::to_string(s._id) + "_" + std::to_string(end_point) +
                                                    "_" + std::to_string(c) + "_" + std::to_string(f._id));
                                }
                            }
                        }
                    }

                    for (int t = 0; t < T; ++t) {
                        model.addConstr(y[{t + 1, end_point}] >= y[{t, end_point}],
                                        "end_monotone_t" + to_string(t));
                    }


                    // Forbid hovering at times when you're at END
                    for (int t = 0; t <= T; ++t) {
                        model.addConstr(h[{t}] <= bigM * (1.0 - y[{t, end_point}]),
                                        "no_hover_at_end_t" + to_string(t));
                    }

                    for (int t = 0; t <= T; ++t) {
                        GRBLinExpr toa_sum = 0;
                        for (const auto &s: S) {
                            for (const auto &c: C) {
                                for (const auto &f: F) {
                                    double toa = get_time_on_air(s.payload, f._id, f.bandwidth * 1000);
                                    model.addConstr(h[{t}] >= d[{t, s._id, c, f._id}] * toa,
                                                    "hoovering_time_" + to_string(t) + "_" +
                                                    to_string(s._id) + "_" + to_string(c) + "_" +
                                                    to_string(f._id));
                                }
                            }
                        }
                    }

                    // (12) linearization of pw * d
                    for (int t = 0; t <= T; ++t) {
                        for (const auto &s: S) {
                            for (const auto &c: C) {
                                for (const auto &f: F) {
                                    // w <= maxPW * d
                                    model.addConstr(w[{t, s._id, c, f._id}] <= 50 * d[{t, s._id, c, f._id}],
                                                    "lin1_pw_d_" + to_string(t) + "_" + to_string(s._id) + "_" +
                                                    to_string(c) + "_" + to_string(f._id));

                                    // w >= minPW * d
                                    model.addConstr(w[{t, s._id, c, f._id}] >= 0 * d[{t, s._id, c, f._id}],
                                                    "lin2_pw_d_" + to_string(t) + "_" + to_string(s._id) + "_" +
                                                    to_string(c) + "_" + to_string(f._id));

                                    // w <= pw
                                    model.addConstr(w[{t, s._id, c, f._id}] <= p[{t, s._id}],
                                                    "lin3_pw_d_" + to_string(t) + "_" + to_string(s._id) + "_" +
                                                    to_string(c) + "_" + to_string(f._id));

                                    // w >= pw - (1-d)*maxPW
                                    model.addConstr(w[{t, s._id, c, f._id}] >=
                                                    p[{t, s._id}] - (1 - d[{t, s._id, c, f._id}]) * 50,
                                                    "lin4_pw_d_" + to_string(t) + "_" + to_string(s._id) + "_" +
                                                    to_string(c) + "_" + to_string(f._id));
                                }
                            }
                        }
                    }

                    // pw >= K * distance(i,j) * d[t,s,c,f] - bigM (1-d[t,s,c,f])

                    for (int t = 0; t <= T; ++t) {
                        for (const auto &s: S) {
                            for (const auto &c: C) {
                                for (const auto &f: F) {
                                    for (auto &i: points) {
                                        double dist_ij = ::euc_distance(s, i);

                                        double K = dist_ij * pow(10.0, (-c_1 * f._id + c_2 - alpha) /
                                                                       10.0);

                                        model.addConstr(p[{t, s._id}] >= K * d[{t, s._id, c, f._id}] -
                                                        bigM * (1 - d[{t, s._id, c, f._id}]),
                                                        "pw_distance_" + to_string(t) + "_" +
                                                        to_string(s._id) + "_" + to_string(c) + "_" +
                                                        to_string(f._id) + "_" + to_string(i._id));
                                    }
                                }
                            }
                        }
                    }

                    // allowed power levels
                    for (int t = 0; t <= T; ++t) {
                        for (const auto &s: S) {
                            // selection constraint: exactly one power level
                            GRBLinExpr sumZ = 0;
                            for (int k = 0; k < (int) PWs.size(); ++k) {
                                z[{t, s._id, k}] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                                                                "z_" + to_string(t) + "_" + to_string(s._id) + "_" +
                                                                to_string(k));
                                sumZ += z[{t, s._id, k}];
                            }
                            model.addConstr(sumZ == 1, "select_pw_t" + to_string(t) + "_s" + to_string(s._id));

                            // link to your existing p[t,s]
                            GRBLinExpr pw_value = 0;
                            for (int k = 0; k < (int) Pw_levels.size(); ++k) {
                                pw_value += Pw_levels[k] * z[{t, s._id, k}];
                            }
                            model.addConstr(p[{t, s._id}] == pw_value,
                                            "link_pw_t" + to_string(t) + "_s" + to_string(s._id));
                        }
                    }


                    // ------------------------------------------------------------
                    // OBJECTIVE: horizontal travel distance + hoovering time
                    // ------------------------------------------------------------
                    GRBLinExpr cost_movement = 0;
                    for (int t = 0; t < T; ++t) {
                        for (auto &i: points)
                            for (auto &j: points) {
                                double dist_ij = ::euc_distance(i, j);
                                double E_ij = calculate_energy(speed, dist_ij);
                                cost_movement += E_ij * x[{t, i._id, j._id}];
                            }
                    }
                    GRBLinExpr cost_hoovering = 0;
                    for (int t = 0; t <= T; ++t) {
                        cost_hoovering += h[{t}] * H_cost;
                    }

                    GRBLinExpr cost_transmission = 0;
                    for (int t = 0; t <= T; ++t) {
                        for (const auto &s: S) {
                            for (const auto &c: C) {
                                for (const auto &f: F) {
                                    double toa = get_time_on_air(s.payload, f._id, f.bandwidth * 1000);
                                    cost_transmission += (0.001 * w[{t, s._id, c, f._id}]) * toa;
                                }
                            }
                        }
                    }

                    //make three decision variable to track the three components of the objective
                    auto var_mov =  model.addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "cost_movement");
                    auto var_hov = model.addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "cost_hoovering");
                    auto var_trs = model.addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "cost_transmission");

                    model.addConstr(var_mov == cost_movement, "def_cost_movement");
                    model.addConstr(var_hov == cost_hoovering, "def_cost_hoovering");
                    model.addConstr(var_trs == cost_transmission, "def_cost_transmission");

                    model.setObjective(cost_movement + cost_hoovering + cost_transmission, GRB_MINIMIZE);

                    // ------------------------------------------------------------
                    // SOLVE
                    // ------------------------------------------------------------
                    model.optimize();

                    // ------------------------------------------------------------
                    // OUTPUT (time-first)
                    // ------------------------------------------------------------
                    cout << "Optimal objective value: " << model.get(GRB_DoubleAttr_ObjVal) << "\n";

                    /*
                    // Positions
                    for (int t = 0; t <= T; ++t) {
                        for (auto &p: points) {
                            if (y[{t, p._id}].get(GRB_DoubleAttr_X) > 0.5) {
                                cout << "y[" << t << "][" << p._id << "] = 1\n";
                            }
                        }
                    }

                    // Horizontal movements
                    for (int t = 0; t < T; ++t) {
                        for (auto &i: points)
                            for (auto &j: points) {
                                if (x[{t, i._id, j._id}].get(GRB_DoubleAttr_X) > 0.5) {
                                    cout << "x[" << t << "][" << i._id << "][" << j._id << "] = 1, distance = "
                                            << ::euc_distance(i, j) << " m\n";
                                }
                            }
                    }

                    // Hoovering times
                    for (int t = 0; t <= T; ++t) {
                        cout << "h[" << t << "] = " << h[{t}].get(GRB_DoubleAttr_X) << " sec\n";
                    }

                    // Offloading tasks
                    for (int t = 0; t <= T; ++t) {
                        for (const auto &s: S) {
                            for (const auto &c: C) {
                                for (const auto &f: F) {
                                    if (d[{t, s._id, c, f._id}].get(GRB_DoubleAttr_X) > 0.5) {
                                        cout << "d[" << t << "][" << s._id << "][" << c << "][" << f._id
                                                << "] = 1\n";
                                    }
                                }
                            }
                        }
                    }

                    // Transmission power levels
                    for (int t = 0; t <= T; ++t) {
                        for (const auto &s: S) {
                            auto val = p[{t, s._id}].get(GRB_DoubleAttr_X);
                            if (val > 0.1) {
                                cout << "p[" << t << "][" << s._id << "] = " << val << " dBm\n";
                            }
                        }
                    }

                    // w values
                    for (int t = 0; t <= T; ++t) {
                        for (const auto &s: S) {
                            for (const auto &c: C) {
                                for (const auto &f: F) {
                                    auto val = w[{t, s._id, c, f._id}].get(GRB_DoubleAttr_X);
                                    if (val > 0.1) {
                                        cout << "w[" << t << "][" << s._id << "][" << c << "][" << f._id
                                                << "] = " << val << " dBm\n";
                                    }
                                }
                            }
                        }
                    }
                    */

                    // save model
                    model.write(instance_folder + "\\\\" + model_file_name);
                    // save solutions
                    model.write(instance_folder + "\\\\" + solution_file_name_sol);
                    model.write(instance_folder + "\\\\" + solution_file_name_json);

                    /*
                    } catch (GRBException &e) {
                        cerr << "GRB Error code = " << e.getErrorCode() << "\n" << e.getMessage() << endl;
                        return 1;
                    } catch (const std::exception &e) {
                        cerr << "Exception: " << e.what() << endl;
                        return 1;
                    } catch (...) {
                        cerr << "Unknown exception" << endl;
                        return 1;
                    }
                     */
                }
            }
        }
    }
    return 0;
}
