#ifndef INC_005_CCNC_DRONE_ENERGY_OPT_SPREADINGFACTOR_H
#define INC_005_CCNC_DRONE_ENERGY_OPT_SPREADINGFACTOR_H

struct SpreadingFactor {
    unsigned int _id;   // Spreading factor value
    double bandwidth;   // Bandwidth in kHz
    double Kbps;         // Bits per second
    double range_min;   // Minimum range in km
    double range_max;   // Maximum range in km

    SpreadingFactor(unsigned int id, double bw, double Kbps, double r_min, double r_max)
            : _id(id), bandwidth(bw), Kbps(Kbps), range_min(r_min), range_max(r_max) {}
};


#endif