#ifndef INC_005_CCNC_DRONE_ENERGY_OPT_LORA_H
#define INC_005_CCNC_DRONE_ENERGY_OPT_LORA_H

#include <cmath>
#include <iostream>

double get_time_on_air(const unsigned int &payload_bytes, const unsigned int &spreading_factor,
                       const double &bandwidth, const unsigned int &coding_rate = 1,
                       const bool &crc_on = true, const bool &header_enabled = true,
                       const unsigned int &preamble_length = 8) {

    if (spreading_factor < 7 || spreading_factor > 12) {
        std::cerr << "Invalid spreading factor. Must be between 7 and 12." << std::endl;
        return -1;
    }

    double bw = bandwidth; // in Hz
    double sf = spreading_factor;

    double cr = coding_rate + 4; // coding rate denominator (e.g., 1 -> 4/5)
    bool crc = crc_on;
    bool ih = !header_enabled;
    bool de = (sf >= 11);

    double tsym = std::pow(2, sf) / bw;
    double tpreamble = (preamble_length + 4.25) * tsym;

    double payload_symb_nb = 8 + std::max(
            std::ceil(
                    (8.0 * payload_bytes - 4.0 * sf + 28 + 16 * crc - 20 * ih) /
                    (4.0 * (sf - 2 * de))
            ) * cr,
            0.0
    );

    double tpayload = payload_symb_nb * tsym;

    double tpacket = tpreamble + tpayload;

    return tpacket;
}

#endif //INC_005_CCNC_DRONE_ENERGY_OPT_LORA_H
