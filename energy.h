#ifndef INC_2025_JOURNAL_DRONE_ENERGY_H
#define INC_2025_JOURNAL_DRONE_ENERGY_H

#include <cmath>

// Default UAV parameters
const double DELTA   = 1.2;
const double RHO     = 1.225;
const double S       = 0.05;
const double A       = 0.2;
const double OMEGA   = 300.0;
const double R       = 0.4;
const double K       = 0.1;
const double W       = 4.0;
const double V0      = 4.03;
const double U_TIP   = 120.0;
const double D0      = 1.0;  // not used directly

// Function to calculate energy consumption for given speed v and distance d
double calculate_energy(double v, double d,
                        double delta = ::DELTA,
                        double rho   = ::RHO,
                        double s     = ::S,
                        double A     = ::A,
                        double Omega = ::OMEGA,
                        double r     = ::R,
                        double k     = ::K,
                        double W     = ::W,
                        double v0    = ::V0,
                        double U_tip = ::U_TIP,
                        double d_0   = ::D0)
{
    // Profile power baseline
    double P0 = delta * (1.0 / 8.0) * rho * s * A * std::pow(Omega, 3) * std::pow(r, 3);

    // Induced power (approx)
    double Pi = (1.0 + k) * std::pow(W, 1.5) / std::sqrt(2.0 * rho * A);

    // Polynomial model coefficients
    double mu1 = P0;
    double mu2 = 3.0 * P0 / std::pow(U_tip, 2);
    double mu3 = Pi;
    double mu4 = 0.0; // optional linear term
    double mu5 = 0.0; // optional distance term

    // Effective propulsion power
    double energy = mu1 + mu2 * std::pow(v, 2) + mu3 + mu4 * v + mu5 * d;

    // Flight time
    double time = d / v;

    // Onboard electronics (3 W constant draw)
    double raspberry_pi_energy = 3.0 * time;

    // Total energy in Joules
    return (energy * d / v) + raspberry_pi_energy;
}


#include <cmath>

inline double hovering_energy_per_sec(
    double mass   = 4.0,      // kg (includes aircraft + battery + package if applicable)
    double g      = 9.81,     // m/s^2
    double rho    = 1.225,    // air density (kg/m^3)
    double A      = 0.503,    // rotor disk area (m^2) - TOTAL for all rotors
    double r      = 0.4,      // rotor radius (m)
    double s      = 0.05,     // rotor solidity
    double delta  = 0.012,    // profile drag coefficient
    double Omega  = 300.0,    // blade angular velocity (rad/s)
    double k      = 0.1       // induced power correction
) {
    double W = mass * g;  // total weight (N)

    // Blade profile power (P0)
    double P0 = (delta / 8.0) * rho * s * A * pow(Omega, 3) * pow(r, 3);

    // Induced power (Pi)
    double Pi = (1.0 + k) * pow(W, 1.5) / sqrt(2.0 * rho * A);

    // Total hovering power (Watts = Joules/second)
    double P_hover = P0 + Pi;

    return P_hover;  // W = J/s
}



#endif //INC_2025_JOURNAL_DRONE_ENERGY_H