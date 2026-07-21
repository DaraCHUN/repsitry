#include <iostream>
#include <cmath>
#include <complex>
#include <iomanip>

constexpr double Z0 = 376.73031355; 
constexpr double C0 = 299792458.0;   
constexpr double PI = 3.14159265358979323846;

enum Polarization { TE, TM };

class ObliquePFMOptimizer {
private:
    int N_layers;
    double f_start, f_stop;
    int f_points;
    Polarization pol;

    std::complex<double> compute_permittivity(double G_n, double f_Hz) {
        if (f_Hz <= 0.0) f_Hz = 1.0e-3;
        double f_GHz = f_Hz / 1.0e9;
        double f_term = std::pow(f_GHz, -0.4);
        double re = 1.0 + (0.039 * G_n * f_term);
        double im = -0.03 * G_n * f_term;
        std::complex<double> eps_sqrt(re, im);
        return eps_sqrt * eps_sqrt; // Eq (6): eps_m = {...}^2
    }

    double compute_objective(const double* X, double theta_deg) {
        double cost = 0.0;
        double f_step = (f_stop - f_start) / (f_points - 1.0);
        
        double theta_rad = theta_deg * PI / 180.0;
        double sin_sq_theta = std::pow(std::sin(theta_rad), 2);
        double cos_theta = std::cos(theta_rad);

        std::complex<double>* Zc = new std::complex<double>[N_layers];
        std::complex<double>* gamma = new std::complex<double>[N_layers];

        for (int step = 0; step < f_points; step++) {
            double freq = f_start + (step * f_step);
            double k0 = (2.0 * PI * freq) / C0;

            for (int i = 0; i < N_layers; i++) {
                double G = X[N_layers + i]; 
                std::complex<double> eps_m = compute_permittivity(G, freq);
                
                // Common root term for oblique incidence: sqrt(eps_m - sin^2(theta))
                std::complex<double> root_term = std::sqrt(eps_m - sin_sq_theta);
                
                // Equations (2), (3), (4), (5)
                if (pol == TE) {
                    Zc[i] = Z0 / root_term;
                } else { // TM
                    Zc[i] = (Z0 * root_term) / eps_m;
                }
                gamma[i] = std::complex<double>(0.0, 1.0) * k0 * root_term; 
            }

            // Impedance at PEC (Layer N)
            std::complex<double> Z_in = Zc[N_layers - 1] * std::tanh(gamma[N_layers - 1] * (X[N_layers - 1] / 1000.0));

            // Backward Recursion to air interface
            for (int n = N_layers - 2; n >= 0; n--) {
                double d_m = X[n] / 1000.0;
                std::complex<double> tanh_gd = std::tanh(gamma[n] * d_m);
                std::complex<double> num = Z_in + (Zc[n] * tanh_gd);
                std::complex<double> den = Zc[n] + (Z_in * tanh_gd);
                Z_in = Zc[n] * (num / den);
            }

            // Equations (7) and (8): Oblique Reflection Coefficients
            std::complex<double> S;
            if (pol == TE) {
                std::complex<double> Z0_TE = Z0 / cos_theta;
                S = (Z_in - Z0_TE) / (Z_in + Z0_TE);
            } else { // TM
                std::complex<double> Z0_TM = Z0 * cos_theta;
                S = (Z_in - Z0_TM) / (Z_in + Z0_TM);
            }
            cost += std::norm(S); 
        }

        delete[] Zc; delete[] gamma;
        return cost / f_points;
    }

    // Stabilized Newton-Raphson step for a fixed angle
    void run_newton_raphson(double* X, double theta_deg) {
        int total_vars = 2 * N_layers;
        double* X_new = new double[total_vars];
        double* grad = new double[total_vars];
        double* hess = new double[total_vars];

        double current_cost = compute_objective(X, theta_deg);
        double mu = 0.1;

        for (int iter = 0; iter < 50; iter++) {
            for (int i = 0; i < total_vars; i++) {
                double h = (i < N_layers) ? 1e-3 : 1e-2; 
                double orig = X[i];
                
                X[i] = orig + h; double cost_plus = compute_objective(X, theta_deg);
                X[i] = orig - h; double cost_minus = compute_objective(X, theta_deg);
                X[i] = orig;     
                
                grad[i] = (cost_plus - cost_minus) / (2.0 * h);
                hess[i] = (cost_plus - 2.0 * current_cost + cost_minus) / (h * h);
            }

            double eta = 1.0; 
            double new_cost = current_cost;
            bool improved = false;

            while (eta > 1e-5) {
                for (int i = 0; i < total_vars; i++) {
                    double step = eta * (grad[i] / (std::abs(hess[i]) + mu));
                    
                    double max_step = (i < N_layers) ? 0.5 : 2.0; 
                    if (step > max_step) step = max_step;
                    if (step < -max_step) step = -max_step;

                    X_new[i] = X[i] - step;

                    // Physical Constraints
                    if (i < N_layers) {
                        if (X_new[i] < 0.1) X_new[i] = 0.1;   
                        if (X_new[i] > 10.0) X_new[i] = 10.0; 
                    } else {
                        if (X_new[i] < 2.0) X_new[i] = 2.0;   
                        if (X_new[i] > 30.0) X_new[i] = 30.0; 
                    }
                }

                new_cost = compute_objective(X_new, theta_deg);
                if (new_cost < current_cost - 1e-9) { 
                    improved = true; break; 
                } else {
                    eta *= 0.5; 
                }
            }

            if (!improved || std::abs(current_cost - new_cost) < 1e-8) break;
            
            for (int i = 0; i < total_vars; i++) X[i] = X_new[i];
            current_cost = new_cost;
        }

        delete[] X_new; delete[] grad; delete[] hess;
    }

public:
    ObliquePFMOptimizer(int layers, double f_min, double f_max, int pts, Polarization p) {
        N_layers = layers;
        f_start = f_min; f_stop = f_max; f_points = pts;
        pol = p;
    }

    // Executes the PFM algorithm's 0.25-degree angular stepping sequence
    void design_oblique_absorber(double* initial_d_mm, double* initial_G, double target_theta_deg) {
        int total_vars = 2 * N_layers;
        double* X = new double[total_vars];

        // 1. Load initial Normal Incidence (0 degrees) baseline values
        for (int i = 0; i < N_layers; i++) {
            X[i] = initial_d_mm[i];
            X[N_layers + i] = initial_G[i];
        }

        std::cout << "--- PFM Oblique Incidence Design Starting ---\n";
        std::cout << "Mode: " << (pol == TE ? "TE Wave" : "TM Wave") << "\n";
        std::cout << "Target Angle: " << target_theta_deg << " degrees\n\n";

        // 2. Step angular increment by 0.25 degrees until target is met
        double current_theta = 0.0;
        double angle_step = 0.25;

        while (current_theta <= target_theta_deg) {
            std::cout << "Optimizing for Incidence Angle: " << current_theta << " deg... ";
            
            // Run Newton method. The solution from the previous angle 
            // acts as the starting condition for the current angle.
            run_newton_raphson(X, current_theta);
            
            std::cout << "Done.\n";

            if (current_theta == target_theta_deg) break;
            current_theta += angle_step;
            if (current_theta > target_theta_deg) current_theta = target_theta_deg; 
        }

        std::cout << "\n=== Final Multi-Layer Parameters at " << target_theta_deg << " Degrees ===\n";
        for (int i = 0; i < N_layers; i++) {
            std::cout << "Layer " << (i + 1) << ":"
                      << " d = " << std::fixed << std::setprecision(4) << X[i] << " m "
                      << " | G = " << X[N_layers + i] << " kg/m^3\n";
        }
        delete[] X;
    }
};

int main() {
    int N = 8, n_points;                                                        // 8-layer ultrawideband design
    float incident_angle; 
    double start_frequency, stop_frequency; 
    // Pre-calculated Normal Incidence (0 deg) baseline values [d in mm, G in kg/m^3]
    double d_initial[8] = {1.93, 2.03, 1.74, 1.50, 1.33, 1.16, 1.23, 1.89};     //
    double G_initial[8] = {17.51, 5.92, 3.10, 1.94, 1.29, 0.92, 0.65, 0.45};    // 
    //
    std::cout << "input incident angle: "       ; std::cin >> incident_angle    ; std::cout << "\n";
    std::cout << "input the start frequency: "  ; std::cin >> start_frequency   ; std::cout << "\n";
    std::cout << "input the stop frequency: "   ; std::cin >> stop_frequency    ; std::cout << "\n";
    std::cout << "input number of points: "     ; std::cin >> n_points          ; std::cout << "\n";
    // Initialize 8-layer absorber for TE Wave incidence, 2 to 18 GHz
    ObliquePFMOptimizer pfm(N, start_frequency * 1.0e9, stop_frequency * 1.0e9, n_points, TE);
    // Sweep and design for an oblique incident angle of 45 degrees
    pfm.design_oblique_absorber(d_initial, G_initial, incident_angle);

    return 0;
}