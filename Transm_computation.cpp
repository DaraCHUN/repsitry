#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>



class PlanarAbsorberSolver {
private:
    // Physical constants:
    const double PI = acos(-1.0), ep0 = 1.0e-9 / (36.0 * PI), mu0 = 4.0 * PI * 1.0e-7, C0 = 1.0 / std::sqrt(ep0 * mu0), Z0 = std::sqrt(mu0 / ep0); 
public:
    // FIX 1 & 2: Computes sqrt(eps_r^n) directly with frequency scaled to GHz
    std::complex<double> compute_sqrt_permittivity(double G_n, double f_Hz) {
        if (f_Hz <= 0.0) f_Hz = 1.0e-3; // Guard against division by zero

        // Convert raw Hz to GHz for empirical power law f^(-0.4)
        double f_GHz = f_Hz / 1.0e9;
        double f_term = std::pow(f_GHz, -0.4);

        double re = 1.0 + (0.039 * G_n * f_term);
        double im = -0.03 * G_n * f_term;
        // Directly return sqrt(eps_r^n) to avoid branch-cut phase errors
        return std::complex<double>(re, im);
    }

    // Solves wave impedance cascade across N layers
    std::complex<double> compute_reflection_coefficient(
        const double* G, 
        const double* d, 
        int N, 
        double freq_Hz
    ) {
        if (N <= 0) return std::complex<double>(0.0, 0.0);

        double k0 = (2.0 * PI * freq_Hz) / C0;

        std::complex<double>* Zc = new std::complex<double>[N];
        std::complex<double>* gamma = new std::complex<double>[N];

        // Stage 1: Equations (29) & (30)
        for (int i = 0; i < N; i++) {
            std::complex<double> sqrt_eps = compute_sqrt_permittivity(G[i], freq_Hz);

            Zc[i] = Z0 / sqrt_eps;                                   // Eq (29)
            gamma[i] = std::complex<double>(0.0, 1.0) * k0 * sqrt_eps; // Eq (30)
        }

        // Stage 2: Equation (32) - PEC Load Condition at Layer N
        std::complex<double> Z_in = Zc[N - 1] * std::tanh(gamma[N - 1] * d[N - 1]);
        
        // Stage 3: Equation (31) - Backward Recursion (Layer N-1 down to 1)
        for (int n = N - 2; n >= 0; n--) {
            std::complex<double> tanh_gd = std::tanh(gamma[n] * d[n]);

            std::complex<double> num = Z_in + (Zc[n] * tanh_gd);
            std::complex<double> den = Zc[n] + (Z_in * tanh_gd);

            Z_in = Zc[n] * (num / den);
        }
        // Stage 4: Equation (33) - Total Reflection Coefficient
        std::complex<double> Gamma_TL = (Z_in - Z0) / (Z_in + Z0);

        // Dynamic memory cleanup
        delete[] Zc;
        delete[] gamma;

        return Gamma_TL;
    }
    // Broadband Sweep Engine
    void run_broadband_sweep(
        const double* G,
        const double* d,
        int N,
        double f_start_Hz,
        double f_stop_Hz,
        int num_pts,
        const char* output_filename
    ) {

        if (num_pts < 2) num_pts = 2;

        std::ofstream outFile(output_filename);
        if (!outFile.is_open()) {
            std::cerr << "Error opening " << output_filename << "\n";
            return;
        }

        double f_step = (f_stop_Hz - f_start_Hz) / (num_pts - 1.0);

        std::cout << "--- Sweeping Broadband Spectrum (" 
                  << f_start_Hz / 1.0e9 << " GHz to " 
                  << f_stop_Hz / 1.0e9 << " GHz) ---\n\n";

        for (int i = 0; i < num_pts; i++) {
            double current_freq = f_start_Hz + (i * f_step);
            
            std::complex<double> Gamma = compute_reflection_coefficient(G, d, N, current_freq);
            
            double mag_Gamma = std::abs(Gamma);
            double reflection_loss_dB = 20.0 * std::log10(mag_Gamma);

            // Output 2-column Cartesian format: "Frequency(GHz) ReflectionLoss(dB)"
            outFile << (current_freq / 1.0e9) << " " << reflection_loss_dB << "\n";
        }

        outFile.close();
        std::cout << "Successfully written broadband simulation to " << output_filename << "\n\n";
    }
};

/* Compute */
int main() {
    int N ;             // Layer absorber
    int num_freq_pts;   //broadband evaluation points
    
    // Layer parameters (0-indexed: Index 0 = Layer 1, Index 1 = Layer 2)

    double* G = new double[]{0.0}; 
    double* d = new double[]{0.0}; 
    //double G = 12.1906; // Carbon content values
    //double d = 0.1189379;       // Layer thicknesses (1.5mm & 1.0mm)
do{
    std::cout << "===========================================\n";
    std::cout << "  Multi-Layer EM Absorber Broadband Solver \n";
    std::cout << "===========================================\n";
    std::cout << "Enter the number of absorber layers (N): ";
    std::cin >> N;
//
    if(N <= 0){
        std::cout << "Invalid layer count. Defaulting to N = 1!\n";
        N = 1;
// 
    }else if(N == 1){
        d = new double[]{0.1265476};                                            // one Layer Carbon Content values
        G = new double[]{12.4096};                                              // one Layer thickness
    }else if(N == 2){
        G = new double[]{4.597653, 24.683648};                                  // two layers Carbon content values
        d = new double[]{0.1224987, 0.0730035};                                 // two layers thicknesses 
    }else if(N == 3){
        G = new double[]{3.6556, 9.4782, 27.6146};                              // three Layers Carbon content values
        d = new double[]{0.0734501, 0.0606653, 0.0557619};                      // three Layers thicknesses 
    }else if(N==4){
       G = new double[]{3.1132, 5.5068, 10.8580, 28.3232};                      // four layers Carbon content values
       d = new double[]{0.0494807, 0.0429712, 0.0475493, 0.0530626};            // four Layers thicknesses (1.5mm & 1.0mm)
    }else if(N==5){
       G = new double[]{2.0373, 3.3446, 5.6553, 10.8990, 28.3334};              // five Layers Carbon content values
       d = new double[]{0.0145521, 0.0411875, 0.0401926, 0.0469192, 0.053018};  // five Layers thicknesses (1.5mm & 1.0mm) 
    }else{ 
       std::cerr <<"The layer must be 1 - 5!\n";                                // return message to the console
    }
}while(N > 5);

    // Broadband Sweep Range: input frequency (S to Ku band)
    double f_start, f_stop  ; // 
    std::cout<<"Input minimum frequency: "; std::cin >> f_start; std::cout      << "\n";
    std::cout<<"Input maximum frequency: "; std::cin >> f_stop; std::cout       << "\n";
    std::cout<<"Input frequency points : "; std::cin >> num_freq_pts; std::cout << "\n";

    PlanarAbsorberSolver solver;

    // Execute broadband computation and save results for Gnuplot
    solver.run_broadband_sweep(G, d, N, f_start * 1e9, f_stop * 1e9, num_freq_pts, "broadband_reflection_loss.txt");
    // Dynamic memory cleanup
    delete[] G;
    delete[] d;
    return 0;
}