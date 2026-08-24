#ifndef COUNTRATES_HH
#define COUNTRATES_HH

#include <vector>
#include <string>
#include <functional>
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <sstream>
#include <limits>
#include <algorithm>

#include "Configuration.hh"

enum class FluxType { PLAW, COMP, SEP, UNIFORM, GALACTIC, TABLE };

struct EnergyRange {
    double Emin;
    double Emax;
};

struct FluxParams {
    // PLAW / COMP
    double A = 0.0;
    double alpha = 0.0;
    double E_piv = 1.0;
    double E_peak = 1.0; // only COMP

    // SEP
    int sep_year = 0;
    int sep_order = 0;
    std::string sep_csv_path; // path to CSV with coefficients

    // Galactic
    double phiMV = 600.0;
    std::string particle = "proton";

    // Table
    std::string table_path;
};

struct RateCounts {
    int tube1Only = 0;
    int tube2Only = 0;
    int both = 0;

    [[nodiscard]] int N1() const { return tube1Only + both; }
    [[nodiscard]] int NTelescope() const { return both; }
};

struct RateResult {
    double area = 0.0;
    double integral = 0.0;     // ∫ flux(E) dE
    double Ndot = 0.0;         // area * integral
    double rate1 = 0.0;        // counter 1 singles, N1/N_gen * Ndot
    double rateTel = 0.0;      // coincidence, Ntel/N_gen * Ndot
    double rateReal1 = 0.0;    // ∫ flux(E) * A_1(E) dE
    double rateRealTel = 0.0;  // ∫ flux(E) * A_tel(E) dE
};

double fluxPLAW(double E, double A, double alpha, double E_piv);

double fluxCOMP(double E, double A, double alpha, double E_piv, double E_peak);

double fluxSEP(double E, int year, int order, const std::string& csvPath);

double fluxTable(double E, const std::string& csvPath);

double fluxUniform(double E);

double fluxGalactic(double E);

double fluxGalactic_cm2_MeV(double E_MeV, double phiMV, const std::string& name);

double J_proton(double E_GeV);

double integrateAdaptiveSimpson(const std::function<double(double)>& f,
                                double a, double b,
                                double rel_tol = 1e-6, int max_depth = 20);


RateResult computeRate(FluxType type,
                       const FluxParams& p,
                       EnergyRange eRange,
                       double A_cm2,
                       int N_histories,
                       const RateCounts& detCounts);

RateResult computeRateReal(FluxType type,
                           const FluxParams& p,
                           EnergyRange eRange,
                           const std::vector<double>& Aeff1,
                           const std::vector<double>& AeffTel,
                           int nBins);


#endif //COUNTRATES_HH
