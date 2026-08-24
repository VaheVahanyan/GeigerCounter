#include "Flux/TableFlux.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <regex>
#include <string>
#include <vector>

using namespace Configuration;

namespace
{
    std::vector<G4double> ExtractNumbers(const std::string& line) {
        static const std::regex re(R"(([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?))");
        std::vector<G4double> out;
        for (std::sregex_iterator it(line.begin(), line.end(), re), end; it != end; ++it) {
            out.push_back(std::stod((*it)[1].str()));
        }
        return out;
    }

    G4double LinearInterp(const G4double x, const G4double x0, const G4double y0,
                          const G4double x1, const G4double y1) {
        const G4double dx = x1 - x0;
        if (!(dx > 0.0)) return y0;
        const G4double t = (x - x0) / dx;
        return y0 + t * (y1 - y0);
    }

    G4double LogLogInterp(const G4double x, const G4double x0, const G4double y0,
                          const G4double x1, const G4double y1) {
        if (!(x > 0.0) || !(x0 > 0.0) || !(x1 > 0.0) || !(y0 > 0.0) || !(y1 > 0.0)) {
            return std::max(0.0, LinearInterp(x, x0, y0, x1, y1));
        }
        const G4double a = std::log(y1 / y0) / std::log(x1 / x0);
        return y0 * std::pow(x / x0, a);
    }

    G4double SegmentIntegral(const G4double E0, const G4double F0,
                             const G4double E1, const G4double F1) {
        if (!(E1 > E0)) return 0.0;

        const G4double f0 = std::max(0.0, F0);
        const G4double f1 = std::max(0.0, F1);

        if (f0 > 0.0 && f1 > 0.0) {
            const G4double a = std::log(f1 / f0) / std::log(E1 / E0);
            if (std::abs(a + 1.0) > 1e-12) {
                return f0 * E0 * (std::pow(E1 / E0, a + 1.0) - 1.0) / (a + 1.0);
            }
            return f0 * E0 * std::log(E1 / E0);
        }

        return 0.5 * (f0 + f1) * (E1 - E0);
    }

    G4double InvertLinearSegment(const G4double E0, const G4double F0,
                                 const G4double E1, const G4double F1,
                                 const G4double target) {
        const G4double dE = E1 - E0;
        if (!(dE > 0.0)) return E0;

        const G4double f0 = std::max(0.0, F0);
        const G4double f1 = std::max(0.0, F1);
        const G4double a = (f1 - f0) / dE;

        if (std::abs(a) < 1e-18) {
            if (f0 <= 0.0) return E0;
            const G4double x = target / f0;
            return std::clamp(E0 + x, E0, E1);
        }

        const G4double D = f0 * f0 + 2.0 * a * target;
        if (D < 0.0) return E0;

        const G4double sD = std::sqrt(D);
        const G4double x1 = (-f0 + sD) / a;
        const G4double x2 = (-f0 - sD) / a;

        const bool ok1 = (x1 >= -1e-12 && x1 <= dE + 1e-12);
        const bool ok2 = (x2 >= -1e-12 && x2 <= dE + 1e-12);

        G4double x = 0.0;
        if (ok1 && ok2) x = std::min(std::max(0.0, x1), std::max(0.0, x2));
        else if (ok1) x = x1;
        else if (ok2) x = x2;
        else x = std::clamp(target / std::max(f0, 1e-300), 0.0, dE);

        return std::clamp(E0 + x, E0, E1);
    }

    G4double InvertPowerLawSegment(const G4double E0, const G4double F0,
                                   const G4double E1, const G4double F1,
                                   const G4double target) {
        if (!(E1 > E0)) return E0;

        const G4double f0 = std::max(0.0, F0);
        const G4double f1 = std::max(0.0, F1);

        if (!(f0 > 0.0 && f1 > 0.0)) {
            return InvertLinearSegment(E0, f0, E1, f1, target);
        }

        const G4double a = std::log(f1 / f0) / std::log(E1 / E0);

        if (std::abs(a + 1.0) > 1e-12) {
            const G4double z = 1.0 + (a + 1.0) * target / (f0 * E0);
            if (z <= 0.0) return E0;
            return std::clamp(E0 * std::pow(z, 1.0 / (a + 1.0)), E0, E1);
        }

        return std::clamp(E0 * std::exp(target / (f0 * E0)), E0, E1);
    }
}

TableFlux::TableFlux() {
    configFile = "../Flux_config/Table_params.txt";
    path = GetParam(configFile, "table_path", "../TableSpectrum/flare_M2.csv");
    particle = GetParam(configFile, "particle", "proton");

    BuildCDF();
}

void TableFlux::BuildCDF() {
    EList.clear();
    FluxList.clear();
    CDF.clear();

    auto setFallback = [&] {
        EList = {1.0 * MeV, 10.0 * MeV};
        FluxList = {1.0, 1.0};
        CDF = {0.0, 9.0 * MeV};
    };

    if (path.empty()) {
        G4Exception("TableFlux::BuildCDF", "NO_PATH",
                    JustWarning, "CSV path is empty. Using trivial 2-point spectrum.");
        setFallback();
        return;
    }

    std::ifstream in(path.c_str());
    if (!in) {
        G4Exception("TableFlux::BuildCDF", "CSV_OPEN_FAIL",
                    JustWarning, ("Cannot open " + path + ", using trivial spectrum.").c_str());
        setFallback();
        return;
    }

    std::vector<Row> rows;
    rows.reserve(2048);

    std::string line;
    while (std::getline(in, line)) {
        const auto nums = ExtractNumbers(line);
        if (nums.size() < 2) continue;

        const G4double E = nums[0] * MeV;
        const G4double flux = nums[1];

        if (!std::isfinite(E) || !std::isfinite(flux)) continue;
        if (!(E > 0.0)) continue;
        if (flux < 0.0) continue;

        rows.push_back({E, flux});
    }
    in.close();

    if (rows.size() < 2) {
        G4Exception("TableFlux::BuildCDF", "CSV_NO_ROWS",
                    JustWarning, "Not enough data rows (need >=2). Using trivial spectrum.");
        setFallback();
        return;
    }

    std::sort(rows.begin(), rows.end(),
              [](const Row& a, const Row& b) {
                  return a.E_MeV < b.E_MeV;
              });

    {
        std::vector<Row> uniq;
        uniq.reserve(rows.size());
        for (const auto& r : rows) {
            if (!uniq.empty() && std::abs(r.E_MeV - uniq.back().E_MeV) <= 1e-12 * MeV) {
                uniq.back().flux = r.flux;
            } else {
                uniq.push_back(r);
            }
        }
        rows.swap(uniq);
    }

    if (rows.size() < 2) {
        G4Exception("TableFlux::BuildCDF", "CSV_NO_UNIQ_ROWS",
                    JustWarning, "Not enough unique energy rows. Using trivial spectrum.");
        setFallback();
        return;
    }

    const G4double dataEmin = rows.front().E_MeV;
    const G4double dataEmax = rows.back().E_MeV;

    G4double lo = std::max(Emin, dataEmin);
    G4double hi = std::min(Emax, dataEmax);

    if (!(hi > lo)) {
        G4Exception("TableFlux::BuildCDF", "CSV_RANGE_INVALID",
                    JustWarning, "Invalid Emin/Emax vs data range, falling back to data bounds.");
        lo = dataEmin;
        hi = dataEmax;
    }

    auto interp_at = [&](const G4double E) -> G4double {
        if (E <= rows.front().E_MeV) return rows.front().flux;
        if (E >= rows.back().E_MeV) return rows.back().flux;

        const auto it = std::upper_bound(
                                         rows.begin(), rows.end(), E,
                                         [](const G4double x, const Row& r) {
                                             return x < r.E_MeV;
                                         });

        const size_t j = static_cast<size_t>(it - rows.begin());
        const size_t i = j - 1;

        const G4double E0 = rows[i].E_MeV;
        const G4double E1 = rows[j].E_MeV;
        const G4double F0 = rows[i].flux;
        const G4double F1 = rows[j].flux;

        return LogLogInterp(E, E0, F0, E1, F1);
    };

    std::vector<G4double> Es;
    std::vector<G4double> Fs;

    Es.reserve(rows.size() + 2);
    Fs.reserve(rows.size() + 2);

    if (lo > rows.front().E_MeV) {
        Es.push_back(lo);
        Fs.push_back(std::max(0.0, interp_at(lo)));
    }

    for (const auto& r : rows) {
        if (r.E_MeV < lo || r.E_MeV > hi) continue;
        Es.push_back(r.E_MeV);
        Fs.push_back(std::max(0.0, r.flux));
    }

    if (hi < rows.back().E_MeV) {
        if (Es.empty() || std::abs(Es.back() - hi) > 1e-12 * MeV) {
            Es.push_back(hi);
            Fs.push_back(std::max(0.0, interp_at(hi)));
        }
    }

    if (Es.size() < 2) {
        G4Exception("TableFlux::BuildCDF", "CSV_RANGE_TOO_NARROW",
                    JustWarning, "Energy range too narrow. Using trivial spectrum.");
        setFallback();
        return;
    }

    EList = Es;
    FluxList = Fs;
    CDF.assign(EList.size(), 0.0);

    long double acc = 0.0L;
    for (size_t i = 0; i + 1 < EList.size(); ++i) {
        const G4double I = SegmentIntegral(EList[i], FluxList[i], EList[i + 1], FluxList[i + 1]);
        if (!std::isfinite(I) || I < 0.0) {
            G4Exception("TableFlux::BuildCDF", "BAD_SEGMENT_INTEGRAL",
                        JustWarning, "Invalid segment integral encountered. Using trivial spectrum.");
            setFallback();
            return;
        }
        acc += static_cast<long double>(I);
        CDF[i + 1] = static_cast<G4double>(acc);
    }

    if (!(acc > 0.0L) || !std::isfinite(static_cast<G4double>(acc))) {
        G4Exception("TableFlux::BuildCDF", "NONPOSITIVE_TOTAL",
                    JustWarning, "Total spectrum integral is not positive. Using trivial spectrum.");
        setFallback();
    }
}

G4double TableFlux::SampleEnergy() {
    if (EList.size() < 2 || FluxList.size() != EList.size() || CDF.size() != EList.size()) {
        return 1.0 * MeV;
    }

    const G4double total = CDF.back();
    if (!(total > 0.0) || !std::isfinite(total)) {
        return EList.front();
    }

    const G4double target = G4UniformRand() * total;

    const auto it = std::upper_bound(CDF.begin(), CDF.end(), target);
    if (it == CDF.begin()) return EList.front();
    if (it == CDF.end()) return EList.back();

    const size_t j = static_cast<size_t>(std::distance(CDF.begin(), it));
    const size_t i = j - 1;

    const G4double local = target - CDF[i];
    return InvertPowerLawSegment(EList[i], FluxList[i], EList[j], FluxList[j], local);
}
