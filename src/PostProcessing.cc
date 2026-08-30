#include "PostProcessing.hh"

namespace fs = std::filesystem;
using namespace Configuration;

PostProcessing::PostProcessing(std::string outputFolderName,
                               std::string particle,
                               const double norm_cm2,
                               const bool isotropic) : outputFolderName(std::move(outputFolderName)),
                                                       eMinMeV(Emin),
                                                       eMaxMeV(Emax),
                                                       particleName(std::move(particle)),
                                                       norm(norm_cm2),
                                                       isotropic(isotropic) {
    ROOT::EnableImplicitMT();
    gROOT->SetBatch(kTRUE);

    gErrorIgnoreLevel = kWarning;

    OpenRootFile();
    PrepareOutputDirs();
}

PostProcessing::~PostProcessing() = default;

void PostProcessing::OpenRootFile() {
    rootFile.reset(TFile::Open(outputFile.c_str(), "READ"));
    if (!rootFile || rootFile->IsZombie()) {
        throw std::runtime_error("Failed to open ROOT file: " + outputFile);
    }
}

void PostProcessing::PrepareOutputDirs() {
    fs::path rootPath(outputFile.data());
    fs::path rootDir = rootPath.parent_path();

    postProcessingDir = (rootDir / "post_processing").string();
    runDir = (fs::path(postProcessingDir) / outputFolderName).string();

    csvDir = (fs::path(runDir) / "CSV").string();
    histogramsDir = (fs::path(runDir) / "Histograms").string();

    fs::create_directories(csvDir);
    fs::create_directories(histogramsDir);
}

TH1* PostProcessing::GetHist(const std::string& histName) {
    TH1* h = nullptr;
    rootFile->GetObject(histName.c_str(), h);
    return h;
}

TH1* PostProcessing::GetHistOrThrow(const std::string& histName) {
    TH1* h = GetHist(histName);
    if (!h) {
        throw std::runtime_error("Histogram not found: " + histName);
    }
    return h;
}

std::string PostProcessing::AreaUnit() const {
    return isotropic ? "cm^2*sr" : "cm^2";
}

std::string PostProcessing::AreaUnitRoot() const {
    return isotropic ? "cm^{2} #upoint sr" : "cm^{2}";
}

double PostProcessing::BinCentre(const double eLow, const double eHigh) {
    if (isLogBin) {
        if (eLow <= 0.0 || eHigh <= 0.0) return 0.0;
        return std::sqrt(eLow * eHigh);
    }
    return 0.5 * (eLow + eHigh);
}

double PostProcessing::EffAreaErrFromCounts(const double n0, const double n, const double effArea) {
    if (n0 <= 0.0 || n <= 0.0) return 0.0;
    const double val = (n0 - n) / (n * n0);
    if (val <= 0.0) return 0.0;
    return effArea * std::sqrt(val);
}

static short GetColorForParticle(const std::string& particleName) {
    if (particleName == "gamma") return kGreen + 2;
    if (particleName == "e-" || particleName == "e+") return kOrange + 7;
    if (particleName == "proton") return kBlue + 1;
    if (particleName == "mu-" || particleName == "mu+") return kMagenta + 1;
    if (particleName == "neutron") return kCyan + 2;
    if (particleName == "alpha") return kRed + 1;
    return kBlack;
}

void PostProcessing::SaveHistPng(const std::string& histName,
                                 const std::string& outPngPath,
                                 const std::string& plotTitle,
                                 const std::string& yTitle,
                                 const bool filled,
                                 const bool log_y,
                                 const bool weighted) {
    TH1* h = GetHist(histName);
    if (!h) return;

    TCanvas canvas("canvas", "canvas", 1920, 1080);
    canvas.SetLogx(isLogBin);
    canvas.SetLogy(log_y and isLogBin);

    const short color = GetColorForParticle(particleName);

    h->SetTitle(plotTitle.c_str());
    h->GetXaxis()->SetTitle("Energy [MeV]");
    h->GetYaxis()->SetTitle(yTitle.c_str());

    if (eMinMeV > eThreshold && eMaxMeV > eMinMeV) {
        h->GetXaxis()->SetRangeUser(eMinMeV, eMaxMeV);
    } else {
        h->GetXaxis()->SetRangeUser(eThreshold, eMaxMeV);
    }

    h->SetLineColor(color);
    h->SetLineWidth(filled ? 1 : 2);

    if (filled) {
        h->SetFillColorAlpha(color, 0.35);
        h->SetFillStyle(1001);
        h->SetMarkerColor(color);
    } else {
        h->SetFillStyle(0);
    }

    if (fluxType != "Uniform" and weighted) h->Scale(1.0, "width");
    h->Draw("HIST");

    canvas.SaveAs(outPngPath.c_str());
}


void PostProcessing::ExportTreeToCsv(const std::string& treeName,
                                     const std::string& csvPath) {
    TTree* tree = nullptr;
    rootFile->GetObject(treeName.c_str(), tree);
    if (!tree) {
        throw std::runtime_error("TTree/NTuple not found: " + treeName);
    }

    auto* leaves = tree->GetListOfLeaves();
    if (!leaves || leaves->GetEntries() == 0) {
        throw std::runtime_error("No leaves found in tree: " + treeName);
    }

    std::ofstream out(csvPath);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output CSV: " + csvPath);
    }

    const int nLeaves = leaves->GetEntries();
    for (int i = 0; i < nLeaves; ++i) {
        auto* leaf = dynamic_cast<TLeaf*>(leaves->At(i));
        out << (leaf ? leaf->GetName() : "unknown");
        if (i + 1 != nLeaves) out << ",";
    }
    out << "\n";

    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        for (int i = 0; i < nLeaves; ++i) {
            auto* leaf = dynamic_cast<TLeaf*>(leaves->At(i));
            if (!leaf) {
                out << "";
                if (i + 1 != nLeaves) out << ",";
                continue;
            }

            if (leaf->InheritsFrom(TLeafC::Class())) {
                auto* lc = dynamic_cast<TLeafC*>(leaf);
                std::string s = lc->GetValueString();

                const bool needQuotes = s.find(',') != std::string::npos ||
                    s.find('"') != std::string::npos ||
                    s.find('\n') != std::string::npos;
                if (needQuotes) {
                    std::string esc;
                    esc.reserve(s.size() + 8);
                    for (char c : s) {
                        if (c == '"') esc += "\"\"";
                        else esc += c;
                    }
                    out << "\"" << esc << "\"";
                } else {
                    out << s;
                }
            } else {
                const int nData = leaf->GetNdata();
                if (nData <= 1) {
                    out << std::setprecision(17) << leaf->GetValue(0);
                } else {
                    std::ostringstream cell;
                    cell << std::setprecision(17);
                    for (int k = 0; k < nData; ++k) {
                        if (k) cell << ";";
                        cell << leaf->GetValue(k);
                    }
                    out << "\"" << cell.str() << "\"";
                }
            }

            if (i + 1 != nLeaves) out << ",";
        }
        out << "\n";
    }

    out.close();
}

void PostProcessing::ExtractNtData() {
    ExportTreeToCsv("edep", (fs::path(csvDir) / "edep.csv").string());
    ExportTreeToCsv("primary", (fs::path(csvDir) / "primary.csv").string());

    if (saveSecondaries) {
        ExportTreeToCsv("event", (fs::path(csvDir) / "event.csv").string());
        ExportTreeToCsv("interactions", (fs::path(csvDir) / "interactions.csv").string());
    }
}

void PostProcessing::SaveResponse() {
    TH1* gen = GetHistOrThrow("genEnergyHist");
    TH1* trig1 = GetHistOrThrow("trigEnergyHist");
    TH1* trigTel = GetHistOrThrow("trigEnergyHistTel");

    if (gen->GetNbinsX() != nBins || trig1->GetNbinsX() != nBins || trigTel->GetNbinsX() != nBins) {
        throw std::runtime_error("Histogram binning mismatch among genEnergyHist/trigEnergyHist/trigEnergyHistTel");
    }

    const std::string responseHist = isotropic ? "sensitivityHist" : "effAreaHist";
    const std::string responseName = isotropic ? "Geometric factor" : "Effective area";

    SaveHistPng("genEnergyHist", (fs::path(histogramsDir) / "genEnergyHist.png").string(),
                "N_{gen} vs Energy", "Counts", true, true, true);

    SaveHistPng("trigEnergyHist", (fs::path(histogramsDir) / "trigEnergyHist.png").string(),
                "N_{trig} vs Energy, counter 1", "Counts", true, true, true);

    SaveHistPng("trigEnergyHistTel", (fs::path(histogramsDir) / "trigEnergyHistTel.png").string(),
                "N_{trig} vs Energy, telescope", "Counts", true, true, true);

    SaveHistPng(responseHist, (fs::path(histogramsDir) / "response_counter1.png").string(),
                responseName + " vs Energy, counter 1", responseName + " [" + AreaUnitRoot() + "]", false);

    SaveHistPng(responseHist + "Tel", (fs::path(histogramsDir) / "response_telescope.png").string(),
                responseName + " vs Energy, telescope", responseName + " [" + AreaUnitRoot() + "]", false);

    const std::string csvName = "A_eff_" + particleName + "_" + geometryType + "_"
                                + GenSurface::DirectionTag() + ".csv";
    const std::string outPath = (fs::path(csvDir) / csvName).string();

    std::ofstream out(outPath);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output CSV: " + outPath);
    }

    out << "# particle: " << particleName << "\n";
    out << "# geometry: " << geometryType << "\n";
    out << "# flux_direction: " << fluxDirection << "\n";
    out << "# theta_deg: " << beamTheta / deg << "\n";
    out << "# phi_deg: " << beamPhi / deg << "\n";
    out << "# normalisation: " << norm << " " << AreaUnit() << "\n";
    out << "# energy unit: MeV, area unit: " << AreaUnit() << "\n";
    out << "E_low,E_high,E_centre,N_gen,N_trig_1,N_trig_tel,A_1,dA_1,A_tel,dA_tel\n";
    out << std::setprecision(17);

    auto* ax = gen->GetXaxis();

    for (int i = 1; i <= nBins; ++i) {
        const double eLow = ax->GetBinLowEdge(i);
        const double eHigh = ax->GetBinUpEdge(i);
        const double eCentre = BinCentre(eLow, eHigh);

        const double n0 = gen->GetBinContent(i);
        const double n1 = trig1->GetBinContent(i);
        const double nTel = trigTel->GetBinContent(i);

        const double a1 = n0 > 0.0 ? norm * n1 / n0 : 0.0;
        const double aTel = n0 > 0.0 ? norm * nTel / n0 : 0.0;

        out << eLow << ","
            << eHigh << ","
            << eCentre << ","
            << n0 << ","
            << n1 << ","
            << nTel << ","
            << a1 << ","
            << EffAreaErrFromCounts(n0, n1, a1) << ","
            << aTel << ","
            << EffAreaErrFromCounts(n0, nTel, aTel) << "\n";
    }

    out.close();
}
