#include "RunAction.hh"

using namespace Configuration;

RunAction::RunAction() {
    const G4String fileName = "GeigerCounter.root";
    analysisManager = new AnalysisManager(fileName);

    auto* mgr = G4AccumulableManager::Instance();
    mgr->Register(tube1Only);
    mgr->Register(tube2Only);
    mgr->Register(bothTubs);
}

RunAction::RunAction(const double Agen_cm2) : EminMeV(Emin),
                                              EmaxMeV(Emax),
                                              area(Agen_cm2) {
    analysisManager = new AnalysisManager(outputFile, nBins);
    if (nBins < 1) {
        throw std::runtime_error("RunAction: nbins must be >= 1");
    }
    if (EminMeV <= 0.0 || EmaxMeV <= 0.0 || EmaxMeV < EminMeV) {
        throw std::runtime_error("RunAction: invalid Emin/Emax");
    }
    G4double E_low = EmaxMeV > EminMeV ? EminMeV : eThreshold;
    if (isLogBin) {
        logEmin = std::log10(E_low);
        logEmax = std::log10(EmaxMeV);
        invDlogE = static_cast<double>(nBins) / (logEmax - logEmin);
    } else {
        invDlinearE = static_cast<double>(nBins) / (EmaxMeV - E_low);
    }

    effArea.assign(nBins, 0.0);

    BookAccumulables();
}

void RunAction::BookAccumulables() {
    auto* mgr = G4AccumulableManager::Instance();

    mgr->Register(tube1Only);
    mgr->Register(tube2Only);
    mgr->Register(bothTubs);

    genCounts.clear();
    trigCounts.clear();
    genCounts.reserve(EminMeV < EmaxMeV ? nBins : 1);
    trigCounts.reserve(nBins);

    for (int i = 0; i < nBins; ++i) {
        trigCounts.emplace_back(0.0);
        mgr->Register(trigCounts.back());
    }
    for (int i = 0; i < (EminMeV < EmaxMeV ? nBins : 1); ++i) {
        genCounts.emplace_back(0.0);
        mgr->Register(genCounts.back());
    }
}

RunAction::~RunAction() {
    delete analysisManager;
}

void RunAction::BeginOfRunAction(const G4Run*) {
    analysisManager->Open();
    auto* mgr = G4AccumulableManager::Instance();
    mgr->Reset();

    totals = {};
    std::fill(effArea.begin(), effArea.end(), 0.0);
}

void RunAction::EndOfRunAction(const G4Run*) {
    auto* mgr = G4AccumulableManager::Instance();
    mgr->Merge();
    if (G4Threading::IsMasterThread()) {
        totals.tube_1 = tube1Only.GetValue();
        totals.tube_2 = tube2Only.GetValue();
        totals.both = bothTubs.GetValue();
        if (EminMeV < EmaxMeV) {
            FillDerivedHists();
        }
    }

    analysisManager->Close();
}

int RunAction::FindBin(double E_MeV) const {
    const double E_low = EmaxMeV > EminMeV ? EminMeV : eThreshold;
    if (E_MeV < E_low || E_MeV >= EmaxMeV) return -1;

    int idx = 0;
    if (isLogBin) {
        const double le = std::log10(E_MeV);
        idx = static_cast<int>((le - logEmin) * invDlogE);
    } else {
        idx = static_cast<int>((E_MeV - E_low) * invDlinearE);
    }

    if (idx < 0) idx = 0;
    if (idx >= nBins) idx = nBins - 1;
    return idx;
}

double RunAction::BinCenterMeV(int i) const {
    const double E_low = EmaxMeV > EminMeV ? EminMeV : eThreshold;

    if (isLogBin) {
        const double ratio = EmaxMeV / E_low;
        const double e1 = E_low * std::pow(ratio, static_cast<double>(i) / nBins);
        const double e2 = E_low * std::pow(ratio, static_cast<double>(i + 1) / nBins);
        return std::sqrt(e1 * e2);
    }
    const double dE = (EmaxMeV - E_low) / nBins;
    return E_low + dE * (static_cast<double>(i) + 0.5);
}

double RunAction::BinWidthMeV(int i) const {
    const double E_low = EmaxMeV > EminMeV ? EminMeV : eThreshold;

    if (isLogBin) {
        const double ratio = EmaxMeV / E_low;
        const double e1 = E_low * std::pow(ratio, static_cast<double>(i) / nBins);
        const double e2 = E_low * std::pow(ratio, static_cast<double>(i + 1) / nBins);
        return e2 - e1;
    }
    return (EmaxMeV - E_low) / nBins;
}

void RunAction::AddGenerated(double E_MeV) {
    const int i = EminMeV < EmaxMeV ? FindBin(E_MeV) : 0;
    if (i < 0) return;

    genCounts[i] += 1.0;

    if (analysisManager and EminMeV < EmaxMeV) {
        analysisManager->FillGenEnergyHist(E_MeV, 1.0);
    }
}


void RunAction::FillDerivedHists() {
    for (int i = 0; i < nBins; ++i) {
        const double nGen = genCounts[i].GetValue();
        double nTrig = trigCounts[i].GetValue();

        const double centerE = BinCenterMeV(i);
        double aEff = 0.0;
        double aEffOpt = 0.0;
        double sens = 0.0;
        double sensOpt = 0.0;
        if (nGen > 0.0) {
            if (nTrig < N_MIN) nTrig = 0;
            aEff = area * (nTrig / nGen);
            sens = aEff;
            sensOpt = aEffOpt;
        }
        effArea[i] = aEff;
        if (fluxDirection.find("isotropic") != std::string::npos) {
            analysisManager->FillSensitivityHist(centerE, sens);
            analysisManager->FillSensitivityOptHist(centerE, sensOpt);
        } else {
            analysisManager->FillEffAreaHist(centerE, aEff);
            analysisManager->FillEffAreaOptHist(centerE, aEffOpt);
        }
    }
}
