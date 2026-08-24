#ifndef RUNACTION_HH
#define RUNACTION_HH

#include <G4UserRunAction.hh>
#include <globals.hh>
#include <G4SystemOfUnits.hh>
#include <G4Accumulable.hh>
#include <G4AccumulableManager.hh>
#include <G4Run.hh>
#include <G4ios.hh>
#include <G4UnitsTable.hh>
#include <Randomize.hh>
#include <G4AnalysisManager.hh>
#include <G4Threading.hh>
#include <array>
#include <iomanip>
#include <sstream>

#include "Sizes.hh"
#include "Configuration.hh"
#include "AnalysisManager.hh"

struct ParticleCounts {
    G4int tube_1 = 0;
    G4int tube_2 = 0;
    G4int both = 0;
};

class RunAction : public G4UserRunAction {
public:
    AnalysisManager *analysisManager;

    RunAction();
    RunAction(double Agen_cm2);
    ~RunAction() override;

    void BeginOfRunAction(const G4Run *) override;
    void EndOfRunAction(const G4Run *) override;

    void AddTube1Only(const G4int v) { tube1Only += v; }
    void AddTube2Only(const G4int v) { tube2Only += v; }
    void AddBoth(const G4int v) { bothTubs += v; }

    void AddGenerated(double E_MeV);
    void AddTriggered(double E_MeV, int channel);

    [[nodiscard]] const ParticleCounts& GetCounts() const { return totals; }

    [[nodiscard]] G4int GetNGenerated() const { return nGenerated; }
    [[nodiscard]] G4int GetN1() const { return totals.tube_1 + totals.both; }
    [[nodiscard]] G4int GetNTelescope() const { return totals.both; }

    [[nodiscard]] const std::vector<double>& GetEffArea() const { return effArea[0]; }
    [[nodiscard]] const std::vector<double>& GetEffAreaTelescope() const { return effArea[1]; }

private:
    G4Accumulable<G4int> tube1Only{0};   // First && !Second
    G4Accumulable<G4int> tube2Only{0};   // !First && Second
    G4Accumulable<G4int> bothTubs{0};   // First && Second
    ParticleCounts totals{};
    G4int nGenerated{0};

    double EminMeV{0.0};
    double EmaxMeV{0.0};
    double area{0.0};

    double logEmin{0.0};
    double logEmax{0.0};
    double invDlogE{0.0};
    double invDlinearE{0.0};

    std::vector<G4Accumulable<G4double>> genCounts;
    std::array<std::vector<G4Accumulable<G4double>>, 2> trigCounts;
    std::array<std::vector<G4double>, 2> effArea;

    [[nodiscard]] int FindBin(double E_MeV) const;
    [[nodiscard]] double BinCenterMeV(int i) const;
    [[nodiscard]] double BinWidthMeV(int i) const;

    void BookAccumulables();
    void FillDerivedHists();
};

#endif //RUNACTION_HH
