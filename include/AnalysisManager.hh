#ifndef ANALYSISMANAGER_HH
#define ANALYSISMANAGER_HH

#include <G4AnalysisManager.hh>
#include <G4AccumulableManager.hh>
#include <G4Accumulable.hh>
#include <G4UnitsTable.hh>
#include <CLHEP/Units/SystemOfUnits.h>
#include <globals.hh>
#include <Sizes.hh>
#include <Configuration.hh>

class AnalysisManager {
public:
    G4String fileName = "GammaDetector";

    explicit AnalysisManager(const std::string&);
    AnalysisManager(const std::string&, int);
    ~AnalysisManager() = default;

    void Open();
    void Close();

    void FillEventRow(G4int eventID, G4int nPrimaries, G4int nInteractions, G4int nEdepHits);

    void FillPrimaryRow(G4int eventID, const G4String& primaryName,
                        G4double E_MeV, const G4ThreeVector& dir,
                        const G4ThreeVector& pos_mm);

    void FillInteractionRow(G4int eventID,
                            G4int trackID, G4int parentID,
                            const G4String& process,
                            const G4String& volumeName,
                            const G4ThreeVector& x_mm,
                            G4int secIndex, const G4String& secName,
                            G4double secE_MeV, const G4ThreeVector& secDir);

    void FillEdepRow(G4int eventID, G4double E0_MeV, G4double edep1_eV, G4double edep2_eV);

    void FillGenEnergyHist(G4double E_MeV, G4double weight = 1.0);
    void FillTrigEnergyHist(G4int channel, G4double E_MeV, G4double weight = 1.0);

    void FillEffAreaHist(G4int channel, G4double E_MeV, G4double value);
    void FillSensitivityHist(G4int channel, G4double E_MeV, G4double value);

private:
    static constexpr G4int nChannels = 2;

    G4int eventNT{-1};
    G4int primaryNT{-1};
    G4int interactionsNT{-1};
    G4int edepNT{-1};

    G4int genEnergyHist{-1};
    G4int trigEnergyHist[nChannels]{-1, -1};
    G4int effAreaHist[nChannels]{-1, -1};
    G4int sensitivityHist[nChannels]{-1, -1};

    G4int nBins{1000};

    void Book();
};


#endif //ANALYSISMANAGER_HH
