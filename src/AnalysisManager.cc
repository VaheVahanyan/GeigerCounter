#include "AnalysisManager.hh"

using namespace Sizes;
using namespace Configuration;

AnalysisManager::AnalysisManager(const std::string& fName) : fileName(fName) {
    Book();
}

AnalysisManager::AnalysisManager(const std::string& fName, const int bins) :
    fileName(fName), nBins(bins) {
    Book();
}


void AnalysisManager::Book() {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType("root");
    analysisManager->SetFileName(fileName);
    analysisManager->SetVerboseLevel(0);
    analysisManager->SetNtupleActivation(true);

#ifdef G4MULTITHREADED
    analysisManager->SetNtupleMerging(true);
#endif
    edepNT = analysisManager->CreateNtuple("edep", "energy deposit per counter, one row per event with a deposit");
    analysisManager->CreateNtupleIColumn("eventID");
    analysisManager->CreateNtupleDColumn("E0_MeV");
    analysisManager->CreateNtupleDColumn("edep1_eV");
    analysisManager->CreateNtupleDColumn("edep2_eV");
    analysisManager->FinishNtuple(edepNT);

    primaryNT = analysisManager->CreateNtuple("primary", "per-primary particles");
    analysisManager->CreateNtupleIColumn("eventID");
    analysisManager->CreateNtupleSColumn("primary_name");
    analysisManager->CreateNtupleDColumn("E_MeV");
    analysisManager->CreateNtupleDColumn("dir_x");
    analysisManager->CreateNtupleDColumn("dir_y");
    analysisManager->CreateNtupleDColumn("dir_z");
    analysisManager->CreateNtupleDColumn("pos_x_mm");
    analysisManager->CreateNtupleDColumn("pos_y_mm");
    analysisManager->CreateNtupleDColumn("pos_z_mm");
    analysisManager->FinishNtuple(primaryNT);

    if (saveSecondaries) {
        interactionsNT = analysisManager->CreateNtuple("interactions",
                                                       "inelastic/compton/photo/conv vertices and secondaries");
        analysisManager->CreateNtupleIColumn("eventID");
        analysisManager->CreateNtupleIColumn("trackID");
        analysisManager->CreateNtupleIColumn("parentID");
        analysisManager->CreateNtupleSColumn("process");
        analysisManager->CreateNtupleSColumn("volume_name");
        analysisManager->CreateNtupleDColumn("x_mm");
        analysisManager->CreateNtupleDColumn("y_mm");
        analysisManager->CreateNtupleDColumn("z_mm");
        analysisManager->CreateNtupleDColumn("t_ns");
        analysisManager->CreateNtupleIColumn("sec_index");
        analysisManager->CreateNtupleSColumn("sec_name");
        analysisManager->CreateNtupleDColumn("sec_E_MeV");
        analysisManager->CreateNtupleDColumn("sec_dir_x");
        analysisManager->CreateNtupleDColumn("sec_dir_y");
        analysisManager->CreateNtupleDColumn("sec_dir_z");
        analysisManager->FinishNtuple(interactionsNT);

        eventNT = analysisManager->CreateNtuple("event", "per-event summary");
        analysisManager->CreateNtupleIColumn("eventID");
        analysisManager->CreateNtupleIColumn("n_primaries");
        analysisManager->CreateNtupleIColumn("n_interactions");
        analysisManager->CreateNtupleIColumn("n_edep_hits");
        analysisManager->FinishNtuple(eventNT);
    }


    if (Emin < Emax) {
        const G4String unit = "MeV";
        G4String logScheme = isLogBin ? "log" : "linear";

        const G4String suffix[nChannels] = {"", "Tel"};
        const G4String label[nChannels] = {"counter 1", "telescope"};

        genEnergyHist = analysisManager->CreateH1("genEnergyHist",
                                                  "N_{gen} vs E",
                                                  nBins, Emin, Emax, unit, "none", logScheme);

        for (G4int c = 0; c < nChannels; ++c) {
            trigEnergyHist[c] = analysisManager->CreateH1("trigEnergyHist" + suffix[c],
                                                          "N_{trig} vs E, " + label[c],
                                                          nBins, Emin, Emax, unit, "none", logScheme);

            if (fluxDirection.find("isotropic") != std::string::npos) {
                sensitivityHist[c] = analysisManager->CreateH1("sensitivityHist" + suffix[c],
                                                               "Sensitivity vs E, " + label[c],
                                                               nBins, Emin, Emax, unit, "none", logScheme);
            } else {
                effAreaHist[c] = analysisManager->CreateH1("effAreaHist" + suffix[c],
                                                           "A_{eff} vs E, " + label[c],
                                                           nBins, Emin, Emax, unit, "none", logScheme);
            }
        }
    }
}

void AnalysisManager::Open() {
    G4AnalysisManager::Instance()->OpenFile(fileName);
}

void AnalysisManager::Close() {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
}

void AnalysisManager::FillEventRow(G4int eventID, G4int nPrimaries, G4int nInteractions, G4int nEdepHits) {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(eventNT, 0, eventID);
    analysisManager->FillNtupleIColumn(eventNT, 1, nPrimaries);
    analysisManager->FillNtupleIColumn(eventNT, 2, nInteractions);
    analysisManager->FillNtupleIColumn(eventNT, 3, nEdepHits);
    analysisManager->AddNtupleRow(eventNT);
}

void AnalysisManager::FillPrimaryRow(G4int eventID, const G4String& primaryName,
                                     G4double E_MeV, const G4ThreeVector& dir,
                                     const G4ThreeVector& pos_mm) {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(primaryNT, 0, eventID);
    analysisManager->FillNtupleSColumn(primaryNT, 1, primaryName);
    analysisManager->FillNtupleDColumn(primaryNT, 2, E_MeV);
    analysisManager->FillNtupleDColumn(primaryNT, 3, dir.x());
    analysisManager->FillNtupleDColumn(primaryNT, 4, dir.y());
    analysisManager->FillNtupleDColumn(primaryNT, 5, dir.z());
    analysisManager->FillNtupleDColumn(primaryNT, 6, pos_mm.x());
    analysisManager->FillNtupleDColumn(primaryNT, 7, pos_mm.y());
    analysisManager->FillNtupleDColumn(primaryNT, 8, pos_mm.z());
    analysisManager->AddNtupleRow(primaryNT);
}

void AnalysisManager::FillInteractionRow(G4int eventID,
                                         G4int trackID, G4int parentID,
                                         const G4String& process,
                                         const G4String& volumeName,
                                         const G4ThreeVector& x_mm,
                                         G4int secIndex, const G4String& secName,
                                         G4double secE_MeV, const G4ThreeVector& secDir) {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(interactionsNT, 0, eventID);
    analysisManager->FillNtupleIColumn(interactionsNT, 1, trackID);
    analysisManager->FillNtupleIColumn(interactionsNT, 2, parentID);
    analysisManager->FillNtupleSColumn(interactionsNT, 3, process);
    analysisManager->FillNtupleSColumn(interactionsNT, 4, volumeName);
    analysisManager->FillNtupleDColumn(interactionsNT, 5, x_mm.x());
    analysisManager->FillNtupleDColumn(interactionsNT, 6, x_mm.y());
    analysisManager->FillNtupleDColumn(interactionsNT, 7, x_mm.z());
    analysisManager->FillNtupleIColumn(interactionsNT, 8, secIndex);
    analysisManager->FillNtupleSColumn(interactionsNT, 9, secName);
    analysisManager->FillNtupleDColumn(interactionsNT, 10, secE_MeV);
    analysisManager->FillNtupleDColumn(interactionsNT, 11, secDir.x());
    analysisManager->FillNtupleDColumn(interactionsNT, 12, secDir.y());
    analysisManager->FillNtupleDColumn(interactionsNT, 13, secDir.z());
    analysisManager->AddNtupleRow(interactionsNT);
}

void AnalysisManager::FillEdepRow(G4int eventID, G4double E0_MeV, G4double edep1_eV, G4double edep2_eV) {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(edepNT, 0, eventID);
    analysisManager->FillNtupleDColumn(edepNT, 1, E0_MeV);
    analysisManager->FillNtupleDColumn(edepNT, 2, edep1_eV);
    analysisManager->FillNtupleDColumn(edepNT, 3, edep2_eV);
    analysisManager->AddNtupleRow(edepNT);
}

void AnalysisManager::FillGenEnergyHist(G4double E_MeV, G4double weight) {
    if (genEnergyHist < 0) return;
    G4AnalysisManager::Instance()->FillH1(genEnergyHist, E_MeV, weight);
}

void AnalysisManager::FillTrigEnergyHist(G4int channel, G4double E_MeV, G4double weight) {
    if (channel < 0 || channel >= nChannels || trigEnergyHist[channel] < 0) return;
    G4AnalysisManager::Instance()->FillH1(trigEnergyHist[channel], E_MeV, weight);
}

void AnalysisManager::FillEffAreaHist(G4int channel, G4double E_MeV, G4double value) {
    if (channel < 0 || channel >= nChannels || effAreaHist[channel] < 0) return;
    G4AnalysisManager::Instance()->FillH1(effAreaHist[channel], E_MeV, value);
}

void AnalysisManager::FillSensitivityHist(G4int channel, G4double E_MeV, G4double value) {
    if (channel < 0 || channel >= nChannels || sensitivityHist[channel] < 0) return;
    G4AnalysisManager::Instance()->FillH1(sensitivityHist[channel], E_MeV, value);
}
