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
    edepNT = analysisManager->CreateNtuple("edep", "energy deposition per sensitive channel");
    analysisManager->CreateNtupleIColumn("eventID");
    analysisManager->CreateNtupleSColumn("det_name");
    analysisManager->CreateNtupleDColumn("edep_MeV");
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

        genEnergyHist = analysisManager->CreateH1("genEnergyHist",
                                                  "N_{gen} vs E",
                                                  nBins, Emin, Emax, unit, "none", logScheme);

        trigEnergyHist = analysisManager->CreateH1("trigEnergyHist",
                                                   "N_{trig} vs E",
                                                   nBins, Emin, Emax, unit, "none", logScheme);

        if (fluxDirection.find("isotropic") != std::string::npos) {
            sensitivityHist = analysisManager->CreateH1("sensitivityHist",
                                                        "Sensitivity vs E",
                                                        nBins, Emin, Emax, unit, "none", logScheme);
        } else {
            effAreaHist = analysisManager->CreateH1("effAreaHist",
                                                    "A_{eff} vs E",
                                                    nBins, Emin, Emax, unit, "none", logScheme);
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

void AnalysisManager::FillEdepRow(G4int eventID, const G4String& det_name, G4double edep_MeV) {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(edepNT, 0, eventID);
    analysisManager->FillNtupleSColumn(edepNT, 1, det_name);
    analysisManager->FillNtupleDColumn(edepNT, 2, edep_MeV);
    analysisManager->AddNtupleRow(edepNT);
}

void AnalysisManager::FillSiPMEventRow(int eventID, int npeC, int npeV, int npeBV) {
    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(SiPMEventNT, 0, eventID);
    analysisManager->FillNtupleIColumn(SiPMEventNT, 1, npeC);
    analysisManager->FillNtupleIColumn(SiPMEventNT, 2, npeV);
    analysisManager->FillNtupleIColumn(SiPMEventNT, 3, npeBV);
    analysisManager->AddNtupleRow(SiPMEventNT);
}

void AnalysisManager::FillSiPMChannelRow(int eventID, const G4String& subdet, int ch, int npe) {
    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(SiPMChannelNT, 0, eventID);
    analysisManager->FillNtupleSColumn(SiPMChannelNT, 1, subdet);
    analysisManager->FillNtupleIColumn(SiPMChannelNT, 2, ch);
    analysisManager->FillNtupleIColumn(SiPMChannelNT, 3, npe);
    analysisManager->AddNtupleRow(SiPMChannelNT);
}

void AnalysisManager::FillPhotonCountRow(G4int eventID,
                                         G4int npeCrystal, G4int npeVeto,
                                         G4int npeBottomVeto) {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(photonsCountNT, 0, eventID);
    analysisManager->FillNtupleIColumn(photonsCountNT, 1, npeCrystal);
    analysisManager->FillNtupleIColumn(photonsCountNT, 2, npeVeto);
    analysisManager->FillNtupleIColumn(photonsCountNT, 3, npeBottomVeto);
    analysisManager->AddNtupleRow(photonsCountNT);
}

void AnalysisManager::FillPhotonRow(G4int eventID, G4int photonID, const G4String& det_name, G4int det_ch,
                                    G4double energy_eV, G4double x_mm, G4double y_mm, G4double z_mm) {
    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(photonsNT, 0, eventID);
    analysisManager->FillNtupleIColumn(photonsNT, 1, photonID);
    analysisManager->FillNtupleSColumn(photonsNT, 2, det_name);
    analysisManager->FillNtupleIColumn(photonsNT, 3, det_ch);
    analysisManager->FillNtupleDColumn(photonsNT, 4, energy_eV);
    analysisManager->FillNtupleDColumn(photonsNT, 5, x_mm);
    analysisManager->FillNtupleDColumn(photonsNT, 6, y_mm);
    analysisManager->FillNtupleDColumn(photonsNT, 7, z_mm);
    analysisManager->AddNtupleRow(photonsNT);
}


void AnalysisManager::FillGenEnergyHist(G4double E_MeV, G4double weight) {
    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillH1(genEnergyHist, E_MeV, weight);
}

void AnalysisManager::FillTrigEnergyHist(G4double E_MeV, G4double weight) {
    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillH1(trigEnergyHist, E_MeV, weight);
}

void AnalysisManager::FillTrigOptEnergyHist(G4double E_MeV, G4double weight) {
    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillH1(trigOptEnergyHist, E_MeV, weight);
}

void AnalysisManager::FillEffAreaHist(G4double E_MeV, G4double value) {
    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillH1(effAreaHist, E_MeV, value);
}

void AnalysisManager::FillEffAreaOptHist(G4double E_MeV, G4double value) {
    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillH1(effAreaOptHist, E_MeV, value);
}

void AnalysisManager::FillSensitivityHist(G4double E_MeV, G4double value) {
    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillH1(sensitivityHist, E_MeV, value);
}

void AnalysisManager::FillSensitivityOptHist(G4double E_MeV, G4double value) {
    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillH1(sensitivityOptHist, E_MeV, value);
}
