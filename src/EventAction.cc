#include "EventAction.hh"

using namespace Sizes;
using namespace Configuration;

EventAction::EventAction(AnalysisManager* an, RunAction* r) : analysisManager(an), run(r) {
}

void EventAction::BeginOfEventAction(const G4Event*) {
    nPrimaries = 0;
    nInteractions = 0;
    nEdepHits = 0;
    hasTube1 = false;
    hasTube2 = false;
    edepTube[0] = 0.0;
    edepTube[1] = 0.0;
}

void EventAction::EndOfEventAction(const G4Event* evt) {
    const int eventID = evt->GetEventID();

    WritePrimaries_(eventID);
    nPrimaries = static_cast<int>(primBuf.size());

    double primaryE_MeV = -1.0;
    if (!primBuf.empty()) {
        primaryE_MeV = primBuf.front().E_MeV;
        if (run) {
            run->AddGenerated(primaryE_MeV);
        }
    }
    primBuf.clear();

    nInteractions = WriteInteractions_(eventID);
    interBuf.clear();

    nEdepHits = ReadEdepFromSD_(evt);

    if (saveSecondaries) {
        analysisManager->FillEventRow(eventID, nPrimaries, nInteractions, nEdepHits);
    }

    if (hasTube1 || hasTube2) {
        analysisManager->FillEdepRow(eventID, primaryE_MeV, edepTube[0] / eV, edepTube[1] / eV);
    }

    if (run) {
        if (hasTube1 && !hasTube2) run->AddTube1Only(1);
        if (!hasTube1 && hasTube2) run->AddTube2Only(1);
        if (hasTube1 && hasTube2) run->AddBoth(1);

        if (hasTube1) run->AddTriggered(primaryE_MeV, 0);
        if (hasTube1 && hasTube2) run->AddTriggered(primaryE_MeV, 1);
    }
}

void EventAction::WritePrimaries_(int eventID) {
    for (const auto& p : primBuf) {
        analysisManager->FillPrimaryRow(eventID, p.name, p.E_MeV, p.dir, p.pos_mm);
    }
}

int EventAction::WriteInteractions_(int eventID) {
    if (saveSecondaries) {
        for (const auto& r : interBuf) {
            analysisManager->FillInteractionRow(eventID,
                                                r.trackID, r.parentID,
                                                r.process, r.volumeName, r.pos_mm,
                                                r.secIndex, r.secName,
                                                r.secE_MeV, r.secDir);
        }
    }
    return static_cast<int>(interBuf.size());
}

int EventAction::ReadEdepFromSD_(const G4Event* evt) {
    auto* hce = evt->GetHCofThisEvent();
    if (!hce) return 0;

    if (edepHCID < 0) {
        edepHCID = G4SDManager::GetSDMpointer()->GetCollectionID("CounterSD/EdepHits");
    }
    if (edepHCID < 0) return 0;

    auto* hc = dynamic_cast<SDHitCollection*>(hce->GetHC(edepHCID));
    if (!hc) return 0;

    const auto N = hc->GetSize();
    for (unsigned j = 0; j < N; ++j) {
        const auto* h = (*hc)[j];
        if (h->volumeID < 0 || h->volumeID > 1) continue;
        edepTube[h->volumeID] += h->edep;
    }

    for (auto& e : edepTube) {
        if (e <= eThreshold) e = 0.0;
    }

    hasTube1 = edepTube[0] > 0.0;
    hasTube2 = edepTube[1] > 0.0;

    return static_cast<int>(N);
}
