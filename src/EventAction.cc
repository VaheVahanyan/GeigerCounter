#include "EventAction.hh"

using namespace Sizes;
using namespace Configuration;

EventAction::EventAction(AnalysisManager* an, RunAction* r) : analysisManager(an), run(r) {
    detMap = {
        {"CounterSD/EdepHits", 0, "GeigerCounter"},
    };
    HCIDs.assign(detMap.size(), -1);
}

void EventAction::BeginOfEventAction(const G4Event*) {
    nPrimaries = 0;
    nInteractions = 0;
    nEdepHits = 0;
    hasTube1 = false;
    hasTube2 = false;
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

    nEdepHits = WriteEdepFromSD_(evt, eventID);

    if (saveSecondaries) {
        analysisManager->FillEventRow(eventID, nPrimaries, nInteractions, nEdepHits);
    }

    if (run and hasTube1 && !hasTube2) run->AddTube1Only(1);
    if (run and !hasTube1 && hasTube2) run->AddTube2Only(1);
    if (run and hasTube1 && hasTube2) run->AddBoth(1);
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

int EventAction::WriteEdepFromSD_(const G4Event* evt, int eventID) {
    auto* hce = evt->GetHCofThisEvent();
    if (!hce) return 0;

    auto* sdm = G4SDManager::GetSDMpointer();

    for (size_t i = 0; i < detMap.size(); ++i) {
        if (HCIDs[i] < 0) {
            const auto& hcName = std::get<0>(detMap[i]);
            HCIDs[i] = sdm->GetCollectionID(hcName);
        }
    }

    int nHitsTotal = 0;

    for (size_t i = 0; i < detMap.size(); ++i) {
        const int hcID = HCIDs[i];
        if (hcID < 0) continue;

        auto* hc = dynamic_cast<SDHitCollection*>(hce->GetHC(hcID));
        if (!hc) continue;

        const auto& det_name = std::get<2>(detMap[i]);

        const auto N = hc->GetSize();
        for (unsigned j = 0; j < N; ++j) {
            auto* h = (*hc)[j];
            double edep_MeV = h->edep / MeV;

            if (edep_MeV <= eThreshold) edep_MeV = 0;
            if (edep_MeV > 0.0) {
                if (det_name == "GeigerCounter") MarkTube1();
                else if (det_name == "GeigerCounter") MarkTube2();
                analysisManager->FillEdepRow(eventID, det_name, edep_MeV);
            }
        }
        nHitsTotal += static_cast<int>(N);
    }

    return nHitsTotal;
}
