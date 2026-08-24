#include "SteppingAction.hh"


void SteppingAction::UserSteppingAction(const G4Step* step) {
    const auto* post = step->GetPostStepPoint();
    const auto* postProc = post->GetProcessDefinedStep();

    const G4ThreeVector x = post->GetPosition() / mm;
    const G4double t = post->GetGlobalTime() / ns;

    const auto* touch = post->GetTouchable();
    G4String volName = "World";
    int copyNo = -1;
    if (touch && touch->GetVolume()) {
        volName = touch->GetVolume()->GetName();
        copyNo = touch->GetVolume()->GetCopyNo();
    }

    const auto* track = step->GetTrack();
    const G4int trackID = track->GetTrackID();
    const G4int parentID = track->GetParentID();

    auto secs = step->GetSecondaryInCurrentStep();
    auto* ea = static_cast<EventAction*>(G4EventManager::GetEventManager()->GetUserEventAction());
    if (ea) {
        if (secs && !secs->empty()) {
            for (size_t i = 0; i < secs->size(); ++i) {
                const auto* sc = (*secs)[i];
                const auto* cp = sc->GetCreatorProcess();

                InteractionRec rec;
                rec.trackID = trackID;
                rec.parentID = parentID;
                rec.process = cp ? cp->GetProcessName() : "unknown";
                rec.volumeName = volName;
                rec.volumeID = copyNo;
                rec.pos_mm = x;
                rec.t_ns = t;
                rec.secIndex = static_cast<G4int>(i);
                rec.secPDG = sc->GetDefinition()->GetPDGEncoding();
                rec.secName = sc->GetDefinition()->GetParticleName();
                rec.secE_MeV = sc->GetKineticEnergy() / MeV;
                rec.secDir = sc->GetMomentumDirection();
                ea->interBuf.emplace_back(std::move(rec));
            }
            return;
        }

        if (postProc) {
            const auto& pname = postProc->GetProcessName();
            if (pname != "Transportation") {
                InteractionRec rec;
                rec.trackID = trackID;
                rec.parentID = parentID;
                rec.process = pname;
                rec.volumeName = volName;
                rec.volumeID = copyNo;
                rec.pos_mm = x;
                rec.t_ns = t;
                rec.secIndex = -1;
                rec.secPDG = 0;
                rec.secName = "";
                rec.secE_MeV = 0.0;
                rec.secDir = G4ThreeVector();
                ea->interBuf.emplace_back(std::move(rec));
            }
        }
    }
}
