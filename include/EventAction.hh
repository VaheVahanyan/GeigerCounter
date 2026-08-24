#ifndef EVENTACTION_HH
#define EVENTACTION_HH

#include <G4UserEventAction.hh>
#include <G4ThreeVector.hh>
#include <G4String.hh>
#include <G4Event.hh>
#include <G4Run.hh>
#include <G4RunManager.hh>
#include <G4SDManager.hh>
#include <G4HCofThisEvent.hh>
#include <G4SystemOfUnits.hh>
#include <cfloat>
#include <vector>

#include "Geometry.hh"
#include "RunAction.hh"
#include "Sizes.hh"
#include "Configuration.hh"
#include "AnalysisManager.hh"
#include "SDHit.hh"
#include "SiPMOpticalSD.hh"

class G4Event;
class Geometry;

struct PrimaryRec {
    int index = 0;
    int pdg = 0;
    G4String name;
    double E_MeV = 0.0;
    G4ThreeVector dir;
    G4ThreeVector pos_mm;  // mm
    double t0_ns = 0.0;    // ns
};

struct InteractionRec {
    int trackID = -1;
    int parentID = -1;
    G4String process;
    G4String volumeName;
    int volumeID = -1;
    G4ThreeVector pos_mm;  // mm
    double t_ns = 0.0;     // ns

    int secIndex = -1;
    int secPDG = 0;
    G4String secName;
    double secE_MeV = 0.0;
    G4ThreeVector secDir;
};

class EventAction : public G4UserEventAction {
public:
    std::vector<PrimaryRec> primBuf;
    std::vector<InteractionRec> interBuf;

    EventAction(AnalysisManager *, RunAction *);
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event *) override;
    void EndOfEventAction(const G4Event *) override;

private:
    void WritePrimaries_(int eventID);
    int WriteInteractions_(int eventID);
    int WriteEdepFromSD_(const G4Event *evt, int eventID);

    void MarkTube1() { hasTube1 = true; }
    void MarkTube2() { hasTube2 = true; }

    AnalysisManager *analysisManager = nullptr;

    std::vector<std::tuple<G4String, int, G4String>> detMap;
    std::vector<int> HCIDs;

    int nPrimaries = 0;
    int nInteractions = 0;
    int nEdepHits = 0;

    RunAction* run = nullptr;
    bool hasTube1 = false;
    bool hasTube2 = false;
};

#endif //EVENTACTION_HH
