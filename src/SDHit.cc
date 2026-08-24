#include "SDHit.hh"


G4ThreadLocal G4Allocator<SDHit> *SDHitAllocator = nullptr;


void SDHit::AddEdep(G4double val) {
    edep += val;
}


G4double SDHit::GetEdep() const {
    return edep;
}


void SDHit::UpdateTmin(G4double t) {
    if (t < tmin) {
        tmin = t;
    }
}