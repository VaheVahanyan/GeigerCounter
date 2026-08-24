#ifndef ACTIONINITIALIZATION_HH
#define ACTIONINITIALIZATION_HH

#include <G4VUserActionInitialization.hh>

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "Configuration.hh"
#include "Geometry.hh"

class ActionInitialization : public G4VUserActionInitialization {
public:
    ActionInitialization(G4double);
    ~ActionInitialization() override = default;

    void BuildForMaster() const override;
    void Build() const override;

private:
    G4double EminMeV;
    G4double EmaxMeV;
    G4double area;
};

#endif //ACTIONINITIALIZATION_HH
