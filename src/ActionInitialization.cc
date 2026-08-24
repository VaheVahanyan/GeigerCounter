#include "ActionInitialization.hh"

using namespace Configuration;

ActionInitialization::ActionInitialization(const G4double a) : EminMeV(Emin),
                                                               EmaxMeV(Emax),
                                                               area(a) {
    if (eThreshold > EmaxMeV) {
        G4Exception("ActionInitialization", "EnergyRange", FatalException,
                    "The energy threshold for a crystal must be less than the maximum value in a given energy range");
    }
}

void ActionInitialization::BuildForMaster() const {
    RunAction* runAct = new RunAction(area);
    SetUserAction(runAct);
}

void ActionInitialization::Build() const {
    RunAction* runAct = new RunAction(area);
    SetUserAction(runAct);

    EventAction* eventAct = new EventAction(runAct->analysisManager, runAct);
    SetUserAction(eventAct);

    PrimaryGeneratorAction* primaryGenerator = new PrimaryGeneratorAction(fluxDirection, fluxType, eThreshold);
    SetUserAction(primaryGenerator);

    if (saveSecondaries) {
        SteppingAction* stepAct = new SteppingAction();
        SetUserAction(stepAct);
    }
}
