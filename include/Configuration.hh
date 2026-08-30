#ifndef CONFIGURATION_HH
#define CONFIGURATION_HH

#include <G4String.hh>
#include <G4Types.hh>
#include <G4SystemOfUnits.hh>


namespace Configuration
{
    inline G4String geometryType{"full"};

    inline G4String fluxType{"Uniform"};
    inline G4String fluxDirection{"isotropic"};
    inline G4double beamTheta{0.};
    inline G4double beamPhi{0.};

    inline G4double eThreshold{0 * MeV};

    inline G4double gasDensity{0.5597 * mg / cm3};
    inline G4double gasTemperature{293.15 * kelvin};
    inline G4double counterRangeCut{1 * um};
    inline G4double cutsEnergyMin{250 * eV};
    inline G4double cutsEnergyMax{100 * GeV};

    inline G4int nBins{1000};
    inline G4String outputFile{"GeigerCounter.root"};
    inline G4bool saveSecondaries{false};

    inline G4double Emin{-1};
    inline G4double Emax{-1};

    inline G4double N_MIN{10};

    inline G4bool isLogBin{false};
}


#endif //CONFIGURATION_HH
