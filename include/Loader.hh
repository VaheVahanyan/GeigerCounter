#ifndef LOADER_HH
#define LOADER_HH

#include <G4GlobalConfig.hh>
#include <G4StepLimiter.hh>
#include <G4UserLimits.hh>
#include <G4UImanager.hh>
#include <G4VisManager.hh>
#include <G4PhysicsListHelper.hh>
#include <FTFP_BERT.hh>
#include <QGSP_BIC.hh>
#include <QBBC.hh>
#include <G4IonPhysics.hh>
#include <G4OpticalPhysics.hh>
#include <G4OpticalParameters.hh>
#include <G4StepLimiterPhysics.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4Types.hh>
#include <G4RadioactiveDecayPhysics.hh>
#include <globals.hh>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <G4VisExecutive.hh>
#include <G4UIExecutive.hh>

#include "Geometry.hh"
#include "Sizes.hh"
#include "Configuration.hh"
#include "ActionInitialization.hh"
#include "CountRates.hh"
#include "PostProcessing.hh"

#ifdef G4MULTITHREADED
#include <G4MTRunManager.hh>
#else
#include <G4RunManager.hh>
#endif


class Loader {
    G4String macroFile;
    int numThreads;
    bool useUI;

    G4double area;
    std::vector<G4double> effArea;
    std::vector<G4double> effAreaOpt;

#ifdef G4MULTITHREADED
    G4MTRunManager *runManager;
#else
    G4RunManager *runManager;
#endif

    G4VisManager *visManager;

public:
    Loader(int argc, char **argv);
    ~Loader();

private:
    std::string configPath;
    G4int tube1Only{};
    G4int tube2Only{};
    G4int both{};
    G4int crystalAndVetoOpt{};

    FluxDir dir{};

    [[nodiscard]] std::string ReadValue(const std::string &, const std::string &) const;
    void SaveConfig() const;
    void RunPostProcessing() const;
};

#endif //LOADER_HH
