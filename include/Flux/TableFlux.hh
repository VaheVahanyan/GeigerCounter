#ifndef TABLEFLUX_HH
#define TABLEFLUX_HH

#include <vector>
#include <string>
#include <fstream>
#include <regex>
#include <algorithm>
#include <limits>
#include <cmath>

#include <G4Types.hh>
#include <G4String.hh>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>

#include "Flux.hh"


class TableFlux : public Flux {
public:
    explicit TableFlux();

private:
    G4String path;

    std::vector<G4double> EList;
    std::vector<G4double> FluxList;;
    std::vector<G4double> CDF;

    void BuildCDF();

    G4double SampleEnergy() override;
};


#endif //TABLEFLUX_HH
