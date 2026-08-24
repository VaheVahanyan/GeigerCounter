#ifndef SEPFLUX_HH
#define SEPFLUX_HH

#include "Flux.hh"

#include <fstream>
#include <regex>

class SEPFlux : public Flux {
public:
    explicit SEPFlux();

private:
    std::string path;
    G4int year{};
    G4int order{};

    std::vector<G4double> EList;
    std::vector<G4double> CDF;

    void BuildCDF();

    G4double SampleEnergy() override;
};


#endif //SEPFLUX_HH
