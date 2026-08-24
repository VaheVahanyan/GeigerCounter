#ifndef UNIFORMFLUX_HH
#define UNIFORMFLUX_HH

#include "Flux.hh"

class UniformFlux : public Flux {
public:
    explicit UniformFlux();

private:
    std::vector<G4String> particles;
    std::vector<G4double> fractions;
    std::vector<G4double> EminVec;
    std::vector<G4double> EmaxVec;

    static std::vector<G4String> Split(const G4String &line);
    static std::vector<G4double> ParseDoubles(const G4String &line);
    size_t SampleIndex() const;

    G4double SampleEnergy() override;
};


#endif //UNIFORMFLUX_HH