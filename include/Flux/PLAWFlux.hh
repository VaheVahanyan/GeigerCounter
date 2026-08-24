#ifndef PLAWFLUX_HH
#define PLAWFLUX_HH

#include "Flux.hh"

class PLAWFlux : public Flux {
public:
    explicit PLAWFlux();

private:
    G4double alpha{};
    G4double Epiv{};

    G4double SampleEnergy() override;
};


#endif //PLAWFLUX_HH
