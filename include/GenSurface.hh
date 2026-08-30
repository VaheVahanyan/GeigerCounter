#ifndef GENSURFACE_HH
#define GENSURFACE_HH

#include <G4String.hh>
#include <G4ThreeVector.hh>
#include <G4SystemOfUnits.hh>
#include <G4PhysicalConstants.hh>

#include "Sizes.hh"


class GenSurface {
public:
    enum class Shape { Sphere, Rectangle };

    static GenSurface For(const G4String& fluxDirection);

    static G4ThreeVector PayloadHalfSize();
    static G4ThreeVector PayloadCentre();
    static G4ThreeVector InstrumentUp() { return {0., 1., 0.}; }
    static G4ThreeVector Arrival(G4double theta, G4double phi);
    static G4String DirectionTag();
    static G4double Margin() { return 5. * mm; }

    [[nodiscard]] G4bool IsIsotropic() const { return shape == Shape::Sphere; }
    [[nodiscard]] G4bool IsHemisphere() const { return hemisphere; }

    [[nodiscard]] G4double Radius() const { return radius; }
    [[nodiscard]] G4double HalfU() const { return halfU; }
    [[nodiscard]] G4double HalfV() const { return halfV; }
    [[nodiscard]] G4double Standoff() const { return standoff; }
    [[nodiscard]] G4double Theta() const { return theta; }
    [[nodiscard]] G4double Phi() const { return phi; }
    [[nodiscard]] const G4ThreeVector& Origin() const { return origin; }

    [[nodiscard]] const G4ThreeVector& Axis() const { return axis; }
    [[nodiscard]] const G4ThreeVector& U() const { return u; }
    [[nodiscard]] const G4ThreeVector& V() const { return v; }

    [[nodiscard]] G4double SPerp_cm2() const;
    [[nodiscard]] G4double ProjectedArea_cm2() const;
    [[nodiscard]] G4double GeomFactor_cm2sr() const;
    [[nodiscard]] G4double Norm_cm2() const;

    [[nodiscard]] G4String ShapeName() const;

private:
    Shape shape{Shape::Sphere};
    G4bool hemisphere{false};

    G4double radius{0.};
    G4double halfU{0.};
    G4double halfV{0.};
    G4double standoff{0.};
    G4double theta{0.};
    G4double phi{0.};

    G4ThreeVector origin{0., 0., 0.};
    G4ThreeVector axis{0., 0., -1.};
    G4ThreeVector u{1., 0., 0.};
    G4ThreeVector v{0., 1., 0.};
};

#endif //GENSURFACE_HH
