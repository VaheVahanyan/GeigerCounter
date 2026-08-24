#include "GenSurface.hh"

using namespace Sizes;


static void BoxWithLidExtent(G4double& xMin, G4double& xMax,
                             G4double& yMin, G4double& yMax,
                             G4double& zMin, G4double& zMax) {
    xMin = -Box::length;
    xMax = Box::length;

    yMin = -Box::width;
    yMax = Box::width;

    zMin = -Box::thickness - Box::height;
    zMax = Box::height + Box::thickness;
}


G4ThreeVector GenSurface::PayloadHalfSize() {
    G4double xMin, xMax, yMin, yMax, zMin, zMax;
    BoxWithLidExtent(xMin, xMax, yMin, yMax, zMin, zMax);

    return {0.5 * (xMax - xMin), 0.5 * (yMax - yMin), 0.5 * (zMax - zMin)};
}


G4ThreeVector GenSurface::PayloadCentre() {
    G4double xMin, xMax, yMin, yMax, zMin, zMax;
    BoxWithLidExtent(xMin, xMax, yMin, yMax, zMin, zMax);

    return {0.5 * (xMax + xMin), 0.5 * (yMax + yMin), 0.5 * (zMax + zMin)};
}


GenSurface GenSurface::For(const G4String& fluxDirection) {
    const G4ThreeVector half = PayloadHalfSize();

    GenSurface s;
    s.origin = PayloadCentre();
    s.standoff = half.mag() + Margin();

    if (fluxDirection.find("isotropic") != std::string::npos) {
        s.shape = Shape::Sphere;
        s.hemisphere = fluxDirection != "isotropic";
        s.radius = half.mag() + Margin();
        return s;
    }

    s.shape = Shape::Rectangle;

    if (fluxDirection == "vertical_up" || fluxDirection == "vertical_down") {
        s.axis = fluxDirection == "vertical_up" ? G4ThreeVector(0., 0., 1.) : G4ThreeVector(0., 0., -1.);
        s.u = G4ThreeVector(1., 0., 0.);
        s.v = G4ThreeVector(0., 1., 0.);
        s.halfU = half.x();
        s.halfV = half.y();
    } else {
        s.axis = G4ThreeVector(-1., 0., 0.);
        s.u = G4ThreeVector(0., 1., 0.);
        s.v = G4ThreeVector(0., 0., 1.);
        s.halfU = half.y();
        s.halfV = half.z();
    }

    return s;
}


G4double GenSurface::SPerp_cm2() const {
    if (shape == Shape::Sphere) {
        return pi * (radius / cm) * (radius / cm);
    }
    return 4.0 * (halfU / cm) * (halfV / cm);
}


G4double GenSurface::GeomFactor_cm2sr() const {
    if (shape != Shape::Sphere) return 0.0;
    const G4double r_cm = radius / cm;
    return (hemisphere ? 2.0 : 4.0) * pi * pi * r_cm * r_cm;
}


G4double GenSurface::Norm_cm2() const {
    return shape == Shape::Sphere ? GeomFactor_cm2sr() : SPerp_cm2();
}


G4String GenSurface::ShapeName() const {
    if (shape == Shape::Sphere) return hemisphere ? "hemisphere" : "sphere";
    return "rectangle";
}
