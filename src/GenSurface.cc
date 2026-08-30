#include "GenSurface.hh"

#include <cmath>
#include <sstream>

#include "Configuration.hh"

using namespace Sizes;
using namespace Configuration;


static void FullExtent(G4double& xMin, G4double& xMax,
                       G4double& yMin, G4double& yMax,
                       G4double& zMin, G4double& zMax) {
    xMin = -Box::length;
    xMax = Box::length;

    yMin = -Box::width;
    yMax = Box::width;

    zMin = -Box::thickness - Box::height;
    zMax = Box::height + Box::thickness;
}


static void TelescopeExtent(G4double& xMin, G4double& xMax,
                            G4double& yMin, G4double& yMax,
                            G4double& zMin, G4double& zMax) {
    const G4double tubeUpperY = Box::width - (Box::thickness + Box::cavity1Gap) * 2 - GeigerCounter::radius;
    const G4double tubeLowerY = Box::width - (Box::thickness + Box::cavity1Gap + Filter::width + Filter::gap * 2) * 2
                                - GeigerCounter::radius * 3;
    const G4double filterY = Box::width - (Box::thickness + GeigerCounter::radius + Box::cavity1Gap + Filter::gap) * 2
                             - Filter::width;
    const G4double retainerY = Box::width - Box::cavity1Width - Box::thickness * 2;
    const G4double lowestTubeY = geometryType == "single" ? tubeUpperY : tubeLowerY;

    xMax = GeigerCounter::length + GeigerCounter::contactLength * 2;
    xMin = -xMax;

    yMax = std::max({tubeUpperY + GeigerCounter::radius, filterY + Filter::width, retainerY + Retainer::width});
    yMin = std::min({lowestTubeY - GeigerCounter::radius, filterY - Filter::width, retainerY - Retainer::width});

    zMax = std::max({GeigerCounter::radius, Filter::height, Retainer::height});
    zMin = -zMax;
}


static void PayloadExtent(G4double& xMin, G4double& xMax,
                          G4double& yMin, G4double& yMax,
                          G4double& zMin, G4double& zMax) {
    if (geometryType == "full") {
        FullExtent(xMin, xMax, yMin, yMax, zMin, zMax);
    } else {
        TelescopeExtent(xMin, xMax, yMin, yMax, zMin, zMax);
    }
}


G4ThreeVector GenSurface::PayloadHalfSize() {
    G4double xMin, xMax, yMin, yMax, zMin, zMax;
    PayloadExtent(xMin, xMax, yMin, yMax, zMin, zMax);

    return {0.5 * (xMax - xMin), 0.5 * (yMax - yMin), 0.5 * (zMax - zMin)};
}


G4ThreeVector GenSurface::PayloadCentre() {
    G4double xMin, xMax, yMin, yMax, zMin, zMax;
    PayloadExtent(xMin, xMax, yMin, yMax, zMin, zMax);

    return {0.5 * (xMax + xMin), 0.5 * (yMax + yMin), 0.5 * (zMax + zMin)};
}


G4ThreeVector GenSurface::Arrival(const G4double theta, const G4double phi) {
    return {std::sin(theta) * std::cos(phi), std::cos(theta), std::sin(theta) * std::sin(phi)};
}


G4String GenSurface::DirectionTag() {
    if (fluxDirection.find("isotropic") != std::string::npos) return fluxDirection;

    std::ostringstream ss;
    ss << fluxDirection << "_t" << beamTheta / deg << "_p" << beamPhi / deg;
    return ss.str();
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
    s.theta = beamTheta;
    s.phi = beamPhi;
    s.axis = (-Arrival(s.theta, s.phi)).unit();

    const G4double ax = std::fabs(s.axis.x());
    const G4double ay = std::fabs(s.axis.y());
    const G4double az = std::fabs(s.axis.z());

    G4ThreeVector e(1., 0., 0.);
    if (ay <= ax && ay <= az) e.set(0., 1., 0.);
    else if (az <= ax && az <= ay) e.set(0., 0., 1.);

    s.u = (e - s.axis * e.dot(s.axis)).unit();
    s.v = s.axis.cross(s.u).unit();

    s.halfU = half.x() * std::fabs(s.u.x()) + half.y() * std::fabs(s.u.y()) + half.z() * std::fabs(s.u.z());
    s.halfV = half.x() * std::fabs(s.v.x()) + half.y() * std::fabs(s.v.y()) + half.z() * std::fabs(s.v.z());

    return s;
}


G4double GenSurface::SPerp_cm2() const {
    if (shape == Shape::Sphere) {
        return pi * (radius / cm) * (radius / cm);
    }
    return 4.0 * (halfU / cm) * (halfV / cm);
}


G4double GenSurface::ProjectedArea_cm2() const {
    if (shape == Shape::Sphere) return SPerp_cm2();

    const G4ThreeVector half = PayloadHalfSize() / cm;
    return 4.0 * (half.y() * half.z() * std::fabs(axis.x())
                  + half.x() * half.z() * std::fabs(axis.y())
                  + half.x() * half.y() * std::fabs(axis.z()));
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
