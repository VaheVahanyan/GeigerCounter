#include "GenSurfaceVis.hh"

#include <algorithm>
#include <cmath>

#include <G4Polyline.hh>
#include <G4Colour.hh>
#include <G4VisAttributes.hh>
#include <G4Box.hh>
#include <G4Sphere.hh>
#include <G4RotationMatrix.hh>
#include <G4Transform3D.hh>

#include "Configuration.hh"

using namespace Configuration;


static void Segment(G4VVisManager* vis, const G4ThreeVector& a, const G4ThreeVector& b, const G4Colour& colour) {
    G4Polyline line;
    line.push_back(a);
    line.push_back(b);
    line.SetVisAttributes(G4VisAttributes(colour));
    vis->Draw(line);
}


static void Loop(G4VVisManager* vis, const std::vector<G4ThreeVector>& pts, const G4Colour& colour) {
    G4Polyline line;
    for (const auto& p : pts) line.push_back(p);
    line.push_back(pts.front());
    line.SetVisAttributes(G4VisAttributes(colour));
    vis->Draw(line);
}


static void Arc(G4VVisManager* vis, const G4ThreeVector& centre, G4double radius,
                const G4ThreeVector& e1, const G4ThreeVector& e2,
                G4double from, G4double to, const G4Colour& colour) {
    constexpr G4int nSteps = 72;
    G4Polyline line;
    for (G4int i = 0; i <= nSteps; ++i) {
        const G4double a = from + (to - from) * i / nSteps;
        line.push_back(centre + radius * (std::cos(a) * e1 + std::sin(a) * e2));
    }
    line.SetVisAttributes(G4VisAttributes(colour));
    vis->Draw(line);
}


void GenSurfaceVis::DrawPayloadBox(G4VVisManager* vis) {
    const G4ThreeVector h = GenSurface::PayloadHalfSize();
    const G4ThreeVector c = GenSurface::PayloadCentre();
    const G4Colour colour = G4Colour::Grey();

    G4ThreeVector corner[8];
    for (G4int i = 0; i < 8; ++i) {
        corner[i] = c + G4ThreeVector((i & 1 ? 1. : -1.) * h.x(),
                                      (i & 2 ? 1. : -1.) * h.y(),
                                      (i & 4 ? 1. : -1.) * h.z());
    }
    for (G4int i = 0; i < 8; ++i) {
        for (const G4int bit : {1, 2, 4}) {
            if (!(i & bit)) Segment(vis, corner[i], corner[i | bit], colour);
        }
    }
}


static void Arrow(G4VVisManager* vis, const G4ThreeVector& tail, const G4ThreeVector& dir,
                  G4double length, const G4ThreeVector& e1, const G4ThreeVector& e2, const G4Colour& colour) {
    const G4ThreeVector tip = tail + dir * length;
    const G4double barb = 0.2 * length;

    Segment(vis, tail, tip, colour);
    Segment(vis, tip, tip - dir * barb + e1 * barb * 0.35, colour);
    Segment(vis, tip, tip - dir * barb - e1 * barb * 0.35, colour);
    Segment(vis, tip, tip - dir * barb + e2 * barb * 0.35, colour);
    Segment(vis, tip, tip - dir * barb - e2 * barb * 0.35, colour);
}


void GenSurfaceVis::DrawRectangle(G4VVisManager* vis, const GenSurface& s) {
    const G4ThreeVector centre = s.Origin() - s.Axis() * s.Standoff();
    const G4ThreeVector du = s.U() * s.HalfU();
    const G4ThreeVector dv = s.V() * s.HalfV();

    const G4RotationMatrix rotation(s.U(), s.V(), s.Axis());
    const G4Box sheet("GenSurfaceSheet", s.HalfU(), s.HalfV(), 0.002 * s.Standoff());
    G4VisAttributes sheetAttributes(G4Colour(1., 1., 0., 0.15));
    sheetAttributes.SetForceSolid(true);
    vis->Draw(sheet, sheetAttributes, G4Transform3D(rotation, centre));

    Loop(vis, {centre - du - dv, centre + du - dv, centre + du + dv, centre - du + dv}, G4Colour::Yellow());

    constexpr G4int nShort = 6;
    const G4double spacing = 2.0 * std::min(s.HalfU(), s.HalfV()) / nShort;
    const G4int nU = std::clamp(static_cast<G4int>(std::lround(2.0 * s.HalfU() / spacing)), 2, 16);
    const G4int nV = std::clamp(static_cast<G4int>(std::lround(2.0 * s.HalfV() / spacing)), 2, 16);
    const G4double length = 0.45 * s.Standoff();

    for (G4int i = 0; i < nU; ++i) {
        for (G4int j = 0; j < nV; ++j) {
            const G4double a = (2.0 * (i + 0.5) / nU - 1.0) * s.HalfU();
            const G4double b = (2.0 * (j + 0.5) / nV - 1.0) * s.HalfV();
            Arrow(vis, centre + s.U() * a + s.V() * b, s.Axis(), length, s.U(), s.V(), G4Colour::Cyan());
        }
    }
}


void GenSurfaceVis::DrawSphere(G4VVisManager* vis, const GenSurface& s) {
    const G4ThreeVector centre = s.Origin();
    const G4double r = s.Radius();
    const G4ThreeVector up = GenSurface::InstrumentUp();
    const G4ThreeVector e1(1., 0., 0.);
    const G4ThreeVector e2(0., 0., 1.);
    const G4Colour colour = G4Colour::Yellow();

    G4double from = 0.;
    G4double to = pi;
    if (s.IsHemisphere()) {
        if (fluxDirection == "isotropic_up") to = halfpi;
        else from = halfpi;
    }

    G4RotationMatrix rotation;
    rotation.rotateX(-halfpi);

    const G4Sphere shell("GenSurfaceShell", 0., r, 0., twopi, from, to - from);
    G4VisAttributes shellAttributes(G4Colour(1., 1., 0., 0.15));
    shellAttributes.SetForceSolid(true);
    vis->Draw(shell, shellAttributes, G4Transform3D(rotation, centre));

    Arc(vis, centre, r, e1, e2, 0., twopi, colour);

    constexpr G4int nMeridians = 12;
    for (G4int i = 0; i < nMeridians; ++i) {
        const G4double psi = twopi * i / nMeridians;
        const G4ThreeVector radial = std::cos(psi) * e1 + std::sin(psi) * e2;
        Arc(vis, centre, r, up, radial, from, to, colour);
    }
}


void GenSurfaceVis::Draw() {
    G4VVisManager* vis = G4VVisManager::GetConcreteInstance();
    if (!vis) return;

    const GenSurface s = GenSurface::For(fluxDirection);

    DrawPayloadBox(vis);
    if (s.IsIsotropic()) DrawSphere(vis, s);
    else DrawRectangle(vis, s);
}


G4VisExtent GenSurfaceVis::Extent() {
    const GenSurface s = GenSurface::For(fluxDirection);
    const G4ThreeVector c = GenSurface::PayloadCentre();

    G4double r = s.Standoff();
    if (!s.IsIsotropic()) {
        r = std::sqrt(s.Standoff() * s.Standoff() + s.HalfU() * s.HalfU() + s.HalfV() * s.HalfV());
    }
    return {G4Point3D(c), r + GenSurface::Margin()};
}
