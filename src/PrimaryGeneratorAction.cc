#include "PrimaryGeneratorAction.hh"


PrimaryGeneratorAction::PrimaryGeneratorAction(G4String fDir, const G4String& fluxType, const G4double cThreshold)
    : particleGun(new G4ParticleGun(1)),
      fluxDirection(std::move(fDir)),
      eCrystalThreshold(cThreshold) {
    genSurface = GenSurface::For(fluxDirection);
    center = genSurface.Origin();

    std::vector<G4String> fluxDirList = {
        "isotropic", "isotropic_up", "isotropic_down", "vertical_up", "vertical_down", "horizontal"
    };
    if (std::find(fluxDirList.begin(), fluxDirList.end(), fluxDirection) == fluxDirList.end()) {
        G4Exception("PrimaryGeneratorAction::GeneratePrimaries", "FluxDirection", FatalException,
                    ("Flux direction is not implemented: " + fluxDirection +
                        ".\nAvailable flux directions: isotropic, isotropic_up, isotropic_down, vertical_up," +
                        " vertical_down, horizontal").c_str());
    }
    std::vector<G4String> fluxTypeList = {"Uniform", "PLAW", "COMP", "SEP", "Galactic", "Table"};
    if (std::find(fluxTypeList.begin(), fluxTypeList.end(), fluxType) == fluxTypeList.end()) {
        G4Exception("PrimaryGeneratorAction::GeneratePrimaries", "FluxType", FatalException,
                    ("Flux type not found: " + fluxType + ".\nAvailable flux types: Uniform, PLAW, SEP, Galactic")
                    .
                    c_str());
    }

    if (fluxType == "Uniform") {
        flux = new UniformFlux();
    } else if (fluxType == "PLAW") {
        flux = new PLAWFlux();
    } else if (fluxType == "COMP") {
        flux = new COMPFlux();
    } else if (fluxType == "SEP") {
        flux = new SEPFlux();
    } else if (fluxType == "Galactic") {
        flux = new GalacticFlux();
    } else if (fluxType == "Table") {
        flux = new TableFlux();
    }
}


PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete particleGun;
}


void PrimaryGeneratorAction::GenerateOnSphere(G4ThreeVector& pos, G4ThreeVector& dir) const {
    G4double u = 0;
    if (fluxDirection == "isotropic") {
        u = 2.0 * G4UniformRand() - 1.0; // cos(theta) ~ U[-1,1]
    } else if (fluxDirection == "isotropic_up") {
        u = G4UniformRand(); // cos(theta) ~ U[0,1]
    } else if (fluxDirection == "isotropic_down") {
        u = -G4UniformRand(); // cos(theta) ~ U[-1,0]
    }
    const G4double phi = 2.0 * M_PI * G4UniformRand();
    const G4double l = std::sqrt(std::max(0.0, 1.0 - u * u));
    const G4ThreeVector rhat(l * std::cos(phi), l * std::sin(phi), u);

    pos = center + genSurface.Radius() * rhat;

    const G4ThreeVector z = rhat.unit();
    const G4ThreeVector a = std::fabs(z.z()) < 0.999 ? G4ThreeVector(0, 0, 1) : G4ThreeVector(1, 0, 0);
    const G4ThreeVector x = z.cross(a).unit();
    const G4ThreeVector y = z.cross(x).unit();

    const G4double ksi = G4UniformRand();
    const G4double sinTh = std::sqrt(ksi);
    const G4double cosTh = std::sqrt(1.0 - ksi);
    const G4double phi2 = 2.0 * M_PI * G4UniformRand();

    const G4ThreeVector v_local(sinTh * std::cos(phi2), sinTh * std::sin(phi2), cosTh);

    dir = -(v_local.x() * x + v_local.y() * y + v_local.z() * z);
    dir = dir.unit();
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* evt) {
    G4ThreeVector x, v;
    if (genSurface.IsIsotropic()) {
        GenerateOnSphere(x, v);
    } else {
        v = genSurface.Axis();
        const G4double a = (2.0 * G4UniformRand() - 1.0) * genSurface.HalfU();
        const G4double b = (2.0 * G4UniformRand() - 1.0) * genSurface.HalfV();
        x = center + genSurface.U() * a + genSurface.V() * b - v * genSurface.Standoff();
    }
    ParticleInfo info = flux->GenerateParticle();

    particleGun->SetParticleDefinition(info.def);
    particleGun->SetParticleEnergy(info.energy);
    particleGun->SetParticlePosition(x);
    particleGun->SetParticleMomentumDirection(v);
    particleGun->SetParticleTime(0.0 * ns);
    particleGun->GeneratePrimaryVertex(evt);

    if (auto* ea = dynamic_cast<EventAction*>(G4EventManager::GetEventManager()->GetUserEventAction())) {
        PrimaryRec rec;
        rec.index = static_cast<int>(ea->primBuf.size());
        rec.pdg = info.pdg;
        rec.name = info.name;
        rec.E_MeV = info.energy / MeV;
        rec.dir = v;
        rec.pos_mm = x / mm;
        rec.t0_ns = 0.0;
        ea->primBuf.emplace_back(std::move(rec));
    }
}
