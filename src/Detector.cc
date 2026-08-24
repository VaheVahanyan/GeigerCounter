#include "Detector.hh"

using namespace Sizes;
using namespace Configuration;


Detector::Detector(G4LogicalVolume* detContLV, G4NistManager* nistMan) : nist(nistMan), detContainerLV(detContLV) {
    DefineMaterials();
    DefineVisual();
}

void Detector::DefineMaterials() {
    auto* elH = nist->FindOrBuildElement("H");
    auto* elC = nist->FindOrBuildElement("C");
    auto* elNe = nist->FindOrBuildElement("Ne");
    auto* elAr = nist->FindOrBuildElement("Ar");
    auto* elBr = nist->FindOrBuildElement("Br");

    G4Material* SiO2 = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    G4Material* polystyrene = nist->FindOrBuildMaterial("G4_POLYSTYRENE");

    // Geiger Tube
    geigerTubeMat = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    filamentMat = nist->FindOrBuildMaterial("G4_W");
    contactMat = nist->FindOrBuildMaterial("G4_Ni");

    {
        const G4int nNe = 494;
        const G4int nAr = 5;
        const G4int nBr2 = 1;
        const G4int nMolecules = nNe + nAr + nBr2;

        const G4double molarMass = (nNe * elNe->GetA() + nAr * elAr->GetA() + nBr2 * 2 * elBr->GetA()) / nMolecules;
        const G4double gasDensity = gasPressure * molarMass / (Avogadro * k_Boltzmann * gasTemperature);

        gasMat = new G4Material("CounterGas", gasDensity, 3, kStateGas, gasTemperature, gasPressure);
        gasMat->AddElementByNumberOfAtoms(elNe, nNe);
        gasMat->AddElementByNumberOfAtoms(elAr, nAr);
        gasMat->AddElementByNumberOfAtoms(elBr, nBr2 * 2);
    }

    // Box
    G4double foamDensity = 0.02 * g / cm3;
    boxMat = new G4Material("Styrofoam", foamDensity, polystyrene);

    // World / Vacuum
    galacticMat = nist->FindOrBuildMaterial("G4_Galactic");
    airMat = nist->FindOrBuildMaterial("G4_Air");

    // Filter
    filterMat = nist->FindOrBuildMaterial("G4_Al");

    // Retainer
    retainerMat = nist->FindOrBuildMaterial("G4_Al");

    // Board
    {
        auto* Epoxy = new G4Material("Epoxy", 1.2 * g / cm3, 2);
        Epoxy->AddElement(elH, 2);
        Epoxy->AddElement(elC, 2);

        boardMat = new G4Material("FiberglassLaminate", 1.86 * g / cm3, 2);
        boardMat->AddMaterial(Epoxy, 0.472);
        boardMat->AddMaterial(SiO2, 0.528);
    }
}


void Detector::DefineVisual() {
    visBox = new G4VisAttributes(G4Color(1.0, 1.0, 1.0));
    visBox->SetForceSolid(true);
    visLid = new G4VisAttributes(G4VisAttributes::GetInvisible());
    visLid->SetForceSolid(false);
    visGeigerTube = new G4VisAttributes(G4Color(0.3125, 0.235, 0.0));
    visGeigerTube->SetForceSolid(true);
    visGas = new G4VisAttributes(G4Color(0.0, 1.0, 1.0, 0.5));
    visGas->SetForceSolid(true);
    visFilament = new G4VisAttributes(G4Color(0.71, 0.67, 0.62));
    visFilament->SetForceSolid(true);
    visContact = new G4VisAttributes(G4Color(0.66, 0.57, 0.0));
    visContact->SetForceSolid(true);
    visFilter = new G4VisAttributes(G4Color(0.65, 0.65, 0.65));
    visFilter->SetForceSolid(true);
    visBoard = new G4VisAttributes(G4Color(1.0, 1.0, 0.0));
    visBoard->SetForceSolid(true);
    visRetainer = new G4VisAttributes(G4Color(0.65, 0.3, 0.05));
    visRetainer->SetForceSolid(true);
}


void Detector::Construct() {
    ConstructBox();
    ConstructGeigerTube();
    ConstructRetainer();
    ConstructFilter();
    ConstructBoard();
}

std::vector<G4LogicalVolume*> Detector::GetSensitiveLV() const {
    return {
        gasLV
    };
}


void Detector::ConstructBox() {
    G4VSolid* boxBase = new G4Box("BoxBase", Box::length, Box::width, Box::height);

    G4VSolid* sub1 = new G4Box("BoxSub1", Box::length - 2 * Box::thickness, Box::cavity1Width, Box::height);
    G4VSolid* boxBase1 = new G4SubtractionSolid("BoxBase1", boxBase, sub1, nullptr,
                                                G4ThreeVector(0, Box::width - Box::thickness * 2 - Box::cavity1Width,
                                                              Box::thickness * 2));

    G4VSolid* sub2 = new G4Box("BoxSub2", Box::cavity2Length, Box::cavity2Width * 2, Box::height);
    G4VSolid* boxBase2 = new G4SubtractionSolid("BoxBase2", boxBase1, sub2, nullptr,
                                                G4ThreeVector(0, Box::width - (Box::thickness + Box::cavity1Width) * 2 -
                                                              Box::cavity2Width, Box::thickness * 2));

    G4VSolid* sub3 = new G4Box("BoxSub3", Box::cavity3Length, Box::cavity3Width, Box::height);
    G4VSolid* boxBase3 = new G4SubtractionSolid("BoxBase3", boxBase2, sub3, nullptr,
                                                G4ThreeVector(0, Box::width - (Box::thickness + Box::cavity1Width +
                                                                  Box::cavity2Width) * 2 - Box::cavity3Width,
                                                              Box::thickness * 2));

    G4VSolid* sub4 = new G4Box("BoxSub4", Box::cavity4Length, Box::cavity4Width, Box::height);
    G4VSolid* boxBase4 = new G4SubtractionSolid("BoxBase4", boxBase3, sub4, nullptr,
                                                G4ThreeVector(-Box::length + Box::thickness * 2 + Box::cavity4Length,
                                                              -Box::width + Box::thickness * 2 + Box::cavity4Width,
                                                              Box::thickness * 2));

    G4VSolid* sub5 = new G4Box("BoxSub5", Box::cavity5Length1,
                               Box::cavity4Width + Box::cavity2Width + Box::cavity3Width, Box::height);
    G4VSolid* boxBase5 = new G4SubtractionSolid("BoxBase5", boxBase4, sub5, nullptr,
                                                G4ThreeVector(Box::length - (Box::thickness + Box::cavity5Length2 +
                                                                  Box::cavity5Length3) * 2 -
                                                              Box::cavity5Length1,
                                                              -Box::width + Box::thickness * 2 + Box::cavity3Width +
                                                              Box::cavity2Width + Box::cavity4Width,
                                                              Box::thickness * 2));

    G4VSolid* sub6 = new G4Box("BoxSub6", Box::cavity5Length2, Box::cavity4Width, Box::height);
    G4VSolid* boxBase6 = new G4SubtractionSolid("BoxBase6", boxBase5, sub6, nullptr,
                                                G4ThreeVector(Box::length - (Box::thickness + Box::cavity5Length2) * 2 -
                                                              Box::cavity5Length2,
                                                              -Box::width + (Box::thickness + Box::cavity5Width1) * 2 +
                                                              Box::cavity4Width, Box::thickness * 2));

    G4VSolid* sub7 = new G4Box("BoxSub7", Box::cavity5Length2, Box::cavity5Width2, Box::height);
    G4VSolid* box = new G4SubtractionSolid("Box", boxBase6, sub7, nullptr,
                                           G4ThreeVector(Box::length - (Box::thickness + Box::cavity5Length3) * 2 -
                                                         Box::cavity5Length2,
                                                         -Box::width + (Box::thickness + Box::cavity5Width1) * 2 +
                                                         Box::cavity5Width2, Box::thickness * 2));

    boxLV = new G4LogicalVolume(box, boxMat, "BoxLV");
    boxLV->SetVisAttributes(visBox);
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -Box::thickness), boxLV, "BoxPVP", detContainerLV, false, 0,
                      true);

    if (withLid) {
        G4VSolid* lid = new G4Box("Lid", Box::length, Box::width, Box::thickness);
        lidLV = new G4LogicalVolume(lid, boxMat, "LidLV");
        lidLV->SetVisAttributes(visLid);
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, Box::height), lidLV, "LidPVP", detContainerLV, false, 0, true);
    }
}

void Detector::ConstructGeigerTube() {
    G4VSolid* tube = new G4Tubs("GeigerTube", 0.0, GeigerCounter::radius, GeigerCounter::length, 0, 360 * deg);
    auto* rotMat = new G4RotationMatrix();
    rotMat->rotateY(270 * deg);
    rotMat->rotateZ(90 * deg);

    geigerTubeLV = new G4LogicalVolume(tube, geigerTubeMat, "GeigerTubeLV");
    geigerTubeLV->SetVisAttributes(visGeigerTube);

    G4VSolid* gas = new G4Tubs("Gas", GeigerCounter::filamentRadius, GeigerCounter::radius - GeigerCounter::thickness,
                               GeigerCounter::length - GeigerCounter::thickness, 0, 360 * deg);
    gasLV = new G4LogicalVolume(gas, gasMat, "GasLV");
    gasLV->SetVisAttributes(visGas);

    G4VSolid* filament = new G4Tubs("Filament", 0.0, GeigerCounter::filamentRadius,
                                    GeigerCounter::length - GeigerCounter::thickness, 0, 360 * deg);
    filamentLV = new G4LogicalVolume(filament, filamentMat, "FilamentLV");
    filamentLV->SetVisAttributes(visFilament);

    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), gasLV, "GasPVP", geigerTubeLV, false, 0, true);
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), filamentLV, "FilamentPVP", geigerTubeLV, false, 0, true);

    ConstructCounterRegion();

    new G4PVPlacement(rotMat,
                      G4ThreeVector(0, Box::width - (Box::thickness + Box::cavity1Gap) * 2 - GeigerCounter::radius, 0),
                      geigerTubeLV, "GeigerTubePVP", detContainerLV, false, 0, true);

    new G4PVPlacement(rotMat,
                      G4ThreeVector(0, Box::width - (Box::thickness + Box::cavity1Gap + Filter::width + Filter::gap * 2)
                                    * 2 - GeigerCounter::radius * 3, 0), geigerTubeLV, "GeigerTubePVP", detContainerLV,
                      false, 1, true);

    G4VSolid* contact = new G4Tubs("Contact", 0, GeigerCounter::contactRadius, GeigerCounter::contactLength, 0,
                                   360 * deg);
    contactLV = new G4LogicalVolume(contact, contactMat, "ContactLV");
    contactLV->SetVisAttributes(visContact);

    new G4PVPlacement(rotMat,
                      G4ThreeVector(GeigerCounter::length + GeigerCounter::contactLength,
                                    Box::width - (Box::thickness + Box::cavity1Gap) * 2 - GeigerCounter::radius, 0),
                      contactLV, "ContactPVP", detContainerLV, false, 0, true);

    new G4PVPlacement(rotMat,
                      G4ThreeVector(GeigerCounter::length + GeigerCounter::contactLength,
                                    Box::width - (Box::thickness + Box::cavity1Gap + Filter::width + Filter::gap * 2)
                                    * 2 - GeigerCounter::radius * 3, 0), contactLV, "ContactPVP", detContainerLV,
                      false, 1, true);

    new G4PVPlacement(rotMat,
                      G4ThreeVector(-(GeigerCounter::length + GeigerCounter::contactLength),
                                    Box::width - (Box::thickness + Box::cavity1Gap) * 2 - GeigerCounter::radius, 0),
                      contactLV, "ContactPVP", detContainerLV, false, 2, true);

    new G4PVPlacement(rotMat,
                      G4ThreeVector(-(GeigerCounter::length + GeigerCounter::contactLength),
                                    Box::width - (Box::thickness + Box::cavity1Gap + Filter::width + Filter::gap * 2)
                                    * 2 - GeigerCounter::radius * 3, 0), contactLV, "ContactPVP", detContainerLV,
                      false, 3, true);
}

void Detector::ConstructFilter() {
    G4VSolid* filter = new G4Box("Filter", Filter::length, Filter::width, Filter::height);
    filterLV = new G4LogicalVolume(filter, filterMat, "FilterLV");
    filterLV->SetVisAttributes(visFilter);

    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, Box::width - (Box::thickness + GeigerCounter::radius + Box::cavity1Gap +
                                        Filter::gap) * 2 - Filter::width, 0), filterLV, "FilterPVP", detContainerLV,
                      false, 0, true);
}

void Detector::ConstructRetainer() {
    G4VSolid* retainerBase = new G4Box("RetainerBase", Retainer::length, Retainer::width, Retainer::height);
    auto* rotMat = new G4RotationMatrix();
    rotMat->rotateY(270 * deg);
    rotMat->rotateZ(90 * deg);

    G4VSolid* sub1 = new G4Tubs("RetainerSub1", 0, GeigerCounter::radius, GeigerCounter::length, 0, 360 * deg);
    G4VSolid* base1 = new G4SubtractionSolid("RetainerBase1", retainerBase, sub1, rotMat,
                                             G4ThreeVector(0, Retainer::width - Box::cavity1Gap * 2 -
                                                           GeigerCounter::radius, 0));
    G4VSolid* base2 = new G4SubtractionSolid("RetainerBase2", base1, sub1, rotMat,
                                             G4ThreeVector(0, -Retainer::width + Box::cavity1Gap * 2 +
                                                           GeigerCounter::radius, 0));

    G4VSolid* sub2 = new G4Box("RetainerSub2", Filter::length, Filter::width, Filter::height);
    G4VSolid* retainer = new G4SubtractionSolid("Retainer", base2, sub2, nullptr, G4ThreeVector(0, 0, 0));

    retainerLV = new G4LogicalVolume(retainer, retainerMat, "RetainerLV");
    retainerLV->SetVisAttributes(visRetainer);
    new G4PVPlacement(nullptr, G4ThreeVector(GeigerCounter::length - Retainer::gap,
                                             Box::width - Box::cavity1Width - Box::thickness * 2, 0), retainerLV,
                      "RetainerPVP", detContainerLV, false, 0, true);
    new G4PVPlacement(nullptr, G4ThreeVector(-GeigerCounter::length + Retainer::gap,
                                             Box::width - Box::cavity1Width - Box::thickness * 2, 0), retainerLV,
                      "RetainerPVP", detContainerLV, false, 1, true);
}

void Detector::ConstructBoard() {
    G4VSolid* board = new G4Box("Board", Board::length, Board::width, Board::height);
    boardLV = new G4LogicalVolume(board, boardMat, "BoardLV");
    boardLV->SetVisAttributes(visBoard);

    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, Box::width - (Box::cavity1Width + Box::cavity2Width + Box::thickness) * 2 -
                                    Box::cavity3Width, 0), boardLV, "BoardPVP", detContainerLV, false, 0, true);
}

void Detector::ConstructCounterRegion() {
    counterRegion = G4RegionStore::GetInstance()->FindOrCreateRegion("CounterRegion");
    counterRegion->AddRootLogicalVolume(geigerTubeLV);

    auto* cuts = new G4ProductionCuts();
    cuts->SetProductionCut(counterRangeCut);
    counterRegion->SetProductionCuts(cuts);
}
