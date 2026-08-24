#include "Geometry.hh"

using namespace Sizes;
using namespace Configuration;


Geometry::Geometry() {
    nist = G4NistManager::Instance();
    zeroRot = new G4RotationMatrix(0, 0, 0);

    detContainerSize = G4ThreeVector(Box::length, Box::width, Box::height + Box::thickness);
    detContainerPos = G4ThreeVector(0, 0, 0);

    worldHalfSize = std::max({Box::length, Box::width, Box::height + Box::thickness}) * 3;

    detContVisAttr = new G4VisAttributes(G4VisAttributes::GetInvisible());
    detContVisAttr->SetForceSolid(false);
}


void Geometry::ConstructDetector() {
    detContainer = new G4Box("DetectorContainer", Box::length, Box::width, Box::height + Box::thickness);
    detContainerLV = new G4LogicalVolume(detContainer, worldMat, "DetectorContainerLV");
    detContainerPVPL = new G4PVPlacement(zeroRot, detContainerPos, detContainerLV, "DetectorContainerPVPL", worldLV,
                                         false, 0, true);
    detContainerLV->SetVisAttributes(detContVisAttr);

    detector = new Detector(detContainerLV, nist);
    detector->Construct();
    std::vector<G4LogicalVolume*> sensitiveLV = detector->GetSensitiveLV();
    geigerCounterLV = sensitiveLV.at(0);
}


G4VPhysicalVolume* Geometry::Construct() {
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::Clean();
    G4LogicalVolumeStore::Clean();
    G4SolidStore::Clean();

    G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(cutsEnergyMin, cutsEnergyMax);

    worldMat = nist->FindOrBuildMaterial("G4_Galactic");

    worldBox = new G4Box("World", worldHalfSize, worldHalfSize, worldHalfSize);
    worldLV = new G4LogicalVolume(worldBox, worldMat, "WorldLV");
    worldPVP = new G4PVPlacement(zeroRot, G4ThreeVector(0, 0, 0), worldLV, "WorldPVPL", nullptr, false, 0, false);
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    ConstructDetector();

    return worldPVP;
}

void Geometry::ConstructSDandField() {
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();

    auto* counterSD = new SensitiveDetector("CounterSD", 0, "GeigerCounter");
    sdManager->AddNewDetector(counterSD);
    geigerCounterLV->SetSensitiveDetector(counterSD);
}
