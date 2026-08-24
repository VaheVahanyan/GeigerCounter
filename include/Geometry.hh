#ifndef GEOMETRY_HH
#define GEOMETRY_HH

#include <G4UserLimits.hh>
#include <G4SystemOfUnits.hh>
#include <G4SolidStore.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4PVParameterised.hh>
#include <G4Box.hh>
#include <G4GeometryManager.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VisAttributes.hh>
#include <G4ProductionCuts.hh>
#include <G4ProductionCutsTable.hh>
#include <G4SubtractionSolid.hh>
#include <G4UnionSolid.hh>
#include <G4VSolid.hh>
#include <G4PhysicalConstants.hh>
#include <G4SDManager.hh>
#include <utility>

#include "SensitiveDetector.hh"
#include "Detector.hh"
#include "Sizes.hh"
#include "Configuration.hh"


class Geometry : public G4VUserDetectorConstruction {
public:
    Geometry();

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

private:
    G4NistManager* nist;
    G4Material* worldMat;
    G4double worldHalfSize;
    G4Box* worldBox;
    G4LogicalVolume* worldLV;
    G4VPhysicalVolume* worldPVP;

    Detector* detector;
    G4ThreeVector detContainerSize;

    G4VSolid* detContainer;
    G4LogicalVolume* detContainerLV;
    G4VPhysicalVolume* detContainerPVPL;
    G4ThreeVector detContainerPos;

    G4LogicalVolume* geigerCounterLV;

    G4RotationMatrix* zeroRot;

    G4VisAttributes* detContVisAttr;

    void ConstructDetector();
};

#endif //GEOMETRY_HH
