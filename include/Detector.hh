#ifndef DETECTOR_HH
#define DETECTOR_HH

#include <G4VisAttributes.hh>
#include <G4SubtractionSolid.hh>
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4SolidStore.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4OpticalSurface.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4UnionSolid.hh>
#include <G4MultiUnion.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4PVParameterised.hh>
#include <G4MaterialPropertiesTable.hh>
#include <G4MaterialPropertyVector.hh>
#include <G4Version.hh>
#include <G4LossTableManager.hh>
#include <G4Element.hh>
#include <G4Box.hh>
#include <G4Cons.hh>
#include <G4Torus.hh>
#include <G4VSolid.hh>
#include <G4Trd.hh>
#include <G4Tubs.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4Region.hh>
#include <G4RegionStore.hh>
#include <G4ProductionCuts.hh>

#include "Sizes.hh"
#include "Configuration.hh"
#include "Utils.hh"


class Detector {
public:
    G4NistManager* nist;

    Detector(G4LogicalVolume*, G4NistManager*);
    ~Detector() = default;

    void Construct();

    [[nodiscard]] std::vector<G4LogicalVolume*> GetSensitiveLV() const;

private:
    void DefineMaterials();
    void DefineVisual();

    G4LogicalVolume* detContainerLV;

    G4VisAttributes* visBox{};
    G4VisAttributes* visLid{};
    G4VisAttributes* visGeigerTube{};
    G4VisAttributes* visGas{};
    G4VisAttributes* visFilament{};
    G4VisAttributes* visContact{};
    G4VisAttributes* visFilter{};
    G4VisAttributes* visRetainer{};
    G4VisAttributes* visBoard{};

    G4Material* boxMat{};
    G4Material* filterMat{};
    G4Material* retainerMat{};
    G4Material* geigerTubeMat{};
    G4Material* gasMat{};
    G4Material* filamentMat{};
    G4Material* contactMat{};
    G4Material* boardMat{};
    G4Material* galacticMat{};
    G4Material* airMat{};

    G4LogicalVolume* boxLV{};
    G4LogicalVolume* lidLV{};
    G4LogicalVolume* filterLV{};
    G4LogicalVolume* geigerTubeLV{};
    G4LogicalVolume* filamentLV{};
    G4LogicalVolume* contactLV{};
    G4LogicalVolume* gasLV{};
    G4LogicalVolume* boardLV{};
    G4LogicalVolume* retainerLV{};

    G4Region* counterRegion{};

    void ConstructBox();
    void ConstructGeigerTube();
    void ConstructCounterRegion();
    void ConstructRetainer();
    void ConstructFilter();
    void ConstructBoard();
};


#endif //DETECTOR_HH
