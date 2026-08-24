#ifndef SIZES_HH
#define SIZES_HH

#include <G4Types.hh>
#include <G4SystemOfUnits.hh>

namespace Sizes
{
    namespace GeigerCounter
    {
        const G4double length = 50 * mm;
        const G4double radius = 9 * mm;
        const G4double thickness = 0.1 * mm;

        const G4double filamentRadius = 0.05 * mm;

        const G4double contactRadius = 3 * mm;
        const G4double contactLength = 2.5 * mm;
    }

    namespace Filter
    {
        const G4double length = 49 * mm;
        const G4double width = 3.5 * mm;
        const G4double height = 10 * mm;

        const G4double gap = 0.5 * mm;
    }

    namespace Box
    {
        const G4double height = 20 * mm;
        const G4double length = 75 * mm;
        const G4double thickness = 5 * mm;

        const G4double cavityDepth = 15 * mm;

        const G4double cavity1Width = 24.5 * mm;
        const G4double cavity1Gap = 1 * mm;

        const G4double cavity2Length = 40 * mm;
        const G4double cavity2Width = 2.5 * mm;

        const G4double cavity3Length = 50 * mm;
        const G4double cavity3Width = 10 * mm;

        const G4double cavity4Length = 20 * mm;
        const G4double cavity4Width = 25 * mm;

        const G4double cavity5Length1 = 30 * mm;
        const G4double cavity5Length2 = 10 * mm;
        const G4double cavity5Length3 = 2.5 * mm;
        const G4double cavity5Width1 = 12.5 * mm;
        const G4double cavity5Width2 = 10 * mm;


        const G4double width = cavity1Width + cavity2Width * 2 + cavity3Width + cavity4Width + 2 * thickness;
    }

    namespace Retainer
    {
        const G4double length = 1 * mm;
        const G4double width = Box::cavity1Width;
        const G4double height = Box::cavityDepth;

        const G4double gap = 5 * mm;
    }

    namespace Board
    {
        const G4double length = Box::cavity3Length;
        const G4double width = 1 * mm;
        const G4double height = Box::cavityDepth;
    }
}

#endif //SIZES_HH
