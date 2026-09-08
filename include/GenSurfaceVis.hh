#ifndef GENSURFACEVIS_HH
#define GENSURFACEVIS_HH

#include <G4VUserVisAction.hh>
#include <G4VVisManager.hh>
#include <G4VisExtent.hh>

#include "GenSurface.hh"


class GenSurfaceVis : public G4VUserVisAction {
public:
    void Draw() override;

    static G4VisExtent Extent();

private:
    static void DrawPayloadBox(G4VVisManager* vis);
    static void DrawRectangle(G4VVisManager* vis, const GenSurface& s);
    static void DrawSphere(G4VVisManager* vis, const GenSurface& s);
};

#endif //GENSURFACEVIS_HH
