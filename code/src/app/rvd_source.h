// What the viewer is shown: the session's meshes and the display's channels,
// in the viewer's terms (design, section 8.1).
//
// The whole boundary between RVD and mesh-viewer is these two adapters, as
// document_source.h is for mesh-dev. Neither copies anything: a MeshView
// points into the session's own arrays, and an OverlayChannel into the
// display's. The viewer reads them once a frame and re-uploads only what a
// revision says has changed.
//
// In M0 there are two meshes, listed in this order on purpose:
//
//   0  the domain  the box. First, because the first mesh the viewer ever sees
//                  fixes its scene origin, and it takes that from a mesh's
//                  triangles -- a mesh without any would put it at zero.
//   1  the seeds   a *carrier*: positions and no triangles. It draws nothing
//                  itself; it holds the vertices the "seeds" channel indexes.
#pragma once

#include "mesh_viewer/mesh_source.h"

#include "app/display.h"
#include "session/rvd_session.h"

class RvdSource : public mesh_viewer::MeshSource
{
  public:
    RvdSource(const RvdSession &session, const Display &display)
        : session_(session), display_(display)
    {
    }

    size_t mesh_count() const override
    {
        return 2;
    }
    mesh_viewer::MeshView mesh(size_t i) const override;

  private:
    const RvdSession &session_;
    const Display &display_;
};

class RvdOverlay : public mesh_viewer::OverlaySource
{
  public:
    explicit RvdOverlay(const Display &display) : display_(display)
    {
    }

    uint64_t revision() const override
    {
        return display_.overlay_revision();
    }
    size_t channel_count() const override
    {
        return 1;
    }
    mesh_viewer::OverlayChannel channel(size_t i) const override;

  private:
    const Display &display_;
};
