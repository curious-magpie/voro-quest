#include "app/rvd_source.h"

mesh_viewer::MeshView RvdSource::mesh(size_t i) const
{
    mesh_viewer::MeshView v;

    if (i == 0)
    {
        // The box's shape never changes, only where its corners are, so it
        // keeps one id and its revision says when they moved. Closed and wound
        // outward, which the viewer's culling and its cutting-plane cap rely
        // on.
        const Box &box = session_.box();
        v.id = display_.domain_id();
        v.name = "box";
        v.revision = box.revision();
        v.positions = {box.corners().data(), box.corners().size()};
        v.triangles = {box.boundary_triangles().data(),
                       box.boundary_triangles().size()};
        v.closed = true;
        v.oriented = true;
        return v;
    }

    // The seeds, as positions only. Not closed: with no triangles there is no
    // cross-section to cap.
    const PointSet &seeds = session_.seeds();
    v.id = display_.seeds_id();
    v.name = "seeds";
    v.revision = seeds.revision();
    v.positions = seeds.positions();
    v.closed = false;
    v.oriented = false;
    return v;
}

mesh_viewer::OverlayChannel RvdOverlay::channel(size_t) const
{
    return {"seeds", mesh_viewer::OverlayPrim::Points, display_.seeds_id(),
            display_.seed_indices()};
}
