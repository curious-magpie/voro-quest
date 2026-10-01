#include "app/display.h"

#include <numeric>

namespace
{

// Translucent enough to see a cloud through, opaque enough to read as a box.
constexpr float kDomainAlpha = 0.2f;

} // namespace

void Display::refresh(const RvdSession &session)
{
    const size_t count = session.seeds().size();
    if (seeds_id_ != 0 && count == seeds_count_)
        // The same number of seeds: the viewer takes their new positions.
        return;

    seeds_id_ = next_id_++;
    seeds_count_ = count;
    seed_indices_.resize(count);
    std::iota(seed_indices_.begin(), seed_indices_.end(), 0u);
    ++overlay_revision_;
}

void apply_first_look(mesh_viewer::Viewer &viewer, const Display &display)
{
    mesh_viewer::Scene &scene = viewer.scene();

    if (mesh_viewer::Layer *box = scene.layer(display.domain_id()))
        box->alpha = kDomainAlpha;

    for (mesh_viewer::DebugOverlay::Item &item : scene.debug_overlay().items())
    {
        if (item.name == "seeds")
        {
            item.enabled = true;
        }
    }
}
