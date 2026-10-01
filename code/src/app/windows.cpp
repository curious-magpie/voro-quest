#include "app/windows.h"

#include <algorithm>

#include "imgui.h"

#include "mesh_viewer/viewer.h"

namespace
{

// Far more than the viewer draws comfortably as points, and a bound on what a
// slip of the keyboard can ask for.
constexpr int kMaxSeeds = 10'000'000;

} // namespace

// ---------------------------------------------------------------------------
// Diagram
// ---------------------------------------------------------------------------

void DiagramWindow::ui(mesh_viewer::Viewer &viewer)
{
    const Box &box = session_.box();
    const PointSet &seeds = session_.seeds();

    ImGui::TextDisabled("M0: a box and its seeds. The diagram arrives in M2.");
    ImGui::Separator();

    ImGui::Text("domain  box  [%g %g %g] .. [%g %g %g]", box.lo().x, box.lo().y,
                box.lo().z, box.hi().x, box.hi().y, box.hi().z);
    ImGui::Text("seeds   %zu, uniform", seeds.size());
    ImGui::TextDisabled("        revision %llu",
                        static_cast<unsigned long long>(seeds.revision()));
    ImGui::Separator();

    // The parameters are edited in place; nothing is drawn until a button asks.
    RvdSession::SeedParams &p = session_.seed_params;
    if (ImGui::InputInt("count", &p.count, 100, 1000))
        p.count = std::clamp(p.count, 0, kMaxSeeds);
    ImGui::InputScalar("random seed", ImGuiDataType_U64, &p.random_seed);

    if (ImGui::Button("regenerate"))
        regenerate(viewer);
    ImGui::SameLine();
    // The same seed always draws the same points, so "a different cloud" is a
    // different seed.
    if (ImGui::Button("next seed"))
    {
        ++p.random_seed;
        regenerate(viewer);
    }
}

void DiagramWindow::regenerate(mesh_viewer::Viewer &viewer)
{
    viewer.defer(
        [this]
        {
            session_.regenerate_seeds();
            // A new count needs a new id for the seeds, and a new channel; the
            // same count needs neither, and the revision carries the new
            // points.
            display_.refresh(session_);
        });
}

// ---------------------------------------------------------------------------
// The windows still to come
// ---------------------------------------------------------------------------

void OptimiseWindow::ui(mesh_viewer::Viewer &)
{
    ImGui::TextDisabled("Lloyd and L-BFGS on the CVT energy: M5.");
    ImGui::TextDisabled("Newton on the weights (transport): M7.");
}

void CellWindow::ui(mesh_viewer::Viewer &)
{
    ImGui::TextDisabled("The selected cell, its neighbours and");
    ImGui::TextDisabled("its replay: M2, with the first diagram.");
}

void RobustnessWindow::ui(mesh_viewer::Viewer &)
{
    ImGui::TextDisabled("Predicate counters: M1.");
    ImGui::TextDisabled("The switches and the validators: M2.");
}
