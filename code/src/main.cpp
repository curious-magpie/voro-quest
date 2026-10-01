// rvd -- Voronoi and restricted Voronoi diagrams you can look inside.
//
//   rvd [--seeds N] [--random-seed S]
//
// Milestone M0 of ../design_document/DESIGN.tex: a box, a uniform cloud of
// seeds inside it, and the four windows the diagrams will live in.
//
// Mouse: drag to orbit, shift-drag or middle-drag to pan, wheel to zoom.
// Keys:  P cutting plane, W surface mode, F frame the scene, Esc quit.
//
// The window, the drawing and the Inspect and View windows are mesh-viewer's.
// What is here is the host:
//
//   RvdSession   the box and the seeds, and nothing about drawing  (session/)
//   Display      the ids and channels the viewer knows them by     (app/)
//   RvdSource,   the two adapters the viewer reads once a frame    (app/)
//   RvdOverlay
//   the windows  Diagram, Optimise, Cell, Robustness               (app/)
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include "mesh_viewer/viewer.h"

#include "app/display.h"
#include "app/rvd_source.h"
#include "app/windows.h"
#include "session/rvd_session.h"

namespace
{

int usage(const char *argv0)
{
    std::fprintf(
        stderr,
        "usage: %s [--seeds N] [--random-seed S]\n"
        "  --seeds N        how many uniform seeds to draw in the box\n"
        "  --random-seed S  what to draw them from; the same S gives the\n"
        "                   same seeds on every machine\n",
        argv0);
    return 1;
}

// A whole non-negative number, or false. strtoull alone accepts "12abc" and
// "-3", neither of which is a count.
bool parse_count(const char *text, unsigned long long &out)
{
    if (!text || !*text || *text == '-')
        return false;
    char *end = nullptr;
    out = std::strtoull(text, &end, 10);
    return *end == '\0';
}

} // namespace

int main(int argc, char **argv)
{
    // What RVD works on, declared before the viewer so that it outlives the
    // viewer's GPU copies of it.
    RvdSession session;
    Display display;

    for (int i = 1; i < argc; ++i)
    {
        unsigned long long value = 0;
        const bool has_value = i + 1 < argc && parse_count(argv[i + 1], value);
        if (std::strcmp(argv[i], "--seeds") == 0 && has_value &&
            value <= 10'000'000)
            session.seed_params.count = int(value);
        else if (std::strcmp(argv[i], "--random-seed") == 0 && has_value)
            session.seed_params.random_seed = value;
        else
            return usage(argv[0]);
        ++i;
    }
    session.regenerate_seeds();
    display.refresh(session);

    mesh_viewer::Viewer viewer;
    std::string err;
    if (!viewer.init({"rvd"}, err))
    {
        std::fprintf(stderr, "rvd: %s\n", err.c_str());
        return 1;
    }

    RvdSource meshes(session, display);
    RvdOverlay overlay(display);
    viewer.set_mesh_source(&meshes);
    viewer.set_overlay_source(&overlay);

    DiagramWindow diagram(session, display);
    OptimiseWindow optimise;
    CellWindow cell;
    RobustnessWindow robustness;
    viewer.add_window("Diagram", diagram, DiagramWindow::kPlacement);
    viewer.add_window("Optimise", optimise, OptimiseWindow::kPlacement);
    viewer.add_window("Cell", cell, CellWindow::kPlacement);
    viewer.add_window("Robustness", robustness, RobustnessWindow::kPlacement);

    // Synced once here, before run() does it again, so that the layers and the
    // channel exist to be given their first look.
    viewer.sync();
    apply_first_look(viewer, display);

    return viewer.run();
}
