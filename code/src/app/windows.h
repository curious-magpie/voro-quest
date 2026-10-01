// RVD's four windows (design, section 8.6). Each is a Panel, registered in
// main.cpp with Viewer::add_window, as mesh-dev's Model window is.
//
//   Diagram     the domain, the seeds, and (from M2) the diagram and its look
//   Optimise    Lloyd and L-BFGS (M5), transport (M7)
//   Cell        the selected cell (M2)
//   Robustness  predicate counters (M1), the switches and the validators (M2)
//
// In M0 only Diagram has controls. The other three are in place, and say what
// they will hold and when, so that the layout -- which ImGui remembers in
// imgui.ini -- is settled from the start.
//
// Which window a control belongs in follows mesh-dev's rule: if the viewer were
// deleted and this were a command-line program, would the control still have
// to exist? Then it is here, in one of these. Otherwise it is the viewer's.
#pragma once

#include "mesh_viewer/panel.h"

#include "app/display.h"
#include "session/rvd_session.h"

class DiagramWindow : public mesh_viewer::Panel
{
  public:
    static constexpr mesh_viewer::WindowPlacement kPlacement{12.0f, 12.0f,
                                                             380.0f, 250.0f};

    // Both borrowed; they must outlive the viewer's run().
    DiagramWindow(RvdSession &session, Display &display)
        : session_(session), display_(display)
    {
    }

    void ui(mesh_viewer::Viewer &viewer) override;

  private:
    // Queued through Viewer::defer, never run from ui(): drawing new seeds
    // changes a mesh, and the viewer's rule is that meshes change between
    // frames, not halfway through drawing the windows that describe them.
    void regenerate(mesh_viewer::Viewer &viewer);

    RvdSession &session_;
    Display &display_;
};

class OptimiseWindow : public mesh_viewer::Panel
{
  public:
    static constexpr mesh_viewer::WindowPlacement kPlacement{12.0f, 274.0f,
                                                             380.0f, 100.0f};

    void ui(mesh_viewer::Viewer &viewer) override;
};

class CellWindow : public mesh_viewer::Panel
{
  public:
    // Anchored to the right edge, under the viewer's View window.
    static constexpr mesh_viewer::WindowPlacement kPlacement{-312.0f, 364.0f,
                                                             300.0f, 100.0f};

    void ui(mesh_viewer::Viewer &viewer) override;
};

class RobustnessWindow : public mesh_viewer::Panel
{
  public:
    static constexpr mesh_viewer::WindowPlacement kPlacement{-312.0f, 476.0f,
                                                             300.0f, 100.0f};

    void ui(mesh_viewer::Viewer &viewer) override;
};
