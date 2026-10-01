# RVD

Voronoi and restricted Voronoi diagrams you can look inside. The design is
`../design_document/DESIGN.tex`; this directory is its code, built milestone by milestone.

**Where it stands: M0.** The build, the host and its four windows, a box domain and a
uniform cloud of seeds inside it. There is no diagram yet: that is M2, after the predicates
of M1.

## Building

```
cmake -S . -B build -G Ninja      # or without -G, for make
cmake --build build -j
cd build && ctest                 # the tests: no display needed
build/rvd                         # the app
build/rvd --seeds 50000 --random-seed 7
```

Two sibling projects are taken by path, and both default to the checkouts next to this one
(`~/phd/my_software/`):

| option | default | what for |
|---|---|---|
| `RVD_MESH_DEV_DIR` | `../../mesh-dev` | mesh-dev's `io/` and `mesh/`, compiled here as `rvd_mesh_dev`, and its test harness `tests/test_support.h` |
| `RVD_VIEWER_DIR` | `ext/mesh-viewer` if that submodule exists, else `../../mesh-viewer` | the viewer |
| `RVD_BUILD_VIEWER` | `ON` | `OFF` builds the core and the tests only, with no GLFW, GL or ImGui |
| `RVD_BUILD_TESTS` | `ON` | the test binaries |

The design has the viewer as a submodule in `ext/mesh-viewer`, as mesh-dev does. That needs
this to be a git repository first; until then the build uses the sibling checkout, and once
the submodule is there it is picked up with no other change.

## Using it

The box is drawn translucent so the seeds inside it show, and the seeds are the Inspect
window's `seeds` channel. The **Diagram** window sets how many seeds to draw and from which
random seed. **regenerate** draws them again, and **next seed** draws a different cloud. The
same random seed gives the same points, bit for bit, on every machine. **Optimise**, **Cell**
and **Robustness** are placeholders that say which milestone fills them.

The viewer's own controls: drag to orbit, shift-drag or middle-drag to pan, wheel to zoom,
**P** for the cutting plane, **W** for the surface mode, **F** to frame the scene, **Esc** to
quit.

## Layout

```
src/
  points/      PointSet (positions, weights, revision), the uniform sampler
  domain/      Box: corners, and its boundary wound outward
  session/     RvdSession: the domain, the seeds, their parameters (no GL)
  app/         the host: Display (ids and channels), the two adapters, the windows
  main.cpp
tests/
  test_domain.cpp   the box, checked by mesh-dev's SurfaceMesh too
  test_points.cpp   the point set, the sampler, the session
```

`rvd_core` is everything outside `app/` and `main.cpp`. The tests link only that.
