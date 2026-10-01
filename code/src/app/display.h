// What the viewer is shown of a session, as of the last refresh: the id each
// mesh goes by, and the indices of the overlay channels.
//
// The viewer uploads a mesh's connectivity once, the first time it sees its id,
// and from then on takes new positions only; a host must never reuse an id
// (mesh_source.h). Most of what RVD shows keeps its shape -- the box always has
// eight corners and twelve triangles, so it keeps one id for good. The seeds
// are the exception, already in M0: their *number* can change, and a vertex
// buffer cannot grow under the id it was uploaded with, so a new count takes a
// new id. From M2 on, every generation's meshes take new ids the same way
// (design, section 8.4).
//
// The channel list only ever grows at its end. The viewer keeps a channel's
// colour and checkbox while the same name stays in the same slot, so a
// channel added by a later milestone goes after the ones already here.
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "mesh_viewer/viewer.h"

#include "session/rvd_session.h"

class Display
{
  public:
    // The box keeps its shape, so it keeps its id. 0 means "no mesh" to the
    // viewer, which is why the ids start at 1.
    static constexpr uint32_t kDomainId = 1;

    // Brings the ids and the channels up to date with `session`. Call it after
    // anything that changes the session and before the viewer next syncs --
    // which is to say, in the same deferred action.
    void refresh(const RvdSession &session);

    uint32_t domain_id() const
    {
        return kDomainId;
    }
    // 0 until the first refresh.
    uint32_t seeds_id() const
    {
        return seeds_id_;
    }

    // The "seeds" channel: every seed, as a point. The viewer draws a point
    // only through a channel, so this is how the cloud is seen at all.
    const std::vector<uint32_t> &seed_indices() const
    {
        return seed_indices_;
    }

    // Changes whenever any channel does, and is 1 or more after the first
    // refresh: 0 is what a fresh viewer has already "seen".
    uint64_t overlay_revision() const
    {
        return overlay_revision_;
    }

  private:
    uint32_t next_id_ = kDomainId + 1;
    uint32_t seeds_id_ = 0;
    size_t seeds_count_ = 0; // how many seeds seeds_id_ was issued for
    std::vector<uint32_t> seed_indices_;
    uint64_t overlay_revision_ = 0;
};

// How RVD's meshes look the first time they appear, set through the viewer's
// public scene rather than left at its defaults. The box is made translucent,
// so that the seeds inside it can be seen, and the seeds channel is switched
// on, since the viewer starts every channel hidden and here the seeds are the
// picture. Call once, after the first sync; the user's changes are never
// touched again.
void apply_first_look(mesh_viewer::Viewer &viewer, const Display &display);
