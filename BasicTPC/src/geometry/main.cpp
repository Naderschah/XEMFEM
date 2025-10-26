#include <gmsh.h>
#include <cmath>
#include <vector>
#include <utility>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <unordered_set>

// ------------ Debug flags ------------
constexpr bool INCLUDE_BACKGROUND = true;
constexpr bool INCLUDE_ANODE      = false;
constexpr bool INCLUDE_CATHODE    = false;

using DimTag = std::pair<int,int>;

//---------------------------- Debugging --------------------------------
// Print any surfaces whose 1D boundary isn't closed
static void reportOpenSurfaceLoops() {
  std::vector<DimTag> surfs;
  gmsh::model::getEntities(surfs, 2);

  int bad = 0;
  for (const auto &s : surfs) {
    // curves on this surface
    std::vector<DimTag> curves;
    gmsh::model::getBoundary({s}, curves,
                             /*combined=*/false,
                             /*oriented=*/true,
                             /*recursive=*/false);

    // count endpoints of all curves
    std::unordered_map<int,int> pointDegree;

    for (const auto &c : curves) {
      if (c.first != 1) continue;
      std::vector<DimTag> pts;
      gmsh::model::getBoundary({c}, pts,
                               /*combined=*/false,
                               /*oriented=*/true,
                               /*recursive=*/false);
      for (const auto &p : pts)
        if (p.first == 0) pointDegree[p.second]++; // end points
    }

    // endpoints in a closed loop come in pairs (even degree)
    bool open = false;
    for (const auto &kv : pointDegree) {
      if (kv.second % 2 != 0) { open = true; break; }
    }
    if (open) {
      bad++;
      std::cout << "[OPEN] surface tag " << s.second
                << " with " << curves.size()
                << " boundary curves (odd-degree point found)\n";
    }
  }

  if (bad == 0) std::cout << "All surface loops appear closed.\n";
}
// Collect unique surface tags bounding a volume
static std::unordered_set<int> surfaceSetOfVolume(int volTag) {
  std::vector<DimTag> bdr, arg = {{3, volTag}};
  gmsh::model::getBoundary(arg, bdr, /*combined=*/true, /*oriented=*/false, /*recursive=*/false);
  std::unordered_set<int> s;
  for (auto &dt : bdr) if (dt.first == 2) s.insert(dt.second);
  return s;
}

// Describe one surface: bbox, type, and which major volume(s) it belongs to
static void describeSurface(int surfTag,
                            const std::unordered_set<int> &bgSurfs,
                            const std::unordered_set<int> &anodeSurfs,
                            const std::unordered_set<int> &cathodeSurfs)
{
  // bbox
  double xmin, ymin, zmin, xmax, ymax, zmax;
  gmsh::model::getBoundingBox(2, surfTag, xmin, ymin, zmin, xmax, ymax, zmax);
  const double cx = 0.5 * (xmin + xmax);
  const double cy = 0.5 * (ymin + ymax);
  const double cz = 0.5 * (zmin + zmax);

  // type (output parameter)
  std::string typ;
  gmsh::model::getType(2, surfTag, typ); // e.g. "Plane Surface", "Cylinder"

  // parents (volumes that own this surface)
  std::vector<int> upDims, upTags;
  gmsh::model::getAdjacencies(2, surfTag, upDims, upTags);
  std::vector<int> parents;
  for (size_t i = 0; i < upDims.size(); ++i)
    if (upDims[i] == 3) parents.push_back(upTags[i]);

  std::cout << "Surface " << surfTag
            << "  type=" << typ
            << "  bbox=[(" << xmin << "," << ymin << "," << zmin << ")→("
            << xmax << "," << ymax << "," << zmax << ")]"
            << "  center=(" << cx << "," << cy << "," << cz << ")\n";

  std::cout << "  Belongs to volumes: ";
  for (int v : parents) std::cout << v << " ";
  std::cout << "\n";

  std::cout << "  In sets: "
            << (bgSurfs.count(surfTag) ? "Background " : "")
            << (anodeSurfs.count(surfTag) ? "Anode " : "")
            << (cathodeSurfs.count(surfTag) ? "Cathode " : "")
            << "\n";
}

// Convenience: print summary for a specific tag (e.g. 8)
static void debugSurfaceTag(int surfTag, int background_vol, int anodeVol, int cathodeVol) {
  auto bgSet      = surfaceSetOfVolume(background_vol);
  auto anodeSet   = surfaceSetOfVolume(anodeVol);
  auto cathodeSet = surfaceSetOfVolume(cathodeVol);
  describeSurface(surfTag, bgSet, anodeSet, cathodeSet);
}
// helper: join tags into comma-separated list
static std::string joinTags(const std::vector<DimTag> &vec) {
  std::ostringstream ss;
  for (size_t i = 0; i < vec.size(); ++i) {
    if (i > 0) ss << ", ";
    ss << vec[i].second;
  }
  return ss.str();
}

// concise debug wrapper for gmsh::model::occ::fragment
static void debugFragment(const std::string &title,
                          const std::vector<DimTag> &objects,
                          const std::vector<DimTag> &tools,
                          std::vector<DimTag> &newEntities,
                          std::vector<std::vector<DimTag>> &provenanceByInput,
                          int tag = -1,
                          bool removeObject = true,
                          bool removeTool = false)
{
  std::cout << "\n=== " << title << " ===\n";
  std::cout << "Objects (" << objects.size() << "): " << joinTags(objects);
  if (removeObject) std::cout << "  [marked for removal]";
  std::cout << '\n';

  std::cout << "Tools (" << tools.size() << "): " << joinTags(tools);
  if (removeTool) std::cout << "  [marked for removal]";
  std::cout << '\n';

  gmsh::model::occ::fragment(objects, tools,
                             newEntities, provenanceByInput,
                             tag, removeObject, removeTool);

  std::cout << "NewEntities (" << newEntities.size() << "): "
            << joinTags(newEntities) << '\n';

  for (size_t i = 0; i < provenanceByInput.size(); ++i)
    std::cout << "Provenance[" << i << "] ("
              << provenanceByInput[i].size() << "): "
              << joinTags(provenanceByInput[i]) << '\n';

  std::cout << "============================================================\n";
}

auto getSurfaces = [](int volTag) {
  std::vector<std::pair<int,int>> vol = {{3, volTag}};
  std::vector<std::pair<int,int>> surfs;
  gmsh::model::getBoundary(vol, surfs, /*combined=*/true, /*oriented=*/false, /*recursive=*/false);

  // Extract the surface tags only without the dim specifier
  std::vector<int> surfTags;
  surfTags.reserve(surfs.size());
  for (auto &dt : surfs)
    if (dt.first == 2) surfTags.push_back(dt.second);
  return surfTags;
};
// Remove any 2D faces that are not on the boundary of a 3D volume
static void removeOrphanSurfaces() {
  // Collect all 2D faces that actually bound some 3D volume
  std::vector<DimTag> vols;
  gmsh::model::getEntities(vols, 3);

  std::unordered_set<int> boundarySurfs;
  for (const auto &v : vols) {
    std::vector<DimTag> bdr;
    gmsh::model::getBoundary({v}, bdr, /*combined=*/true, /*oriented=*/false, /*recursive=*/false);
    for (const auto &dt : bdr) if (dt.first == 2) boundarySurfs.insert(dt.second);
  }

  // Find *all* 2D faces
  std::vector<DimTag> surfs;
  gmsh::model::getEntities(surfs, 2);

  // Any surface not in boundarySurfs is orphaned → remove it
  std::vector<DimTag> toRemove;
  toRemove.reserve(surfs.size());
  for (const auto &s : surfs) {
    if (!boundarySurfs.count(s.second)) toRemove.push_back(s);
  }

  if (!toRemove.empty()) {
    std::cout << "Removing orphan surfaces: ";
    for (auto &s : toRemove) std::cout << s.second << " ";
    std::cout << "\n";
    gmsh::model::occ::remove(toRemove, /*recursive=*/true);
  }
}

//---------------------------- AnodeWires --------------------------------

auto makeWire = [](double xc, double yc, double zc,
                   double length, double radius)
{
  return gmsh::model::occ::addCylinder(
      xc - length * 0.5, yc, zc,
      length, 0.0, 0.0,  // oriented along X
      radius, -1, 2.0 * M_PI);
};

// Create all wire cylinders for one electrode (anode or cathode).
// Returns vector of volume tags.
static std::vector<int> makeElectrodeWires(int n_wires,
                                           double wire_diameter,
                                           double inner_radius,
                                           double zBase)
{
  std::vector<int> wireTags;
  if (n_wires <= 0) return wireTags;

  const double rWire = wire_diameter * 0.5;
  const double R     = inner_radius;

  // Feasibility: need room for n circles of radius r with uniform surface gaps
  // c = (2R + 2r) / (n + 1) must be > 2r  =>  R > n*r
  if (R <= n_wires * rWire) {
    // Not enough room to place n wires with equal end/neighbor gaps
    return wireTags;
  }

  // Center-to-center spacing that yields equal end gap and inter-wire surface gap
  const double c = (2.0 * R + 2.0 * rWire) / (n_wires + 1.0);

  // Symmetric center positions with spacing c
  std::vector<double> yPositions;
  yPositions.reserve(n_wires);
  const double y0 = -0.5 * (n_wires - 1) * c;
  for (int i = 0; i < n_wires; ++i) yPositions.push_back(y0 + i * c);

  // Helper to make a cylinder along X
  auto makeWire = [](double xc, double yc, double zc, double length, double radius) {
    return gmsh::model::occ::addCylinder(xc - length * 0.5, yc, zc,
                                         length, 0.0, 0.0, radius,
                                         /*tag=*/-1);
  };

  // Compute chord length so each wire is fully contained inside the inner radius
  for (double y : yPositions) {
    // Safety: by construction |y| <= R - (c - 2r) - something; but keep guard
    double delta = std::max(0.0, std::abs(y) - rWire);
    double chord = 2.0 * std::sqrt(std::max(0.0, R * R - delta * delta));
    if (chord > 0.0) {
      int tag = makeWire(0.0, y, zBase, chord, rWire);
      wireTags.push_back(tag);
    }
  }

  return wireTags;
}

static int makeRingElectrode(double xCenter,
                             double yCenter,
                             double zBase,
                             double thickness,
                             double innerR,
                             double outerR)
{
  if(outerR <= innerR || innerR <= 0.0 || thickness <= 0.0)
    throw std::runtime_error("makeRingElectrode: invalid parameters");

  // Build outer & inner solid cylinders
  const int outerTag = gmsh::model::occ::addCylinder(xCenter, yCenter, zBase,
                                                     0.0, 0.0, thickness,
                                                     outerR, -1);
  const int innerTag = gmsh::model::occ::addCylinder(xCenter, yCenter, zBase,
                                                     0.0, 0.0, thickness,
                                                     innerR, -1);
  // Subtract inner from outer to make a ring (hollow cylinder)
  std::vector<DimTag> result;
  std::vector<std::vector<DimTag>> outMap; // unused, but required by overload
  gmsh::model::occ::cut({{3, outerTag}}, {{3, innerTag}}, result, outMap, /*tag*/ -1, /*removeObject=*/true, /*removeTool=*/true);

  if(result.empty())
    throw std::runtime_error("makeRingElectrode: OCC cut produced no volume");

  // Usually a single volume; if multiple, return the first (caller may inspect others if desired)
  return result.front().second;
}
// Create a full electrode: ring + all its wires merged into one volume.
// Returns the resulting merged (unioned) volume tag.
static int makeParallelWireElectrodeAssembly(double xCenter,
                         double yCenter,
                         double zBase,
                         double zBaseWires,
                         double ring_thickness,
                         double inner_radius,
                         double outer_radius,
                         int n_wires,
                         double wire_diameter)
{
  // --- Create ring shell ---
  int ringTag = makeRingElectrode(xCenter, yCenter, zBase,
                                  ring_thickness, inner_radius, outer_radius);

  // --- Create parallel wire grid ---
  std::vector<int> wireTags = makeElectrodeWires(
      n_wires, wire_diameter, inner_radius, zBaseWires);

  if (wireTags.empty()) return ringTag;

  // --- Fuse ring + wires ---
  using DimTag = std::pair<int,int>;
  std::vector<DimTag> out;
  std::vector<std::vector<DimTag>> map;
  std::vector<DimTag> wires;
  for (int t : wireTags) wires.emplace_back(3, t);

  gmsh::model::occ::fuse({{3, ringTag}}, wires, out, map,
                         /*tag=*/-1, /*removeObject=*/true, /*removeTool=*/true);

  if (out.empty())
    throw std::runtime_error("makeElectrode: fuse produced no volume");

  return out.front().second;
}

int main() {
  gmsh::initialize();
  gmsh::model::add("tpc_occ");

  // ---- Constants ----
  const bool USE_ANODE   = true; // set true to include anode
  const bool USE_CATHODE = true;  // set true to include cathode

  constexpr double DriftRegionHeight   = 1.0;
  constexpr int    n_wires             = 11;
  constexpr double wire_diameter       = 0.002;
  constexpr double inner_radius        = 0.2;
  constexpr double outer_radius        = 0.25;
  constexpr double anode_wire_height   = 0.95;
  constexpr double cathode_wire_height = 0.05;
  constexpr double ring_thickness      = 0.05;

  constexpr int anode_BC_index   = 1000;
  constexpr int cathode_BC_index = 1001;
  constexpr int TPC_Volume_index = 2001;

  const double H = DriftRegionHeight + 2*ring_thickness;
  const double R = outer_radius * 1.1;
  const double rWire = wire_diameter * 0.5;

  // ---------------------  Build Geometry --------------------- 
  // Outer Volume
  int background_vol = gmsh::model::occ::addCylinder(0, 0, 0, 
                                                     0, 0, H, 
                                                     R);

  std::vector<DimTag> tools;                    // what we pass to fragment
  std::vector<std::pair<std::string,int>> elec; // {name, tag} in the same order as 'tools'

  // Ring Electrodes with wires
  if (USE_ANODE) {
    const double anode_zBase       = anode_wire_height - ring_thickness;
    const double anode_wireCenterZ = anode_zBase + rWire + 0.005; // tiny nudge up
    int anodeRing = makeParallelWireElectrodeAssembly( 0.0, 0.0, anode_zBase, anode_wireCenterZ, ring_thickness, inner_radius, outer_radius, n_wires, wire_diameter);
    //int anodeRing =  makeRingElectrode(0.0, 0.0, anode_zBase, ring_thickness, inner_radius, outer_radius);
    tools.emplace_back(3, anodeRing);
    elec.push_back({"Anode", anodeRing});
  }
  
  // cathode (wires at top extent: z_max - rWire)
  if (USE_CATHODE) {
    const double cathode_zBase       = cathode_wire_height;
    const double cathode_wireCenterZ = cathode_zBase + ring_thickness - rWire - 0.005; // tiny nudge down
    int cathodeRing = makeParallelWireElectrodeAssembly(
        0.0, 0.0, cathode_zBase, cathode_wireCenterZ,
        ring_thickness, inner_radius, outer_radius,
        n_wires, wire_diameter);
    tools.emplace_back(3, cathodeRing);
    elec.push_back({"Cathode", elec.empty() ? tools.back().second : tools.back().second}); // keep order
  }


  // Remove from background Volume
  std::vector<DimTag> newEntities;
  std::vector<std::vector<DimTag>> provenanceByInput;
  
  if (!tools.empty()) {
    debugFragment("Partition background with rings",
                  {{3, background_vol}}, tools,
                  newEntities, provenanceByInput,
                  /*tag*/ -1, /*removeObject=*/true, /*removeTool=*/false);
  
    // Track new background tag (bucket 0 is the object list)
    if (!provenanceByInput.empty() && !provenanceByInput[0].empty())
      background_vol = provenanceByInput[0][0].second;
  
    // Update each electrode’s tag from its provenance bucket
    // (bucket i+1 corresponds to tools[i])
    for (size_t i = 0; i < elec.size(); ++i) {
      if (provenanceByInput.size() > i + 1 && !provenanceByInput[i + 1].empty())
        elec[i].second = provenanceByInput[i + 1][0].second;
    }
  } else {
    std::cout << "No electrodes selected; skipping fragment.\n";
  }

  // ---------------------  Synchronize --------------------- 
  gmsh::model::occ::removeAllDuplicates();
  removeOrphanSurfaces();
  gmsh::model::occ::synchronize();
  reportOpenSurfaceLoops();
  // ------------  Create Physical Groups & Save --------------------- 
  gmsh::model::addPhysicalGroup(3, {background_vol}, TPC_Volume_index);
  gmsh::model::setPhysicalName(3, TPC_Volume_index, "TPC_Volume");
  
  // helper to get surface tags from a volume
  auto getSurfaces = [](int volTag) {
    std::vector<DimTag> bdr;
    gmsh::model::getBoundary({{3, volTag}}, bdr, /*combined=*/true, /*oriented=*/false, /*recursive=*/false);
    std::vector<int> out; out.reserve(bdr.size());
    for (auto &dt : bdr) if (dt.first == 2) out.push_back(dt.second);
    return out;
  };
  
  // Add phys surfaces for whichever electrodes were included
  for (auto &e : elec) {
    const std::vector<int> s = getSurfaces(e.second);
    if (s.empty()) {
      std::cout << "[warn] No surfaces for " << e.first << " (tag " << e.second << ")\n";
      continue;
    }
    if (e.first == "Anode") {
      gmsh::model::addPhysicalGroup(2, s, anode_BC_index);
      gmsh::model::setPhysicalName(2, anode_BC_index, "Anode");
    } else if (e.first == "Cathode") {
      gmsh::model::addPhysicalGroup(2, s, cathode_BC_index);
      gmsh::model::setPhysicalName(2, cathode_BC_index, "Cathode");
    }
  }

  gmsh::option::setNumber("General.Terminal", 1);
  gmsh::option::setNumber("General.Verbosity", 5); // more detail

  // Pure Cad (does not include physical groups)
  gmsh::write("tpc_occ.brep");
  // Internal Rep (includes physical groups)
  // gmsh::write("tpc_occ.geo_unrolled");

  // -------------------- Meshing ---------------------
  gmsh::option::setNumber("Mesh.SaveAll", 0);
  gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
  gmsh::option::setNumber("Mesh.Optimize", 1);
  //debugSurfaceTag(/*surfTag=*/8, background_vol, anodeRing, cathodeRing);
  //debugSurfaceTag(/*surfTag=*/9, background_vol, anodeRing, cathodeRing);

  // We generate the Maximum characteristic length based on the smallest geometry
  double r = wire_diameter * 0.5;
  double c = (2.0*inner_radius + 2.0*r) / (n_wires + 1.0);
  double gap = std::max(0.0, c - 2.0*r);
  double lc_circ = 2.0 * M_PI * r / 20.0;
  double lc_rad  = r / 3.0;
  double lc_gap  = (gap > 0.0) ? gap / 2.0 : lc_circ; // fallback if degenerate
  double lcWire  = std::max(1e-6, std::min({lc_circ, lc_rad, lc_gap}));
  
  gmsh::option::setNumber("Mesh.CharacteristicLengthMax", lcWire);
  
  gmsh::model::mesh::generate(3); 

  gmsh::write("tpc_occ.msh");
  std::cout << "Created mesh file" << std::endl;
  gmsh::finalize();
}
