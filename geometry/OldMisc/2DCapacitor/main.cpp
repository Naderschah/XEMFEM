#include <gmsh.h>
#include <cmath>
#include <vector>
#include <utility>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <unordered_set>
#include <memory>
#include <filesystem>

#include "Config.h"
#include "embedded_config.h" 

int main(int argc, char *argv[]) {
  // config path 
  auto cfg = std::make_shared<const Config>(
    Config::LoadFromString(embedded_config::kConfigYaml)
  );

  gmsh::initialize();
  gmsh::model::add(cfg->geometry_id);

  // Background size and alignment
  double bg_W = 100.0, bg_H = 100.0;
  double bg_cx = 0.0,  bg_cy = 0.0;

  // Capacitor sizes
  double cp_W = 10.0, cp_H = 3.0;

  // Capacitor placement (independent of background)
  double y_center = 0.0;
  double gap      = 7.0;              // distance BETWEEN PLATE CENTERS

  double y_bot_c = y_center - 0.5 * gap;
  double y_top_c = y_center + 0.5 * gap;

  // Sanity: dielectric must have positive height => gap > cp_H
  if (gap <= cp_H) {
    throw std::runtime_error("gap must be > cp_H so the dielectric has positive height");
  }

  // Background
  int background = gmsh::model::occ::addRectangle(
    bg_cx - 0.5 * bg_W, bg_cy - 0.5 * bg_H, 0.0, bg_W, bg_H
  );
  double roundedRadius = 1.0;
  // Plates (centered at y_bot_c / y_top_c)
  int bottom_capacitor = gmsh::model::occ::addRectangle(
    -0.5 * cp_W, y_bot_c - 0.5 * cp_H, 0.0, cp_W, cp_H, -1, roundedRadius
  );
  int top_capacitor = gmsh::model::occ::addRectangle(
    -0.5 * cp_W, y_top_c - 0.5 * cp_H, 0.0, cp_W, cp_H, -1, roundedRadius
  );

  // Dielectric fills from top of bottom plate to bottom of top plate
  double diel_y0 = y_bot_c + 0.5 * cp_H;   // top face of bottom plate
  double diel_h  = gap - cp_H;             // (center gap) - total plate thickness (cp_H)
  int dielectric = gmsh::model::occ::addRectangle(
    -0.5 * cp_W, diel_y0, 0.0, cp_W, diel_h
  );

  // Fragment 
  gmsh::vectorpair object = {{2, background}};
  gmsh::vectorpair tools  = {
    {2, dielectric}, {2, bottom_capacitor}, {2, top_capacitor}
  };

  gmsh::vectorpair outDimTags;
  std::vector<gmsh::vectorpair> outMap;
  gmsh::model::occ::fragment(object, tools, outDimTags, outMap, -1, true, true );

  // Also need to cut the area
  object.clear();
  object.push_back(outMap[0][0]);

  tools.clear();
  tools.push_back(outMap[1][0]);   // dielectric
  tools.push_back(outMap[2][0]);   // bottom capacitor
  tools.push_back(outMap[3][0]);   // top capacitor
  std::vector<gmsh::vectorpair> outMapBckgrndCut;
  gmsh::model::occ::fragment(object, tools, outDimTags, outMapBckgrndCut, -1, true, false );

  gmsh::model::occ::synchronize();

 // ======================= Materials ==========================

  // Extract background 
  std::vector<int> backgroundSurf;
  for (const auto &dt : outMapBckgrndCut[0]) {      
    if (dt.first == 2) backgroundSurf.push_back(dt.second);
  }

  // Assign Material
  if (!backgroundSurf.empty()) {
    int physBackground = gmsh::model::addPhysicalGroup(
        2, backgroundSurf, cfg->materials.find("air")->second.id);
    gmsh::model::setPhysicalName(2, physBackground, "air");
  }

  // Dielectric 
  std::vector<int> dielectricSurf;
  for (const auto &dt : outMap[1]) {      
    if (dt.first == 2) dielectricSurf.push_back(dt.second);
  }

  // Assign Material
  if (!dielectricSurf.empty()) {
    int physDielectric = gmsh::model::addPhysicalGroup(
        2, dielectricSurf, cfg->materials.find("dielectric")->second.id);
    gmsh::model::setPhysicalName(2, physDielectric, "dielectric");
  }

  // ==================== Boundary Conditions ========================

  // Helper: get boundary curves of given surface(s)
  auto boundaryCurvesOf = [&](const std::vector<int> &surfTags, bool combined) {
    gmsh::vectorpair in, out;
    in.reserve(surfTags.size());
    for (int s : surfTags) in.emplace_back(2, s);

    gmsh::model::getBoundary(in, out,
                            /*combined=*/combined,
                            /*oriented=*/false,
                            /*recursive=*/false);

    std::vector<int> curves;
    curves.reserve(out.size());
    for (auto &dt : out)
      if (dt.first == 1) curves.push_back(dt.second);
    std::sort(curves.begin(), curves.end());
    curves.erase(std::unique(curves.begin(), curves.end()), curves.end());
    return curves;
  };

  // Capacitors 
  std::vector<int> topSurf;
  for (const auto &dt : outMap[2]) {      
    if (dt.first == 2) topSurf.push_back(dt.second);
  }
  std::vector<int> topCurves = boundaryCurvesOf(topSurf, /*combined=*/false);

  // Assign Boundary 
  if (!topCurves.empty()) {
    int physBCTop = gmsh::model::addPhysicalGroup(
        1, topCurves, cfg->boundaries.find("TopPlate")->second.bdr_id);
    gmsh::model::setPhysicalName(1, physBCTop, "TopPlate");
  }

  std::vector<int> bottomSurf;
  for (const auto &dt : outMap[3]) {      
    if (dt.first == 2) bottomSurf.push_back(dt.second);
  }
  std::vector<int> bottomCurves = boundaryCurvesOf(bottomSurf, /*combined=*/false);

  if (!bottomCurves.empty()) {
    int physBCBottom = gmsh::model::addPhysicalGroup(
        1, bottomCurves, cfg->boundaries.find("BottomPlate")->second.bdr_id);
    gmsh::model::setPhysicalName(1, physBCBottom, "BottomPlate");
  }

  // ==================== Outer Background Boundary ========================

  // Now get the *outer* background perimeter only (no internal edges)
  std::vector<int> outerCurves;
  outerCurves = boundaryCurvesOf(backgroundSurf, /*combined=*/true);

  int physOuterBG = gmsh::model::addPhysicalGroup(
        1, outerCurves, cfg->boundaries.find("OuterBoundary")->second.bdr_id);
  gmsh::model::setPhysicalName(1, physOuterBG, "OuterBoundary");

  // -================================ Finish Up ================================
  gmsh::option::setNumber("General.Terminal", 1);
  gmsh::option::setNumber("General.Verbosity", 5);

  gmsh::write(std::filesystem::path(cfg->mesh.path).replace_extension(".brep").string());

  gmsh::option::setNumber("Mesh.SaveAll", 0);
  gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
  gmsh::option::setNumber("Mesh.Optimize", 1);
  gmsh::option::setNumber("General.NumThreads", cfg->compute.threads.num);


  double characteristicLength = 0.5;
  gmsh::option::setNumber("Mesh.CharacteristicLengthMax", characteristicLength);

  gmsh::model::mesh::generate(2);

  gmsh::write(cfg->mesh.path); 
  std::cout << "Created mesh file\n";
  
  gmsh::finalize();
  return 0; 
}
