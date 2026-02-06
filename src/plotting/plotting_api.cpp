#include "plotting_api.h"
// -----------------------------------------------------------------------------
// API entry
// -----------------------------------------------------------------------------
int make_plot_api(Config cfg)
{
  // First check if we had multy or single run
  const std::filesystem::path root = cfg.save_path;
  bool single_run = false;
  // --- detect single run by presence of .msh in root ---
  for (const auto& e : std::filesystem::directory_iterator(root)) {
      if (e.is_regular_file() && e.path().extension() == ".msh") {
          single_run = true;
          break;
      }
  }
  if (single_run) {
    std::cout << "Print Single Line" << std::endl;
    std::filesystem::path pvd = root / "Simulation" / "Simulation.pvd";
    return _make_plots(pvd, cfg.debug.debug);
  }

  // --- multiple runs ---
  for (const auto& e : std::filesystem::directory_iterator(root)) {
    if (!e.is_directory()) continue;

    const std::string name = e.path().filename().string();
    if (name.rfind("run_", 0) != 0) continue;

    std::filesystem::path pvd = e.path() / "Simulation" / "Simulation.pvd";

    _make_plots(pvd, cfg.debug.debug);
  }
  return 0;
}