#include "plotting_api.h"
#include "path_handler.h"
// -----------------------------------------------------------------------------
// API entry
// -----------------------------------------------------------------------------
int make_plot_api(Config cfg)
{
    const auto targets = targets_from_save_root(cfg);

    // Optional: keep your message behavior
    if (targets.size() == 1 && targets.front() == std::filesystem::path(cfg.save_path)) {
        std::cout << "Print Single Line" << std::endl;
    }

    for (const auto& run_dir : targets)
    {
        const std::filesystem::path pvd = run_dir / "Simulation" / "Simulation.pvd";
        _make_plots(pvd, cfg.debug.debug);
    }

    return 0;
}