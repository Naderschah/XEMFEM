#include "plotting_api.h"
#include "path_handler.h"
// -----------------------------------------------------------------------------
// API entry
// -----------------------------------------------------------------------------
int make_plot_api(Config cfg)
{
    const auto targets = targets_from_save_root(cfg);

    for (const auto& run_dir : targets)
    {
        const std::filesystem::path pvd = run_dir / "Simulation" / "Simulation.pvd";
        _make_plots(pvd, cfg.debug.debug);
    }

    return 0;
}