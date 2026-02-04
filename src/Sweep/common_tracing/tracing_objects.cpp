#include "tracing_objects.h"

void DumpElectronPathsCSV(const Config                          &cfg,
                          const std::vector<ElectronTraceResult> &out_results)
{
    if (!cfg.debug.dumpdata) { return; }

    std::cout << "[DEBUG] Dumping Electron Paths in CIV" << std::endl;
    namespace fs = std::filesystem;
    fs::path outdir(cfg.save_path);
    fs::path outfile = outdir / "electron_paths_debug.csv";
    std::ofstream ofs(outfile);
    ofs << "# id, step, r, z, exit_code\n";
    ofs << std::setprecision(17);

    for (std::size_t i = 0; i < out_results.size(); ++i)
    {
        const auto &res = out_results[i];
        const auto &pts = res.points;

        for (std::size_t k = 0; k < pts.size(); ++k)
        {
            const mfem::Vector &x = pts[k];
            const double r = x(0);
            const double z = x(1);

            ofs << i << "," << k << "," << r << "," << z << ","
                << static_cast<int>(res.exit_code) << "\n";
        }
    }
}
