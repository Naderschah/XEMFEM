#include "cmdLineInteraction.h"

/*
https://stackoverflow.com/questions/865668/parsing-command-line-arguments-in-c
*/


namespace cli {

InputParser::InputParser(int& argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
        tokens.emplace_back(argv[i]);
    }
    // subcommand is only argv[1] if it exists and is not a flag
    if (argc >= 2) {
        std::string a1 = argv[1] ? std::string(argv[1]) : std::string();
        if (!a1.empty() && a1[0] != '-') {
            cmd_token = a1;
        }
    }
}


bool InputParser::has(const std::string& opt) const {
  return std::find(tokens.begin(), tokens.end(), opt) != tokens.end();
}

std::optional<std::string> InputParser::get(const std::string& opt) const {
  auto it = std::find(tokens.begin(), tokens.end(), opt);
  if (it != tokens.end() && ++it != tokens.end())
    return *it;
  return std::nullopt;
}
std::optional<std::string> InputParser::subcommand() const {
    return cmd_token;
}

// Remove the first positional token that does not start with '-' from argv/argc,
// so downstream parsers see the same flags as before.
void InputParser::strip_subcommand(int& argc, char** argv)
{
    if (argc < 2) return;
    const char* s = argv[1];
    if (!s || s[0] == '-') {
        // No subcommand present
        return;
    }
    // Remove argv[1]
    for (int j = 1; j < argc - 1; ++j) {
        argv[j] = argv[j + 1];
    }
    --argc;
}

std::filesystem::path to_absolute(const std::filesystem::path& p) {
  std::filesystem::path abs = p.is_absolute() ? p : std::filesystem::absolute(p);
  std::error_code ec;
  std::filesystem::path norm = std::filesystem::weakly_canonical(abs, ec);
  return ec ? abs : norm;
}

void print_usage(const char* prog) {
  std::cerr
    << "Usage:\n"
    << "  " << prog << " [sim|plot|metrics|interpolate] [options]\n\n"
    << "Notes:\n"
    << "  - If no subcommand is provided, 'sim' is assumed.\n\n"
    << "Subcommands:\n"
    << "  sim          Run simulation (single/sweep/opt determined by config)\n"
    << "  metrics      Compute optimization metrics only (uses config)\n"
    << "  plot         Generate plots (forwards args to plotter)\n"
    << "  interpolate  Generate Grid interpolation of mesh output\n\n"
    << "Common options (sim/metrics):\n"
    << "  -c, --config   Path to config file\n"
    << "  -h, --help     Show this help\n\n";
}

} // namespace cli
