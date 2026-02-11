#pragma once

#include <algorithm>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>
#include <iostream>

#include <mpi.h>
#include <streambuf>
#include <unordered_map>
#include <sstream>

namespace cli {

namespace fs = std::filesystem;

/// Lightweight command line argument parser
/// Example usage:
///     cli::InputParser args(argc, argv);
///     if (args.has("-h") || args.has("--help")) ...
///     auto config = args.get("-c").value_or("default.yaml");
///     auto cmd = args.subcommand().value_or("sim");
class InputParser {
public:
    InputParser(int& argc, char** argv);

    /// Returns true if a flag (e.g. "-h") is present
    bool has(const std::string& opt) const;

    /// Returns the value following a flag (e.g. "-c config.yaml")
    std::optional<std::string> get(const std::string& opt) const;

    /// Returns the first positional token that does not start with '-'
    /// (e.g. "sim", "plot", "metrics")
    std::optional<std::string> subcommand() const;

    /// Removes the first positional token that does not start with '-'
    /// from argv/argc so downstream parsers see only flags/options
    static void strip_subcommand(int& argc, char** argv);

private:
    std::vector<std::string> tokens;

    std::optional<std::string> cmd_token;
};

/// Converts a given path to an absolute, normalized path
fs::path to_absolute(const fs::path& p);

/// Prints standard usage info for XEMFEM
void print_usage(const char* prog);

} // namespace cli

