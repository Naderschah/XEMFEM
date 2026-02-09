#pragma once
#include <string>
#include <vector>
#include <utility>

// For bookkeeping of sweeps
struct RunRecord
{
    std::string run_dir_name;
    std::vector<std::pair<std::string, std::string>> params; // (label, value)
};

// (path, value) for replacement in config  
struct Assignment {
    std::string path; 
    std::string value;
};

// (YAML) Path helpers 
static inline std::vector<std::string> split_path(const std::string &path);

// -----------------------------------------------------------------------------
// String to type helpers
// -----------------------------------------------------------------------------
template<typename T>
T from_string(const std::string &s);

template<>
inline int from_string<int>(const std::string &s)
{
    return std::stoi(s);
}

template<>
inline double from_string<double>(const std::string &s)
{
    return std::stod(s);
}

template<>
inline bool from_string<bool>(const std::string &s)
{
    if (s == "1" || s == "true" || s == "True" || s == "TRUE")  return true;
    if (s == "0" || s == "false" || s == "False" || s == "FALSE") return false;
    throw std::runtime_error("Cannot parse bool from '" + s + "'");
}

template<>
inline std::string from_string<std::string>(const std::string &s)
{
    return s;
}

template<typename T>
std::string to_string_any(const T &v)
{
    return std::to_string(v);
}

template<>
inline std::string to_string_any<std::string>(const std::string &v)
{
    return v;
}

template<>
inline std::string to_string_any<bool>(const bool &v)
{
    return v ? "true" : "false";
}

// -----------------------------------------------------------------------------
// Config field access by path
// -----------------------------------------------------------------------------
