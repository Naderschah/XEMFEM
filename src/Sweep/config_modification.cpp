#include <string>     
#include <vector>     
#include <utility>    
#include <stdexcept>
#include "Config.h"
#include "config_modification.h"

static inline std::vector<std::string> split_path(const std::string &path)
{
    std::vector<std::string> parts;
    std::string current;
    for (char c : path)
    {
        if (c == '.')
        {
            if (!current.empty())
            {
                parts.push_back(current);
                current.clear();
            }
        }
        else
        {
            current.push_back(c);
        }
    }
    if (!current.empty()) { parts.push_back(current); }
    return parts;
}