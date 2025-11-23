// parallel/runtime.h
#pragma once

#include "Config.h"
#include "mfem.hpp"

namespace parallel
{
    void init_environment(const Config &cfg,
                          int &argc, char **&argv);
} // namespace parallel
