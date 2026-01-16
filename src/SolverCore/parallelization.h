// parallel/runtime.h
#pragma once

#include "Config.h"
#include "mfem.hpp"

namespace parallel
{
    void init_environment(Config &cfg,
                          int &argc, char **&argv);
} // namespace parallel
