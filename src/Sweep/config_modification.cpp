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

// -----------------------------------------------------------------------------
// Config field access by path
// -----------------------------------------------------------------------------
void set_cfg_value_from_string(Config &cfg, const std::string &path, const std::string &value)
{
    auto parts = split_path(path);
    if (parts.empty())
        throw std::runtime_error("Empty sweep path : " + path);

    // Boundaries may have any name 
    if (parts[0] == "boundaries")
    {
        if (parts.size() != 3)
            throw std::runtime_error("Boundary path must be specified as boundaries.BC_name.value");

        const std::string &bname = parts[1];
        const std::string &field = parts[2];

        auto it = cfg.boundaries.find(bname);
        if (it == cfg.boundaries.end())
            throw std::runtime_error("No such boundary: " + bname);

        Boundary &b = it->second;
        
        // set the value 
        if (field == "type")    { b.type   = value;                      return; }
        if (field == "value")   { b.value  = from_string<double>(value); return; }
        // err on bdr_id (or when they choose a different field)
        throw std::runtime_error("Can not sweep on: " + field);
    }

    // Same for materials
    if (parts[0] == "materials")
    {
        if (parts.size() != 3)
            throw std::runtime_error("Material path must be specified as boundaries.name.epsilon_r");

        const std::string &mname = parts[1];
        const std::string &field = parts[2];

        auto it = cfg.materials.find(mname);
        if (it == cfg.materials.end())
            throw std::runtime_error("No such material: " + mname);

        Material &m = it->second;

        if (field == "id")        { m.id        = from_string<int>(value);    return; }
        if (field == "epsilon_r") { m.epsilon_r = from_string<double>(value); return; }

        throw std::runtime_error("Unknown material field: " + field);
    }

    // Top-level / nested structs
    if (parts[0] == "solver")
    {
        if (parts.size() != 2)
            throw std::runtime_error("Solver path must be solver.<field>");

        const std::string &field = parts[1];
        SolverSettings &s = cfg.solver;

        if (field == "axisymmetric")                 { s.axisymmetric = from_string<bool>(value);           return; }
        if (field == "axisymmetric_r0_bd_attribute") { s.axisymmetric_r0_bd_attribute = from_string<int>(value); return; }
        if (field == "order")                        { s.order        = from_string<int>(value);            return; }
        if (field == "assembly_mode")                { s.assembly_mode = value;                             return; }
        if (field == "solver")                       { s.solver        = value;                             return; }
        if (field == "precond")                      { s.precond       = value;                             return; }
        if (field == "atol")                         { s.atol          = from_string<double>(value);        return; }
        if (field == "rtol")                         { s.rtol          = from_string<double>(value);        return; }
        if (field == "maxiter")                      { s.maxiter       = from_string<int>(value);           return; }
        if (field == "printlevel")                   { s.printlevel    = from_string<int>(value);           return; }

        throw std::runtime_error("Unknown solver field: " + field);
    }

    if (parts[0] == "debug")
    {
        if (parts.size() != 2)
            throw std::runtime_error("Debug path must be debug.<field>");

        const std::string &field = parts[1];
        DebugSettings &d = cfg.debug;

        if (field == "debug")   { d.debug   = from_string<bool>(value); return; }
        if (field == "dry_run") { d.dry_run = from_string<bool>(value); return; }

        throw std::runtime_error("Unknown debug field: " + field);
    }

    if (parts[0] == "mesh")
    {
        if (parts.size() != 2)
            throw std::runtime_error("Mesh path must be mesh.<field>");

        const std::string &field = parts[1];
        MeshSettings &m = cfg.mesh;

        if (field == "path") { m.path = value; return; }

        throw std::runtime_error("Unknown mesh field: " + field);
    }

    if (parts[0] == "compute")
    {
        if (parts.size() != 3)
            throw std::runtime_error("Compute path must be compute.<subgroup>.<field>");

        const std::string &sub   = parts[1];
        const std::string &field = parts[2];

        if (sub == "mpi")
        {
            MPISettings &m = cfg.compute.mpi;
            if (field == "enabled")                  { m.enabled  = from_string<bool>(value);  return; }
            if (field == "ranks")                    { m.ranks    = from_string<int>(value);   return; }
            if (field == "ranks_auto")               { m.ranks_auto = from_string<bool>(value);return; }
            if (field == "repartition_after_refine") { m.repartition_after_refine = from_string<bool>(value); return; }
            throw std::runtime_error("Unknown compute.mpi field: " + field);
        }
        if (sub == "threads")
        {
            ThreadsSettings &t = cfg.compute.threads;
            if (field == "enabled")  { t.enabled  = from_string<bool>(value);  return; }
            if (field == "num")      { t.num      = from_string<int>(value);   return; }
            if (field == "num_auto") { t.num_auto = from_string<bool>(value);  return; }
            if (field == "affinity") { t.affinity = value;                     return; }
            throw std::runtime_error("Unknown compute.threads field: " + field);
        }
        if (sub == "device")
        {
            DeviceRuntime &d = cfg.compute.device;
            if (field == "type")     { d.type    = value;                      return; }
            if (field == "id")       { d.id      = from_string<int>(value);    return; }
            if (field == "id_auto")  { d.id_auto = from_string<bool>(value);   return; }
            if (field == "per_rank") { d.per_rank = from_string<int>(value);   return; }
            throw std::runtime_error("Unknown compute.device field: " + field);
        }

        throw std::runtime_error("Unknown compute subgroup: " + sub);
    }

    // You can extend this with more root-level things if needed, or explicitly
    // reject paths you consider "not sensible" to sweep over.
    throw std::runtime_error("No setter implemented for path: " + path);
}
