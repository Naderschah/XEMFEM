// Compile with: mpicxx -O2 -std=c++17 parmesh_to_paraview_with_MPIRank.cpp -I/home/felix/mfem-install/include -L/home/felix/mfem-install/lib -lmfem  -lHYPRE -fopenmp -lmetis -o parmesh_to_paraview

#include "mfem.hpp"

#include <mpi.h>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

namespace fs = std::filesystem;

static std::string RankFileName(const std::string &prefix, int rank)
{
   std::ostringstream os;
   os << prefix << "." << std::setw(6) << std::setfill('0') << rank;
   return os.str();
}

// Accept either:
//   sim_results/simulation_pmesh
// or:
//   sim_results/simulation_pmesh.000000
// and normalize to the shared prefix:
//   sim_results/simulation_pmesh
static std::string NormalizeInputPrefix(const std::string &arg)
{
   if (arg.size() > 7)
   {
      const std::string suffix = arg.substr(arg.size() - 7);
      if (suffix.size() == 7 &&
          suffix[0] == '.' &&
          std::isdigit(static_cast<unsigned char>(suffix[1])) &&
          std::isdigit(static_cast<unsigned char>(suffix[2])) &&
          std::isdigit(static_cast<unsigned char>(suffix[3])) &&
          std::isdigit(static_cast<unsigned char>(suffix[4])) &&
          std::isdigit(static_cast<unsigned char>(suffix[5])) &&
          std::isdigit(static_cast<unsigned char>(suffix[6])))
      {
         return arg.substr(0, arg.size() - 7);
      }
   }
   return arg;
}

// Split an output path like:
//   out/simulation_partition
// into:
//   prefix_path = "out"
//   collection_name = "simulation_partition"
static void SplitOutputPath(const std::string &out_arg,
                            std::string &prefix_path,
                            std::string &collection_name)
{
   fs::path p(out_arg);
   prefix_path = p.parent_path().string();
   collection_name = p.filename().string();

   if (collection_name.empty())
   {
      collection_name = "paraview_output";
   }
}

int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);

   int rank = 0, size = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if (argc != 3)
   {
      if (rank == 0)
      {
         std::cerr
            << "Usage:\n"
            << "  " << argv[0] << " <input_prefix|input_rank_file> <output_path>\n\n"
            << "Examples:\n"
            << "  " << argv[0]
            << " sim_results/simulation_pmesh out/simulation_partition\n"
            << "  " << argv[0]
            << " sim_results/simulation_pmesh.000000 out/simulation_partition\n";
      }
      MPI_Finalize();
      return 1;
   }

   const std::string input_prefix = NormalizeInputPrefix(argv[1]);
   const std::string input_file   = RankFileName(input_prefix, rank);

   std::ifstream mesh_in(input_file);
   if (!mesh_in)
   {
      std::cerr << "[rank " << rank << "] Failed to open input mesh file: "
                << input_file << '\n';
      MPI_Abort(MPI_COMM_WORLD, 2);
      return 2;
   }

   // MFEM parallel mesh reader: each rank reads its own file/stream.
   mfem::ParMesh pmesh(MPI_COMM_WORLD, mesh_in);

   // Piecewise-constant discontinuous scalar field, one value per element.
   mfem::L2_FECollection fec(0, pmesh.Dimension());
   mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
   mfem::ParGridFunction partition_id(&pfes);

   // Set every local element to the MPI rank.
   partition_id = static_cast<mfem::real_t>(rank);

   std::string prefix_path, collection_name;
   SplitOutputPath(argv[2], prefix_path, collection_name);

   mfem::ParaViewDataCollection dc(collection_name, &pmesh);
   if (!prefix_path.empty())
   {
      dc.SetPrefixPath(prefix_path);
   }

   // Binary VTU output.
   dc.SetDataFormat(mfem::VTKFormat::BINARY);
   dc.SetHighOrderOutput(false);   // safer/default for plain partition coloring
   dc.SetLevelsOfDetail(1);
   dc.SetCycle(0);
   dc.SetTime(0.0);

   dc.RegisterField("partition_id", &partition_id);
   dc.Save();

   if (rank == 0)
   {
      std::cout << "Wrote ParaView dataset: " << argv[2] << '\n';
      std::cout << "Open the generated .pvtu or .pvd in ParaView and color by "
                   "\"partition_id\".\n";
   }

   MPI_Finalize();
   return 0;
}