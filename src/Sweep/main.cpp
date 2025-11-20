
#include <memory>
#include <filesystem>
#include <string>

#include "solver_api.h"
#include "Config.h"
#include "cmdLineParser.h"


int main(int argc, char *argv[])
{
  """
  Logic here should  be read the config 

  parse destination path and wheter or not to generate plots from cmd line?  


  Need to add basic sweep logic here 
  """
  cli::InputParser args(argc, argv);
  if (args.has("-h") || args.has("--help")) {
      cli::print_usage(argv[0]);
      return 0;
  }

  // config path 
  auto config_str_opt = args.get("-c");
  if (!config_str_opt) config_str_opt = args.get("--config");
  if (!config_str_opt) {
      std::cerr << "Error: missing required argument -c/--config\n";
      cli::print_usage(argv[0]);
      return 1;
  }
  auto config_path = cli::to_absolute(*config_str_opt);
  if (!std::filesystem::exists(config_path)) {
      std::cerr << "Error: config file not found: " << config_path << "\n";
      return 1;
  }

  auto cfg = std::make_shared<const Config>(
      Config::Load(config_path)
  );

  """
  Implement function that receives a config and sweep entries, 
  it sets the first value of its sweep entry and calls itself, with the config and remaining sweep entries
  If no more sweep entries remain, it computes 
  This should call every possible combination right?

  Also implement the suoptions in it I am not sure how to do the global optimization
  """
  do_sweep_over_param(base_conf, sweeps_remain) 
}