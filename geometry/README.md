# Geometry Implementation in SALOME 

Geometry is implemented in CAD library SALOME, [website](https://www.salome-platform.org/?page_id=327), using their python TUI interface. 

Geometry may be generated on the command line using `./make_mesh.sh`.
Alternatively if new geometry is to be implemented one can specify `./make_mesh.sh gui` to start salome grapical and run the script (`build_geometry.py`) with Ctrl+t inside.  If you would like to overwrite the default path set an env parameter inside salome gui or specify the --config_path commandline flag for autogen (needs to be docker relative ie we mount the XEMFEM directory at /work). 

Note that if one is debugging post geometry partitioning one must restart salome after each test as partitioning determinism changes due to a new study not being equivalent to resetting the internal software state.
