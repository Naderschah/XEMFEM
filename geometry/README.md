# Geometry Implementation in SALOME 

Geometry is implemented in CAD library SALOME, [website](https://www.salome-platform.org/?page_id=327), using their python TUI interface. 

Geometry may be generated on the command line using `./make_mesh.sh`.
Alternatively if new geometry is to be implemented one can specify `./make_mesh.sh gui` to start salome grapical and run the script (`build_geometry.py`) with Ctrl+t inside.  If you would like to overwrite the default path set an env parameter inside salome gui or specify the --config_path commandline flag for autogen (needs to be docker relative ie we mount the XEMFEM directory at /work). 

`make_mesh.sh` now builds a local packaged SALOME image on demand from the official archive. The defaults target `SALOME 9.15.0`. Override them with `SALOME_VERSION`, `SALOME_URL`, or `SALOME_IMAGE` if you need to pin a different package or tag.

Note that if one is debugging post geometry partitioning one must restart salome after each test as partitioning determinism changes due to a new study not being equivalent to resetting the internal software state.



## In case a geometry is saved after Sketching 


To open enter the path to the produced folder
```
from ModelAPI import ModelAPI_Session

session = ModelAPI_Session.get()
session.closeAll()
ok = session.load("/work/geometry/mesh/slice_142.50deg_after_sketch.shaper")
print("loaded:", ok)
```
