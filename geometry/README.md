# Geometry Implementation in SALOME 

Geometry is implemented in CAD library SALOME, [website](https://www.salome-platform.org/?page_id=327), using their python TUI interface. 

To generate the geometry and mesh it the python script build_geometry.py needs to be run from within the GUI, to do so simply hit Ctrl+t and select the script. 

Debugging note, Salome needs to be restarted in between tests if moving past the partition step, for some reason in Salome 9.9 the order of partition results objects changes on being run a second time, likely due to the python interpreter hanging on to some variables changing the global scope, this affects only the manual renaming indices. 

The easiest way to run SALOME is via docker
```
# Normal GPU 
docker run --rm \
  --user $(id -u):$(id -g) \
  -e DISPLAY=$DISPLAY \
  -e QT_X11_NO_MITSHM=1 \
  -e XAUTHORITY=/home/felix/.Xauthority \
  -v /tmp/.X11-unix:/tmp/.X11-unix:ro \
  -v $XAUTHORITY:/home/felix/.Xauthority:ro \
  -v /home/felix:/home/felix \
  --shm-size=2g --ipc=host --net=host \
  trophime/salome:9.9.0-focal salome

# Nvidia GPU
docker run --rm \
  --user $(id -u):$(id -g) \
  --device nvidia.com/gpu=all \
  -e NVIDIA_VISIBLE_DEVICES=all \
  -e NVIDIA_DRIVER_CAPABILITIES=all \
  -e DISPLAY=$DISPLAY \
  -e QT_X11_NO_MITSHM=1 \
  -e XAUTHORITY=/home/felix/.Xauthority \
  -v /tmp/.X11-unix:/tmp/.X11-unix:ro \
  -v $XAUTHORITY:/home/felix/.Xauthority:ro \
  -v /home/felix:/home/felix \
  --shm-size=2g --ipc=host --net=host \
  trophime/salome:9.9.0-focal salome
```
After meshing one has a .med file, this unfortunately needs converting via gmsh, as this is not available inside the GUI it is made available via the MFEM environment. 

From there the following commands must be run:

```
gmsh -0 mesh.med -format msh2 -o mesh22.msh
awk '/\$PhysicalNames/{flag=1;next}/\$EndPhysicalNames/{flag=0}flag' mesh22.msh \
  > phys_map.txt
python3 make_elements.py
# And then merge config and config_autogen.py
```
The first converts the mesh, the second creates a list of boundary and element names and their attribute ids the third generates the yaml file for XENONnT to register these. Voltages are assigned automatically and need adjusting, the Field Shaping rings are automatically computed by the internal circuit solver. 