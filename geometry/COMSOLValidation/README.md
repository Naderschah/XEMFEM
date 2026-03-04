This folder holds configuration files and matched COMSOL files to validate the function of XEMFEM.

This is required as matching the more complex TPC geometry has prooved prohibitatively complex.

The examples include in 2D only (quicker to implement in both frameworks)

1. Parallel plate capacitor in volume – verifies uniform-field accuracy, fringing fields, and sensitivity to outer boundary conditions.

2. Eccentric coaxial cylinder – tests asymmetric gap handling, surface charge redistribution, and peak-field localization due to misalignment.

3. Two opposite-polarity wires – probes small-gap field amplification and resolution of near-contact electrostatic interactions.

4. Conducting wedge facing a plane – evaluates solver behavior near true corner singularities and field-scaling with distance from the edge.

To replicate the comparison you will want at least 30GB of hard drive space available, around 2-3 hours, of course COMSOL with the electrostatics license (im using 6.3), and something around 16GB of RAM. Your device will likely freeze during if your not well above this.

```bash 
for cfg_path in COMSOLValidation/*/; do ./geometry/make_mesh.sh tui --config_path "$cfg_path"; done
```

Then post processing
```bash
# Enter the XEMFEM shell
./run_shell
cd geometry
for cfg in COMSOLValidation/*/config.yaml; do python3 postprocess_mesh_soft.py -c "$cfg"; done
```
run the simulations (do `mpirun -np 5` if you dont want all core processes running but only 5, depending on CPU speed your device may not be usable for other tasks when running with all)
```bash
cd ../build
# Make if required
cmake ../src
make 
# Actually run
for dir in ../geometry/COMSOLValidation/*/; do mpirun ./XEMFEM -c "$dir/config.yaml"; done
```
Lastly produce the grid interpolations (this will take quite long)
```bash
for dir in ../geometry/COMSOLValidation/*/; do mpirun ./XEMFEM interpolate -c "$dir/config.yaml"; done
```
Now we also want to compare directly with COMSOL, this will generate nastran bdf files that COMSOL can actually import, note that this will convert the AMR mesh not the source mesh 
```bash
for dir in ../geometry/COMSOLValidation/*/; do python3 ../geometry/COMSOLValidation/createCOMSOLmesh.py -c "$dir/config.yaml"; done
```
Now produce the COMSOL comparison (it will get deleted if XEMFEM is run after due to the path its saved in), adjust export paths for your filesystem. The COMSOL mesh based results can be easily made, for the mesh specific imports I cant seem to find a way to keep the boundaries preselected, so import the mesh (Nastran format .bdf) tick import as linear otherwise it makes them curved elements, go to materials select them again, same for boundaries, and make sure COMSOL boundary condition *type and value* exactly match the XEMFEM config (for these validation configs `BC_Outer` is typically Dirichlet/Ground at 0 V, not Zero Charge). Remember to change the export name.

And to plot stay in the XEMFEM env and run the following

```bash
python -m pip install --no-cache-dir --force-reinstall \
  "numpy<2" \
  "matplotlib<4" \
  "h5py" \
  "pyvista"

# Strict mode (default):
# - requires COMSOL and XEMFEM meshgrid coordinates to match within --coord-tol
# - prints grid size + coordinate offset diagnostics
python3 ../geometry/COMSOLValidation/make_plots.py /work/sim_results/ --plot-dpi 600 --imshow-interpolation bilinear --axes-aspect equal
```
The plots will be placed in the respective directories.
