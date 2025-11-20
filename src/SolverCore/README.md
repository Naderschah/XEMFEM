## Electrostatics Solver Core


Implements the main simulation suite, everything else is just wrappers










### Boundary Condition Reference
From ChatGPT:

| **Type** | **Mathematical Condition** | **Physical Meaning** | **Implication if used on outer boundary** |
|-----------|-----------------------------|-----------------------|--------------------------------------------|
| **Dirichlet** | $$ \phi = \phi_D \quad \text{on } \Gamma_D $$ | Fixed potential (e.g., conducting surface or grounded boundary). | Encloses the domain in a conducting shell at constant potential — all field lines terminate there. Often used to mimic a distant grounded enclosure. |
| **Neumann** | $$ \varepsilon \nabla \phi \cdot \mathbf{n} = g \quad \text{on } \Gamma_N $$ | Prescribed normal electric flux or surface charge density. For $$ g = 0 $$ → insulated or symmetry plane. | No electric field crosses the outer boundary — acts as a perfect insulator or symmetry plane. Reflects field lines; not a true “open” boundary. |
| **Robin** | $$ \alpha \phi + \beta \varepsilon \nabla \phi \cdot \mathbf{n} = r \quad \text{on } \Gamma_R $$ | Mixed BC: models resistive/impedance surfaces or approximate “open” (decaying) boundaries. | Allows partial field penetration; for $$ \partial \phi / \partial n = -\phi/R $$ it mimics free-space decay — a practical “open domain” approximation. |
| **Periodic** | $$ \phi(\mathbf{x}_1) - \phi(\mathbf{x}_2) = \text{const} $$ | Field repeats periodically between opposite faces. | Simulates an infinite lattice of repeating cells; not an isolated structure — outer boundary effectively connects back to the opposite side. |
| **No BCs (default)** | $$ \varepsilon \nabla \phi \cdot \mathbf{n} = 0 \quad \text{on } \partial \Omega $$ | Homogeneous Neumann on all edges → insulated cavity; solution defined up to an arbitrary constant. | Implicitly assumes a closed, perfectly insulating box. Field lines cannot exit; potential is arbitrary up to a constant (singular system). |




### Debugging


##### SIGSEGV in Solver while checking boundaries 
This

```
[tower:02289] Signal: Segmentation fault (11)
[tower:02289] Signal code: Address not mapped (1)
[tower:02289] Failing at address: 0x2331e1
[tower:02289] [ 0] /lib/x86_64-linux-gnu/libc.so.6(+0x42520)[0x7f0d96a85520]
[tower:02289] [ 1] ./SOLVER(+0x198f1a)[0x55a809d85f1a]
[tower:02289] [ 2] ./SOLVER(+0x1bd7ef)[0x55a809daa7ef]
[tower:02289] [ 3] ./SOLVER(+0xb16c4)[0x55a809c9e6c4]
[tower:02289] [ 4] ./SOLVER(+0xb0dc1)[0x55a809c9ddc1]
[tower:02289] [ 5] ./SOLVER(+0xa54bf)[0x55a809c924bf]
[tower:02289] [ 6] /lib/x86_64-linux-gnu/libc.so.6(+0x29d90)[0x7f0d96a6cd90]
[tower:02289] [ 7] /lib/x86_64-linux-gnu/libc.so.6(__libc_start_main+0x80)[0x7f0d96a6ce40]
[tower:02289] [ 8] ./SOLVER(+0xa4845)[0x55a809c91845]
```

means you computed the BC lines seperate from the surface mesh, I encountered this when doing both a 1D meshing algo and NETGEN_1D2D which caused both netgen and the algo to create 1D meshes on boundaries while meshing the BC's as well instead of assigning them after the meshing is complete (NOTE: I am not sure this is the actual reason as I never explicitly inspected the mesh for this but it seems like the most reasonable conclusion). 

Extra Context: This arrose specifically when Meshing Geom is the union of LXe GXe PTFE and BC's where all faces have their own topologically distinct faces, and BCs have topological distinct boundaries. The latter point is the problem as they need to be identical in the topology. This is diagnosed when moving the compute prior to the BC's and seeing no boundaries defined in the resulting mesh files. 


##### Part of the geometry not meshing 

I encountered this when using a custom 1D algorithm to get rid of the no 1D algo set while using NETGEN 1D 2D, I chose a very small non adaptive local line length (1e-5) and I think it did not mesh due to the size inconsistencies of NETGEN2D and my local line logic. Sizes shoulbe matched in this case. TODO Update here how matching should be done  



