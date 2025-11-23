## ELectrostatics with MFEM

Compile solver: 

cmake /path/to/src && make 

Run solver 

TPC --config /path/to/config.yaml --mesh /path/to/mesh

Generate Geometry:



### Simulation configuration



## Software Stack


MFEM: 
- Version 4.8.0 
- Highest currently available

Hypre:
- version 2.33.0 
- version geq 3.0.0 does not work with MFEM 
- Previous versions 
  - tried with 2.32.0 : HypreBoomer AMG reports 1 OpenMP thread regardless of set number. Either checks in a funny way, or bugged  

cuda_toolkit:
- Not yet implemented, will be required for GPU (cuda)
- Versions unknown 



# TODO

## Geometry and config files 

Need to think about how ot easily supply configs and source dirs

Need to implement neumann Robin and periodic boundary conditions 

Should think about axisymmetric

Also really need to fix Electric Field 


## General

use pragma stop with ifndef 

Wall boundary conditions 

Also implement misc boundary conditions

Need to fix indentation 

Gotta stop doing auto - I dont think this one will happen 


Geometry 
- Need a few generic functions for fragment and return objects in some reasonable way 
- Need a function for select boundary and return 