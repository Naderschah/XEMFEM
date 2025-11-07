## ELectrostatics with MFEM

Compile solver: 

cmake /path/to/src && make 

Run solver 

TPC --config /path/to/config.yaml --mesh /path/to/mesh

Generate Geometry:



### Simulation configuration



# TODO

## Geometry and config files 

Need to think about how ot easily supply configs and source dirs

Need to implement neumann Robin and periodic boundary conditions 

Should think about axisymmetric

Also really need to fix Electric Field 
- Rather fix it up a little its pure GPT at the moment  


## General

use pragma stop with ifndef 

Wall boundary conditions 

Also implement misc boundary conditions

Need to fix indentation 

Gotta stop doing auto 


Geometry 
- Need a few generic functions for fragment and return objects in some reasonable way 
- Need a function for select boundary and return 