## Geometries

To generate a simulation mesh gmsh is used. Gmsh has its own scripting language that was explored in all the .geo.in files where CMake was used to substitute in variable expressions. This language works for simple 2D problems allbeit very rough, the syntax is limited and since there is a Cpp language for it its not worthwhile using. Nonetheless at the bottom are my notes on this language.

### Logic

Goemetries and solver must share parameters with one another. Mainly relating to surface and volume ID's to assign boundary conditions, to achieve this the solver and indeces are predefined in a config.yaml file. Geometry relevant fields are 
```config.yaml
materials: 
  Name:
    attr_id: int
    epsilon_r: float

boundaries:
  Name:
    bdr_id: int
    type: dirichlet
    value: float
``` 
Start your ID's high the used range incremenets with every CAD operation, so we want to be well outside the range (but also not obscenely large as we have to make a vector the size of max(bdr_id) and one of max(attr_id)).
Each material and boundary is iterated over and its respective parameters are set. Note that the ID is the only reference between the mesh and the solver. 

The only boundary condition currently implemented is dirichlet.


### Geometry in Cpp

The best source on the language [is the header file](https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/api/gmsh.h), its fairly short and there are also tutorials. Avoid Booleans as much as possible but for simplicity it should be fine to fragment objects out of the background volume. Do not use boolean difference but use boolean fragment which merges the surfaces of the fragmented volumes (which is important for the ID assignment)





### Geo Scripting Language

TODO Include cleaned up bits from README_FromCapacitor2D