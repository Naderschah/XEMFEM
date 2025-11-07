## Electrostatics Solver







### Boundary Condition Reference
From ChatGPT:

| **Type** | **Mathematical Condition** | **Physical Meaning** | **Implication if used on outer boundary** |
|-----------|-----------------------------|-----------------------|--------------------------------------------|
| **Dirichlet** | $$ \phi = \phi_D \quad \text{on } \Gamma_D $$ | Fixed potential (e.g., conducting surface or grounded boundary). | Encloses the domain in a conducting shell at constant potential — all field lines terminate there. Often used to mimic a distant grounded enclosure. |
| **Neumann** | $$ \varepsilon \nabla \phi \cdot \mathbf{n} = g \quad \text{on } \Gamma_N $$ | Prescribed normal electric flux or surface charge density. For $$ g = 0 $$ → insulated or symmetry plane. | No electric field crosses the outer boundary — acts as a perfect insulator or symmetry plane. Reflects field lines; not a true “open” boundary. |
| **Robin** | $$ \alpha \phi + \beta \varepsilon \nabla \phi \cdot \mathbf{n} = r \quad \text{on } \Gamma_R $$ | Mixed BC: models resistive/impedance surfaces or approximate “open” (decaying) boundaries. | Allows partial field penetration; for $$ \partial \phi / \partial n = -\phi/R $$ it mimics free-space decay — a practical “open domain” approximation. |
| **Periodic** | $$ \phi(\mathbf{x}_1) - \phi(\mathbf{x}_2) = \text{const} $$ | Field repeats periodically between opposite faces. | Simulates an infinite lattice of repeating cells; not an isolated structure — outer boundary effectively connects back to the opposite side. |
| **No BCs (default)** | $$ \varepsilon \nabla \phi \cdot \mathbf{n} = 0 \quad \text{on } \partial \Omega $$ | Homogeneous Neumann on all edges → insulated cavity; solution defined up to an arbitrary constant. | Implicitly assumes a closed, perfectly insulating box. Field lines cannot exit; potential is arbitrary up to a constant (singular system). |


