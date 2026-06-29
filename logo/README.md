## Generate the logo

Generate the mesh
```
python xemfem_logo_mesh.py --out xemfem_logo.msh --geo xemfem_logo.geo_unrolled --show
```

Render into png

```
python render_logo.py xemfem_logo.msh   --out xemfem_logo_mesh_render.png   --line-source mesh   --fill-alpha 0.03   --line-width 0.35   --dpi 1200   --width-in 10
```
