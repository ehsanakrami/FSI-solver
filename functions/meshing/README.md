*) Meshing of a specific test-case geometry.

A vertical solid rectangular plate is immersed in a box of fluid, joined at the base level. 
The plate's height can be smaller, equal or larger than the fluid's.

Meshing is done using block elements. 
For pressure DOFs, linear 8-node blocks are used.
For velocity DOFs, quadratic 20-node blocks are used.

Gauss quadrature schemes are available with 8 or 27 quadrature points.

Connectivity tables are available, for both the 8-node and the 20-node meshes, for the solid and fluid domains.
2D connectivity tables (4-node and 8-node rectangular surface elements) are available for each face of the solid and fluid domains,
and for the whole fluid-structure interface.

*) Sample figures

The fluid domain's outer dimensions are 1m in length, 0.5m in width and height.
The solid domain is roughly 0.1m, 0.3m and 0.67m in those directions, respectively.

sample_mesh_8.png  : 540 elements for the solid domain, 6120 elements for the fluid domain, 7880 nodes.
sample_mesh_20.png : 540 elements for the solid domain, 6120 elements for the fluid domain, 30234 nodes. Extra nodes per element are represented with extra lines forming subtriangles on each face.
