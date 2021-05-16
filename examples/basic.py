from math import pi, cos, sin
from compas.utilities import pairwise
from compas.datastructures import Mesh
from compas_kagome.kagome import Kagome

from compas_view2.app import App
from compas.geometry import Point

# Kagome pattern for a simple patch with one singularity...
n = 5
vertices = [[cos(i / n * 2 * pi), sin(i / n * 2 * pi), 0.0] for i in range(n)] + [[0.0, 0.0, 0.0]]
faces = [[i, j, n] for i, j in pairwise(list(range(n)) + [0])]
coarse_tri_mesh = Mesh.from_vertices_and_faces(vertices, faces) # hard-coded coarse triangulated mesh

# ...or Kagome pattern from a polyhedron
# coarse_tri_mesh = Mesh.from_polyhedron(20) # 4 or 20 for tetrahedron (4 faces) and icosahedron (20 faces)

kagome_mesh = Kagome.from_mesh(coarse_tri_mesh, k = 3) # play with density value k

viewer = App()
for fkey in kagome_mesh.singularities():
	viewer.add(Point(*kagome_mesh.face_centroid(fkey)))
mesh = Mesh.from_vertices_and_faces(*kagome_mesh.to_vertices_and_faces())
viewer.add(mesh, show_edges=True)
viewer.run()