from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from compas.datastructures import Mesh

from compas.geometry import subtract_vectors
from compas.geometry import normalize_vector
from compas.geometry import scale_vector
from compas.geometry import cross_vectors

from compas.datastructures import trimesh_subdivide_loop

from compas.datastructures import mesh_conway_ambo

from compas.utilities import pairwise
from compas.utilities import window

__all__ = ['Kagome']


class Kagome(Mesh):

	def __init__(self):
		super(Kagome, self).__init__()
		self.polyedge_data = None

	### from ###

	@classmethod
	def from_mesh(cls, coarse_mesh, k = 1, fixed_boundary = True):
		if k > 0:
			fixed = coarse_mesh.vertices_on_boundary() if fixed_boundary else None
			dense_mesh = Mesh.from_vertices_and_faces(*trimesh_subdivide_loop(coarse_mesh, k, fixed).to_vertices_and_faces())
		else:
			dense_mesh = coarse_mesh
		vertices, faces = mesh_conway_ambo(dense_mesh).to_vertices_and_faces()
		return cls.from_vertices_and_faces(vertices, faces)

	### singularities ###

	def store_polyedge_data(self):
		self.polyedge_data = self.polyedges()

	def singularities(self):

		singular_faces = []
		for fkey in self.faces():
			singular_hex = all([len(self.face_vertices(nbr)) == 3 for nbr in self.face_neighbors(fkey)]) and len(self.face_vertices(fkey)) != 6
			singular_tri = all([len(self.face_vertices(nbr)) == 6 for nbr in self.face_neighbors(fkey)]) and len(self.face_vertices(fkey)) != 3
			if singular_hex or singular_tri:
				singular_faces.append([self.vertex_coordinates(vkey) for vkey in self.face_vertices(fkey) + self.face_vertices(fkey)[: 1]])
		return singular_faces

	def negative_singularities(self):

		return [fkey for fkey in self.faces() if all([len(self.face_vertices(nbr)) == 3 for nbr in self.face_neighbors(fkey)]) and len(self.face_vertices(fkey)) > 6]

	def negative_polygons(self):

		return [[self.vertex_coordinates(vkey) for vkey in self.face_vertices(fkey)] for fkey in self.negative_singularities()]

	### polyedges ###

	def vertex_opposite_vertex(self, u, v):

		if self.is_edge_on_boundary(u, v):
			
			if self.vertex_degree(v) == 2:
				return None
			
			else:
				return [nbr for nbr in self.vertex_neighbors(v, ordered = True) if nbr != u and self.is_edge_on_boundary(v, nbr)][0]
		
		elif self.is_vertex_on_boundary(v):
			return None

		else:
			nbrs = self.vertex_neighbors(v, ordered = True)
			return nbrs[nbrs.index(u) - 2]

	def polyedge(self, u0, v0):

		polyedge = [u0, v0]

		while len(polyedge) <= self.number_of_vertices():

			# end if closed loop
			if polyedge[0] == polyedge[-1]:
				break

			# get next vertex accros four-valent vertex
			w = self.vertex_opposite_vertex(*polyedge[-2 :])

			# flip if end of first extremity
			if w is None:
				polyedge = list(reversed(polyedge))
				# stop if end of second extremity
				w = self.vertex_opposite_vertex(*polyedge[-2 :])
				if w is None:
					break

			# add next vertex
			polyedge.append(w)

		return polyedge

	def polyedges(self):

		polyedges = []

		edge_visited = {(u, v): False for u, v in self.edges()}

		for edge in self.edges():
			if edge_visited[edge]:
				continue
			u0, v0 = edge
			# collect new polyedge
			polyedges.append(self.polyedge(u0, v0))

			# remove collected edges
			for u, v in pairwise(polyedges[-1]):
				#if (u, v) in edge_visited:
				edge_visited[(u, v)] = True
				#elif (v, u) in edge_visited:
				edge_visited[(v, u)] = True

		return polyedges

	def polyline(self, u, v):

		return [self.vertex_coordinates(vkey) for vkey in self.polyedge(u0, v0)]

	def polylines(self):

		return [[self.vertex_coordinates(vkey) for vkey in polyedge] for polyedge in self.polyedge_data]

	def polyline_frames(self):
		polylines_frames = []
		for polyedge in self.polyedge_data:
			polyline_frames = []
			for i, u in enumerate(polyedge):
				#if end
				if i == len(polyedge) - 1:
					# if closed
					if polyedge[0] == polyedge[-1]:
						v = polyedge[1]
					else:
						v = polyedge[i - 1]
				else:
					v = polyedge[i + 1]
				x = self.vertex_normal(u)
				y = normalize_vector(subtract_vectors(self.vertex_coordinates(v), self.vertex_coordinates(u)))
				if i == len(polyedge) - 1 and polyedge[0] != polyedge[-1]:
					y = scale_vector(y, -1)
				z = cross_vectors(x, y)
				polyline_frames.append([x, y, z])
			polylines_frames.append(polyline_frames)
		return polylines_frames

	### weave ###

	def polyedge_weaving(self):

		edge_to_polyedge_index = {}
		for i, polyedge in enumerate(self.polyedge_data):
			for u, v in pairwise(polyedge):
				edge_to_polyedge_index[(u, v)] = i
				edge_to_polyedge_index[(v, u)] = i

		vertex_to_polyege_offset = {vkey: {} for vkey in self.vertices()}
		for fkey in self.faces():
			if len(self.face_vertices(fkey)) == 3:
				for u, v, w in window(self.face_vertices(fkey) + self.face_vertices(fkey)[:2], n = 3):
					vertex_to_polyege_offset[v].update({edge_to_polyedge_index[(u, v)]: +1, edge_to_polyedge_index[(v, w)]: -1})
			else:
				for u, v, w in window(self.face_vertices(fkey) + self.face_vertices(fkey)[:2], n = 3):
					vertex_to_polyege_offset[v].update({edge_to_polyedge_index[(u, v)]: -1, edge_to_polyedge_index[(v, w)]: +1})

		polyedge_weave = []
		for i, polyedge in enumerate(self.polyedge_data):
			polyedge_weave.append([vertex_to_polyege_offset[vkey][i] for vkey in polyedge])

		return polyedge_weave


# ==============================================================================
# Main
# ==============================================================================

if __name__ == '__main__':

	import compas
	from compas_kagome.colouring import kagome_polyline_colouring
	from compas_plotters import MeshPlotter

	vertices = [
		[0., 0., 0.],
		[1., 0., 0.],
		[1., 1., 0.],
		[0., 1., 0.],
		[0.5, 0.5, 0.],	
	]

	faces = [
		[0, 1, 4],
		[1, 2, 4],
		[2, 3, 4],
		[3, 0, 4],
	]

	lines = [
		([0., 0., 0.],[1., 0., -1.]),
		([0., 0., 0.],[-1., 0., -1.]),
		([0., 0., 0.],[0., 1., -1.]),
		([0., 0., 0.],[0., -1., -1.]),
		([0., 0., 0.],[0., 0., 1.]),
		]

	#mesh = Mesh.from_skeleton(lines)
	mesh = Mesh.from_vertices_and_faces(vertices, faces)
	kagome = Kagome.from_mesh(mesh, k = 2, fixed_boundary = True)

	kagome.store_polyedge_data()

	plotter = MeshPlotter(kagome)
	plotter.draw_vertices(radius=.005)
	plotter.draw_edges()
	plotter.draw_faces()
	plotter.show()

	# print(kagome.negative_singularities())
	# print(kagome.singularities())
	# print(kagome.polyline_frames())
	# kagome.polyedges()
	kagome_polyline_colouring(kagome)
	kagome.polyedge_weaving()
