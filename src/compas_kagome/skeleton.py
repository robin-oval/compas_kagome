from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from math import floor

from compas.datastructures import Network
from compas.geometry import Polyline

from compas.geometry import convex_hull, discrete_coons_patch

from compas.geometry import add_vectors
from compas.geometry import subtract_vectors
from compas.geometry import normalize_vector
from compas.geometry import scale_vector
from compas.geometry import distance_point_point
from compas.geometry import midpoint_line

from compas.datastructures import mesh_unweld_vertices
from compas.datastructures import meshes_join_and_weld

from compas.datastructures import trimesh_subdivide_loop

from compas.utilities import geometric_key


__all__ = [
	'trimesh_skeleton'
	]


def insert_triface_on_boundary(mesh, bdry_vkey, nbr_vkey):

	# topology
	fkey = mesh.halfedge[nbr_vkey][bdry_vkey]
	faces = mesh.vertex_faces(bdry_vkey, ordered=True)
	i = faces.index(fkey)
	left_faces = faces[i:]
	x, y, z = mesh.vertex_coordinates(bdry_vkey)
	attr = {'x': x, 'y': y, 'z': z}
	new_vkey = mesh.add_vertex(attr_dict=attr)
	for face in left_faces:
		new_vertices = [vkey if vkey != bdry_vkey else new_vkey for vkey in mesh.face_vertices(face)]
		mesh.delete_face(face)
		mesh.add_face(new_vertices, fkey=face)

	new_vkeys = [bdry_vkey, new_vkey]

	# geometry
	for vkey in new_vkeys:
		for nbr in mesh.vertex_neighbors(vkey):
			if mesh.is_edge_on_boundary(vkey, nbr) and nbr != nbr_vkey:
				x, y, z = mesh.edge_point(vkey, nbr, t=.33)
				mesh.vertex[vkey]['x'] = x
				mesh.vertex[vkey]['y'] = y
				mesh.vertex[vkey]['z'] = z
				break

	# add new face
	mesh.add_face([nbr_vkey] + new_vkeys)

	return new_vkey


def trimesh_skeleton(cls, lines, radius=1):

	network = Network.from_lines(lines)

	tube_extremities = {}

	nodes = []
	for vkey in network.nodes():
		if len(network.adjacency[vkey]) > 1:
			
			points = {nbr: network.edge_point(vkey, nbr, t = float(radius) / network.edge_length(vkey, nbr)) for nbr in network.adjacency[vkey]}
			idx_to_key = {i: key for i, key in enumerate(network.adjacency[vkey])}
			faces = convex_hull(list(points.values()))
			faces = [[idx_to_key[idx] for idx in face] for face in faces]
			mesh = cls.from_vertices_and_faces(points, faces)
			nodes.append(mesh)

			meshes = []
			
			for fkey in mesh.faces():
				vertices = [mesh.edge_midpoint(u, v) for u, v in mesh.face_halfedges(fkey)]
				faces = [[0,1,2]]
				meshes.append(cls.from_vertices_and_faces(vertices, faces))

			for vkey_2 in mesh.vertices():
				tops = []
				bottoms = []
				n = normalize_vector(subtract_vectors(mesh.vertex_coordinates(vkey_2), network.node_coordinates(vkey)))
				for i in range(len(mesh.vertex_neighbors(vkey_2))):
					pt_0 = mesh.edge_midpoint(vkey_2, mesh.vertex_neighbors(vkey_2, ordered = True)[i - 1])
					bottoms.append(pt_0)
					pt_1 = mesh.edge_midpoint(vkey_2, mesh.vertex_neighbors(vkey_2, ordered = True)[i])
					pt_2 = midpoint_line([pt_0, pt_1])
					pt_2 = add_vectors(scale_vector(n, distance_point_point(pt_0, pt_1)), pt_2)
					tops.append(pt_2)
					vertices = [pt_0, pt_2, pt_1]
					faces = [[0,1,2]]
					meshes.append(cls.from_vertices_and_faces(vertices, faces))
				for i in range(len(tops)):
					vertices = [tops[i - 1], tops[i], bottoms[i]]
					faces = [[0,1,2]]
					meshes.append(cls.from_vertices_and_faces(vertices, faces))

				tube_extremities[(vkey, vkey_2)] = tops
					
			mesh = meshes_join_and_weld(meshes)

			nodes.append(mesh)

	all_nodes = meshes_join_and_weld(nodes)

	
	# leaf node ring
	for u in network.nodes():
		if len(network.adjacency[u]) == 1:
			v = next(iter(network.adjacency[u].keys()))
			ring_v = tube_extremities[(v, u)]
			ring_u = [add_vectors(pt, network.edge_vector(v, u)) for pt in ring_v[::-1]]
			tube_extremities[(u, v)] = ring_u

	beams = []
	for u, v in network.edges():

		geom_key_map = {geometric_key(all_nodes.vertex_coordinates(vkey)): vkey for vkey in all_nodes.vertices()}
		if len(tube_extremities[(u, v)]) != len(tube_extremities[(v, u)]):
			if len(tube_extremities[(u, v)]) < len(tube_extremities[(v, u)]):
				a, b = u, v
			else:
				a, b = v, u
			n = len(tube_extremities[(b, a)]) - len(tube_extremities[(a, b)])
			for i in range(n):
				bdry_vkey = geom_key_map[geometric_key(tube_extremities[(a, b)][0])]
				nbr_vkey = None
				for nbr in all_nodes.vertex_neighbors(bdry_vkey):
					if not all_nodes.is_edge_on_boundary(bdry_vkey, nbr):
						nbr_vkey = nbr
				k = tube_extremities[(a, b)].index(all_nodes.vertex_coordinates(bdry_vkey))
				new_vkey = insert_triface_on_boundary(all_nodes, bdry_vkey, nbr_vkey)
				tube_extremities[(a, b)].insert(k + 1 - len(tube_extremities[(a, b)]), all_nodes.vertex_coordinates(new_vkey))

		if len(tube_extremities[(u, v)]) == len(tube_extremities[(v, u)]):
			n = len(tube_extremities[(u, v)])
			l = network.edge_length(u, v) - 2 * radius
			m = (floor(l / radius) + 1) * 2
			pt_uv = tube_extremities[(u, v)]
			pt_vu = list(reversed(tube_extremities[(v, u)]))
			dmin = -1
			imin = None
			for i in range(n):
				distance = sum([distance_point_point(pt_uv[j], pt_vu[i + j - len(pt_vu)]) for j in range(n)])
				if dmin < 0 or distance < dmin:
					dmin = distance
					imin = i
			pt_vu = [pt_vu[imin + j - len(pt_vu)] for j in range(n)]
			
			# ab = pt_uv# + pt_uv[0:]
			# dc = pt_vu# + pt_vu[0:]
			# line = Polyline([ab[0], dc[0]])
			# ad = [line.point(i / (m - 1)) for i in range(m)]
			# bc = ad
			# vertices, faces = discrete_coons_patch(ab, bc, dc, ad)
			# tri_faces = []
			# for (a, b, c, d) in faces:
			# 	tri_faces += [[a, b, c], [a, c, d]] # reverse?
			# beams.append(Mesh.from_vertices_and_faces(vertices, tri_faces))
			
			array = [pt_uv]
			for i in range(int(m)):
				polygon = []
				for j in range(int(n)):
					u = pt_uv[j]
					v = pt_vu[j]
					polygon.append(add_vectors(scale_vector(u, (float(m) - 1 - float(i))/float(m - 1)), scale_vector(v, float(i)/float(m - 1))))
				array.append(polygon)
			array.append(pt_vu)
			for i in range(1, int(m + 2)):
				for j in range(int(n)):
					vertices = [array[i - 1][j - 1], array[i - 1][j], array[i][j - 1], array[i][j]] # create staggered pattern?
					faces = [[3, 1, 0], [2, 3, 0]]
					beams.append(cls.from_vertices_and_faces(vertices, faces))


	# return all_nodes
	# #return meshes_join_and_weld(beams)
	return meshes_join_and_weld([all_nodes] + beams)


# ==============================================================================
# Main
# ==============================================================================

if __name__ == '__main__':

	from compas_kagome.kagome import Kagome
	from compas.datastructures import Mesh
	from compas.datastructures import mesh_weld, mesh_conway_ambo, mesh_smooth_area
	from compas.geometry import Point
	from compas_view2.app import App

	# a = [0, 0, 0]
	# b = [0, 0, 1]
	# c = [1, 0, -1]
	# d = [-1, 0, -1]
	# e = [0, 1, -1]
	# f = [0, -1, -1]
	# g = [1, 0, 2]
	# h = [-1, 0, 2]
	# i = [0, 1, 2]
	# j = [0, -1, 2]
	# lines = [[a, b], [a, c], [a, d], [a, e], [a, f], [b, g], [b, h], [b, i], [b, j]]

	# a = [0, 0, 0]
	# b = [0, 0, 1]
	# c = [1, 0, -1]
	# d = [-1, 0, -1]
	# e = [0, 1, -1]
	# f = [0, -1, -1]
	# g = [0, 1, 2]
	# h = [1, -1, 2]
	# i = [-1, -1, 2]
	# lines = [[a, b], [a, c], [a, d], [a, e], [a, f], [b, g], [b, h], [b, i]]

	lines = [[[0.0, 4.0, 0.0], [0.0, 6.0, 2.0]], [[0.0, 4.0, 0.0], [0.0, 6.0, -2.0]], [[-3.8042260651806141, 1.2360679774997898, 0.0], [-5.706339097770921, 1.8541019662496847, 2.0]], [[-3.8042260651806141, 1.2360679774997898, 0.0], [-5.706339097770921, 1.8541019662496847, -2.0]], [[-2.351141009169893, -3.2360679774997894, 0.0], [-3.5267115137548393, -4.8541019662496838, 2.0]], [[-2.351141009169893, -3.2360679774997894, 0.0], [-3.5267115137548393, -4.8541019662496838, -2.0]], [[2.3511410091698921, -3.2360679774997898, 0.0], [3.5267115137548384, -4.8541019662496847, 2.0]], [[2.3511410091698921, -3.2360679774997898, 0.0], [3.5267115137548384, -4.8541019662496847, -2.0]], [[3.8042260651806146, 1.2360679774997889, 0.0], [5.7063390977709219, 1.8541019662496834, 2.0]], [[3.8042260651806146, 1.2360679774997889, 0.0], [5.7063390977709219, 1.8541019662496834, -2.0]], [[0.0, 4.0, 0.0], [-3.8042260651806141, 1.2360679774997898, 0.0]], [[-3.8042260651806141, 1.2360679774997898, 0.0], [-2.351141009169893, -3.2360679774997894, 0.0]], [[-2.351141009169893, -3.2360679774997894, 0.0], [2.3511410091698921, -3.2360679774997898, 0.0]], [[2.3511410091698921, -3.2360679774997898, 0.0], [3.8042260651806146, 1.2360679774997889, 0.0]], [[3.8042260651806146, 1.2360679774997889, 0.0], [0.0, 4.0, 0.0]]]
	#lines = [[[0.0, 4.0, 0.0], [0.0, 6.0, 2.0]], [[0.0, 4.0, 0.0], [0.0, 6.0, -2.0]], [[-3.8042260651806141, 1.2360679774997898, 0.0], [-5.706339097770921, 1.8541019662496847, 2.0]], [[-3.8042260651806141, 1.2360679774997898, 0.0], [-5.706339097770921, 1.8541019662496847, -2.0]], [[-2.351141009169893, -3.2360679774997894, 0.0], [-3.5267115137548393, -4.8541019662496838, 2.0]], [[-2.351141009169893, -3.2360679774997894, 0.0], [-3.5267115137548393, -4.8541019662496838, -2.0]], [[2.3511410091698921, -3.2360679774997898, 0.0], [3.5267115137548384, -4.8541019662496847, 2.0]], [[2.3511410091698921, -3.2360679774997898, 0.0], [3.5267115137548384, -4.8541019662496847, -2.0]], [[3.8042260651806146, 1.2360679774997889, 0.0], [5.7063390977709219, 1.8541019662496834, 2.0]], [[3.8042260651806146, 1.2360679774997889, 0.0], [5.7063390977709219, 1.8541019662496834, -2.0]], [[0.0, 4.0, 0.0], [-3.8042260651806141, 1.2360679774997898, 0.0]], [[-3.8042260651806141, 1.2360679774997898, 0.0], [-2.351141009169893, -3.2360679774997894, 0.0]], [[-2.351141009169893, -3.2360679774997894, 0.0], [2.3511410091698921, -3.2360679774997898, 0.0]], [[2.3511410091698921, -3.2360679774997898, 0.0], [3.8042260651806146, 1.2360679774997889, 0.0]], [[3.8042260651806146, 1.2360679774997889, 0.0], [0.0, 4.0, 0.0]], [[-2.351141009169893, -3.2360679774997894, 0.0], [0.0, 0.0, 3.0]], [[-2.351141009169893, -3.2360679774997894, 0.0], [0.0, 0.0, -3.0]], [[2.3511410091698921, -3.2360679774997902, 0.0], [0.0, 0.0, 3.0]], [[3.8042260651806146, 1.2360679774997887, 0.0], [0.0, 0.0, 3.0]], [[8.8817841970012523e-16, 3.9999999999999996, 0.0], [0.0, 0.0, 3.0]], [[-3.8042260651806137, 1.2360679774997911, 0.0], [0.0, 0.0, 3.0]], [[2.3511410091698921, -3.2360679774997902, 0.0], [0.0, 0.0, -3.0]], [[3.8042260651806146, 1.2360679774997887, 0.0], [0.0, 0.0, -3.0]], [[8.8817841970012523e-16, 3.9999999999999996, 0.0], [0.0, 0.0, -3.0]], [[-3.8042260651806137, 1.2360679774997911, 0.0], [0.0, 0.0, -3.0]], [[0.0, 0.0, -3.0], [0.0, 0.0, -5.0]], [[0.0, 0.0, 5.0], [0.0, 0.0, 3.0]]]

	mesh = trimesh_skeleton(Kagome, lines, radius=2)
	mesh = trimesh_subdivide_loop(mesh, k=3)
	mesh = mesh_conway_ambo(mesh)
	fixed = [vkey for bdry in mesh.vertices_on_boundaries() for vkey in bdry]
	mesh_smooth_area(mesh, fixed=fixed, kmax=50, damping=0.5)
	viewer = App()
	#for fkey in mesh.singularities():
	#	viewer.add(Point(*mesh.face_centroid(fkey)))
	vertices, faces = mesh.to_vertices_and_faces()
	mesh = Mesh.from_vertices_and_faces(vertices, faces)
	viewer.add(mesh, show_edges=True)
	viewer.run()

