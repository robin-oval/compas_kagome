from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from compas.datastructures import Network

from compas.geometry import convex_hull

from compas.geometry import add_vectors
from compas.geometry import subtract_vectors
from compas.geometry import normalize_vector
from compas.geometry import scale_vector
from compas.geometry import distance_point_point
from compas.geometry import midpoint_line

from compas.datastructures import meshes_join_and_weld

from compas.datastructures import trimesh_subdivide_loop


__all__ = [
	'trimesh_skeleton'
	]


def trimesh_skeleton(cls, lines, radius=1):

	network = Network.from_lines(lines)

	tube_extremities = {}

	nodes = []
	for vkey in network.nodes():
		if len(network.adjacency[vkey]) > 1:
			
			points = [network.edge_point(vkey, nbr, t = float(radius) / network.edge_length(vkey, nbr)) for nbr in network.adjacency[vkey]]
			faces = convex_hull(points)
			mesh = cls.from_vertices_and_faces(points, faces)
			
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
				#print network.vertex_neighbors(vkey), network.vertex_neighbors(vkey)[vkey_2]
				#tube_extremities[(vkey, network.adjacency[vkey][vkey_2])] = tops
					
			mesh = meshes_join_and_weld(meshes)
			
			#dense_mesh = trimesh_subdivide_loop(mesh, k = 3)
			
			nodes.append(mesh)

	return nodes[0]

	meshes_2 = []
	for u, v in network.edges():
		if len(network.vertex_neighbors(u)) > 1 and len(network.vertex_neighbors(v)) > 1:
			#print len(tube_extremities[(u, v)])
			#print len(tube_extremities[(v, u)])
			if len(tube_extremities[(u, v)]) == len(tube_extremities[(v, u)]):
				n = len(tube_extremities[(u, v)])
				l = network.edge_length(u, v) - 2 * radius
				m = math.floor(l / radius) + 1
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
				array = [pt_uv]
				for i in range(int(m)):
					polygon = []
					for j in range(int(n)):
						u = pt_uv[j]
						v = pt_vu[j]
						polygon.append(add_vectors(scale_vector(u, (float(m) - 1 - float(i))/float(m - 1)), scale_vector(v, float(i)/float(m - 1))))
					array.append(polygon)
				array.append(pt_vu)
				#print len(array), len(array[0]), len(array[1]), len(array[2]), len(array[3])
				for i in range(int(n)):
					for j in range(int(m)):
						vertices = [array[i - 1][j - 1], array[i - 1][j], array[i][j]]
						faces = [[0, 1, 2]]
						meshes_2.append(Mesh.from_vertices_and_faces(vertices, faces))

	vertices, faces = join_and_weld_meshes(meshes_2)

	#meshes_2 = rs.AddMesh(vertices, faces)

	meshes = []
	for node in nodes:
		vertices, faces = node.to_vertices_and_faces()
		meshes.append(rs.AddMesh(vertices, faces))


# ==============================================================================
# Main
# ==============================================================================

if __name__ == '__main__':

	from compas.datastructures import Mesh

	a = [0, 0, 0]
	b = [0, 0, 1]
	c = [1, 0, -1]
	d = [-1, 0, -1]
	e = [0, 1, -1]
	f = [0, -1, -1]
	lines = [[a, b], [a, c], [a, d], [a, e], [a, f]]
	
	mesh = trimesh_skeleton(Mesh, lines, radius=.1)
