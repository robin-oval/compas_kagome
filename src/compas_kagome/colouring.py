from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from compas.datastructures import Network

from compas.geometry import centroid_points

from compas.topology import vertex_coloring

from compas.utilities import pairwise


__all__ = [
	'kagome_polyedge_colouring',
	'kagome_polyline_colouring'
	]


def kagome_polyedge_colouring(kagome):

	polyedges = kagome.polyedge_data

	edge_to_polyedge_index = {vkey: {} for vkey in kagome.vertices()}
	for i, polyedge in enumerate(polyedges):
		for u, v in pairwise(polyedge):
			edge_to_polyedge_index[u][v] = i
			edge_to_polyedge_index[v][u] = i

	vertices = [centroid_points([kagome.vertex_coordinates(vkey) for vkey in polyedge]) for polyedge in polyedges]

	edges = []
	for idx, polyedge in enumerate(polyedges):
		for vkey in polyedge:
			for vkey_2 in kagome.vertex_neighbors(vkey):
				idx_2 = edge_to_polyedge_index[vkey][vkey_2]
				if idx_2 != idx and idx < idx_2 and (idx, idx_2) not in edges:
					edges.append((idx, idx_2))

	polyedge_network = Network.from_nodes_and_edges(vertices, edges)

	key_to_colour = vertex_coloring(polyedge_network.adjacency)

	return [key_to_colour[key] for key in sorted(key_to_colour.keys())]

def kagome_polyline_colouring(kagome):

	return {tuple([tuple(kagome.vertex_coordinates(vkey)) for vkey in polyedge]): colour for polyedge, colour in kagome_polyedge_colouring(kagome).items()}


# ==============================================================================
# Main
# ==============================================================================

if __name__ == '__main__':
	pass
