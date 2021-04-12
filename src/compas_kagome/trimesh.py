from __future__ import print_function
from __future__ import absolute_import
from __future__ import division


__all__ = ['TriMesh']


class TriMesh(Mesh):

	def __init__(self):
		super(TriMesh, self).__init__()

	### singularities ###

	def singularities(self):
		return [vkey for vkey in self.vertices() if (self.is_vertex_on_boundary(vkey) and self.vertex_degree(vkey) != 4) or (not self.is_vertex_on_boundary(vkey) and self.vertex_degree(vkey) != 6)]

	def singularity_points(self):
		return [self.vertex_coordinates(vkey) for vkey in self.singularities()]


# ==============================================================================
# Main
# ==============================================================================

if __name__ == '__main__':
	pass
