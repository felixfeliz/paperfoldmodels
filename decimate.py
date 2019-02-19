import openmesh as om
import numpy as np

#mesh = om.TriMesh()

#om.read_trimesh(mesh, "models/icosahedron.obj")

mesh = om.read_trimesh('models/bunny.stl')

mesh.triangulate()

om.write_mesh('original.off',mesh)

number = 0
for face in mesh.faces():
    number = number + 1

print ("No. triangles: " + str(number))

decimater = om.TriMeshDecimater(mesh)
moduleHandleQuad = om.TriMeshModQuadricHandle()
#moduleHandle = om.TriMeshModEdgeLengthHandle()

#decimater.add(moduleHandle)
decimater.add(moduleHandleQuad)


#om.TriMeshModEdgeLength.set_edge_length()
#decimater.module(moduleHandle).set_edge_length(3.5)
#decimater.module(moduleHandle).set_binary(True)
decimater.module(moduleHandleQuad).set_binary(False)

#print(decimater.module(moduleHandle).is_binary())
#decimater.module(moduleHandle).set_aspect_ratio(0.3)
#decimater.module(moduleHandle).set_error_tolerance_factor(0.1)



decimater.initialize()
print(decimater.is_initialized())
a = decimater.decimate_to_faces(0,120)
print(a)
mesh.garbage_collection()
om.write_mesh('reduced.obj',mesh)

print ("No. triangles: " + str(mesh.n_faces()))