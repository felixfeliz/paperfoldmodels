import numpy as np
import openmesh as om

from time import sleep

def findTree(vertex, forest):
    """Find the index of the tree in forest that contains vertex"""
    for i in range(len(forest)):  # type: int
        for element in forest[i]:
            if element == vertex:
                return i
    return None


def mergeTrees(vertex1, vertex2, forest):
    index1 = findTree(vertex1, forest)
    index2 = findTree(vertex2, forest)
    for element in forest[index2]:
        forest[index1].append(element)
    forest.pop(index2)


def unrollTree(face, lastHalfEdge, unrolledLastHalfEdge, oppositeUnrolledVertex, halfEdgeTree, mesh, unfoldedMesh,
               isfoldingEdge):
    #print("Processing face " + str(face.idx()))

    # First unroll the current face
    # Get the two known unrolled Vertices
    unrolledFromVertex = unfoldedMesh.to_vertex_handle(unrolledLastHalfEdge)
    unrolledToVertex = unfoldedMesh.from_vertex_handle(unrolledLastHalfEdge)

    # Get the edge lengths
    lastHalfEdgeInFace = mesh.opposite_halfedge_handle(lastHalfEdge)
    secondHalfEdgeInFace = mesh.next_halfedge_handle(lastHalfEdgeInFace)
    thirdHalfEdgeInFace = mesh.next_halfedge_handle(secondHalfEdgeInFace)
    edgelengths = [mesh.calc_edge_length(lastHalfEdgeInFace), mesh.calc_edge_length(secondHalfEdgeInFace),
                   mesh.calc_edge_length(thirdHalfEdgeInFace)]

    # Find the third unrolled Point
    # print([unfoldedMesh.point(unrolledFromVertex),unfoldedMesh.point(unrolledToVertex), edgelengths ])
    [newUnrolledVertex0, newUnrolledVertex1] = getThirdPoint(unfoldedMesh.point(unrolledFromVertex),
                                                             unfoldedMesh.point(unrolledToVertex), edgelengths[0],
                                                             edgelengths[1], edgelengths[2])
    # Check which one is on the opposite side of the edge from  oppositeUnrolledVertex
    lastUnrolledEdgeVector = unfoldedMesh.point(unrolledToVertex) - unfoldedMesh.point(unrolledFromVertex)
    newUnrolledVertexVector0 = newUnrolledVertex0 - unfoldedMesh.point(unrolledFromVertex)

    newUnrolledVertexVector1 = newUnrolledVertex1 - unfoldedMesh.point(unrolledFromVertex)
    oppositeVector = unfoldedMesh.point(oppositeUnrolledVertex) - unfoldedMesh.point(unrolledFromVertex)

    if np.dot(np.cross(lastUnrolledEdgeVector, newUnrolledVertexVector0),
              np.cross(lastUnrolledEdgeVector, oppositeVector)) <= 0:
        newUnrolledVertex = unfoldedMesh.add_vertex(newUnrolledVertex0)
    else:
        newUnrolledVertex = unfoldedMesh.add_vertex(newUnrolledVertex1)

    # Make the face
    unfoldedMesh.add_face(unrolledFromVertex, unrolledToVertex, newUnrolledVertex)

    # Make the last edge a folding edge
    lastEdge = unfoldedMesh.edge_handle(unrolledLastHalfEdge)
    isFoldingEdge[lastEdge.idx()] = True

    #print("Unfolded face " + str(face.idx()) + ".")
    # om.write_mesh('unfolding' + str(face.idx()) + ".off", unfoldedMesh)

    # Now find the neighbours
    # Get the unrolled halfEdge
    secondUnrolledHalfEdge = unfoldedMesh.next_halfedge_handle(
        unfoldedMesh.opposite_halfedge_handle(unrolledLastHalfEdge))
    thirdUnrolledHalfEdge = unfoldedMesh.next_halfedge_handle(secondUnrolledHalfEdge)

    # Check the two other half edges
    if secondHalfEdgeInFace in halfEdgeTree:
        # Get the face
        neighbourFace = mesh.face_handle(mesh.opposite_halfedge_handle(secondHalfEdgeInFace))
        unrollTree(neighbourFace, secondHalfEdgeInFace, secondUnrolledHalfEdge, unrolledFromVertex, halfEdgeTree, mesh,
                   unfoldedMesh, isFoldingEdge)
    if thirdHalfEdgeInFace in halfEdgeTree:
        # Get the face
        neighbourFace = mesh.face_handle(mesh.opposite_halfedge_handle(thirdHalfEdgeInFace))
        unrollTree(neighbourFace, thirdHalfEdgeInFace, thirdUnrolledHalfEdge, unrolledToVertex, halfEdgeTree, mesh,
                   unfoldedMesh, isFoldingEdge)


def findLeafIndex(forest, spanningTree, mesh):
    for i in range(len(forest)):
        outerEdges = []

        for edge in spanningTree:
            dualVertex1 = mesh.face_handle(mesh.halfedge_handle(edge, 0))
            dualVertex2 = mesh.face_handle(mesh.halfedge_handle(edge, 1))
            if (dualVertex1 in forest[i]) != (dualVertex2 in forest[i]):
                outerEdges.append(edge)

        if len(outerEdges) == 1:
            # Find the parent
            dualVertex1 = mesh.face_handle(mesh.halfedge_handle(outerEdges[0], 0))
            dualVertex2 = mesh.face_handle(mesh.halfedge_handle(outerEdges[0], 1))
            if dualVertex1 in forest[i]:
                parent = dualVertex1
            else:
                parent = dualVertex2
            return [True, i, outerEdges[0], parent]
    return [False, -1, -1, -1]


def getThirdPoint(v0, v1, l01, l12, l20):
    v2rotx = (l01 ** 2 + l20 ** 2 - l12 ** 2) / (2 * l01)
    v2roty0 = np.sqrt((l01 + l20 + l12) * (l01 + l20 - l12) * (l01 - l20 + l12) * (-l01 + l20 + l12)) / (2 * l01)

    v2roty1 = - v2roty0

    theta = np.arctan2(v1[1] - v0[1], v1[0] - v0[0])
    # print(v2rotx, v2roty0)
    # print(theta)

    v2trans0 = np.array(
        [v2rotx * np.cos(theta) - v2roty0 * np.sin(theta), v2rotx * np.sin(theta) + v2roty0 * np.cos(theta), 0])
    v2trans1 = np.array(
        [v2rotx * np.cos(theta) - v2roty1 * np.sin(theta), v2rotx * np.sin(theta) + v2roty1 * np.cos(theta), 0])
    return [v2trans0 + v0, v2trans1 + v0]


def writeSVG(filename, mesh, isFoldingEdge):
    #Get the bounding box
    firstpoint = mesh.point(mesh.vertex_handle(0))
    xmin = firstpoint[0]
    xmax = firstpoint[0]
    ymin = firstpoint[1]
    ymax = firstpoint[1]
    for vertex in mesh.vertices():
        coordinates = mesh.point(vertex)
        if(coordinates[0] < xmin):
            xmin = coordinates[0]
        if(coordinates[0] > xmax):
            xmax = coordinates[0]
        if (coordinates[1] < ymin):
            ymin = coordinates[1]
        if (coordinates[1] > ymax):
            ymax = coordinates[1]
    boxSize = np.maximum(np.abs(xmax - xmin), np.abs(ymax - ymin))
    #xmin = np.minimum(xmin,ymin)
    #ymin = np.minimum(xmin,ymin)

    strokewidth = 0.002 * boxSize
    dashLength = 0.002*boxSize
    spaceLength = 0.003*boxSize

    #dashLength = 0.2
    #spaceLength = 0.1


    file = open(filename, 'w')
    file.write("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n")
    file.write("<svg width=\"29.6cm\" height=\"29.6cm\" viewBox = \""  + str(xmin) + " " + str(ymin) + " " + str(boxSize) + " " + str(boxSize) + "\" version = \"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n")

    #file.write("<rect x=\"" + str(xmin) + "\" y=\"" +str(ymin) + "\" width=\"100%\" height=\"100%\"/>\n")

    file.write("<g stroke = \"black\" stroke-width = \"" + str(strokewidth) + "\" >\n")
    for edge in mesh.edges():
        # Get the two points
        he = mesh.halfedge_handle(edge, 0)
        vertex0 = mesh.point(mesh.from_vertex_handle(he))
        vertex1 = mesh.point(mesh.to_vertex_handle(he))

        file.write("<line x1=\"" + str(vertex0[0]) + "\" y1=\"" + str(vertex0[1]) + "\" x2 = \"" + str(
            vertex1[0]) + "\" y2 = \"" + str(vertex1[1]) + "\"")

        if isFoldingEdge[edge.idx()]:
            file.write(" stroke-dasharray=\"" + str(dashLength) + " " + str(spaceLength) + "\"")
            #file.write(" stroke-dasharray=\"none\"")
            #file.write("style=\"stroke-dasharray:%.2f" %dashLength + ",%.2f" % spaceLength + ";stroke-dashoffset:0\"")
            #file.write(" style=\"fill:none;stroke:#000000;stroke-width:0.31999999;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:1;stroke-dasharray:1.28, 2.56;stroke-dashoffset:0;stroke-opacity:1\"")

        file.write(" />\n")

    file.write("</g>\n")
    file.write("</svg>")
    file.close()



#FILENAME = 'models/icosahedron.obj'
# FILENAME = 'original.off'
#FILENAME = 'reduced.off'
#FILENAME = 'models/polyhedron.obj'
FILENAME = 'models/kndC.obj'

# Import the mode
mesh = om.read_trimesh(FILENAME)

# Berechne die Anzahl der Flächen, Kanten und Ecken, sowie die Längen der längsten kürzesten Kante
numEdges = mesh.n_edges()
numVertices = mesh.n_vertices()
numFaces = mesh.n_faces()
minLength = 1000
maxLength = 0
for edge in mesh.edges():
    edgelength = mesh.calc_edge_length(edge)
    if edgelength < minLength:
        minLength = edgelength
    if edgelength > maxLength:
        maxLength = edgelength

print("Das Modell hat " + str(numEdges) + " Kanten, " + str(numFaces) + " Flächen und " + str(
    numVertices) + " Ecken.\nDie längste Kante hat Länge " + str(
    maxLength) + " und die kürzeste Kante hat Länge " + str(minLength) + ".")

# Compute the weights
# iterate over all edges
minPerimeterEdgeWeights = np.empty(numEdges)
flatAngleWeights = np.empty(numEdges)
dihedralAngleWeights = np.empty(numEdges)

#phi = -np.pi
#theta = 0.9*np.pi
#cutVector = np.array([np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta), np.cos(theta)])
#print(cutVector)
cutVector = np.array([1.0,1.0,1.0])

cutVector = cutVector / np.linalg.norm(cutVector)
for edge in mesh.edges():
    edgelength = mesh.calc_edge_length(edge)
    # Schneide kürzere Kanten um Umfang zu minimieren
    minPerimeterEdgeWeights[edge.idx()] = 1.0 - (edgelength - minLength) / (maxLength - minLength)

    #dihedralAngle = mesh.calc_dihedral_angle(edge)
    #if dihedralAngle >= 0:
    #    minPerimeterEdgeWeights[edge.idx()] = dihedralAngle
    #else:
    #    minPerimeterEdgeWeights[edge.idx()] = 10 - dihedralAngle
            #np.abs(mesh.calc_dihedral_angle(edge))

    edgeVector = mesh.calc_edge_vector(edge)
    flatAngleWeights[edge.idx()] = np.abs(np.dot(cutVector,edgeVector)) / np.linalg.norm(edgeVector)

weights = 0 *minPerimeterEdgeWeights + 1*flatAngleWeights


# Kruzkal's algorithm

# Sort the weights and find the order of the indices
sortedIndices = np.argsort(weights)
# Make ordered list of (dual) edges
sorteddualEdges = []
for i in range(numEdges):
    sorteddualEdges.append(mesh.edge_handle(sortedIndices[i]))

# Initialize the forest with all dual-vertices (i.e. faces)
spanningTree = []
forest = []
for face in mesh.faces():
    forest.append([face])
# Perform the step as long as there is more than one tree in the forest
i = 0
while len(forest) > 1:
    # Get the two faces adjacent to the ith edge
    dualVertex1 = mesh.face_handle(mesh.halfedge_handle(sorteddualEdges[i], 0))
    dualVertex2 = mesh.face_handle(mesh.halfedge_handle(sorteddualEdges[i], 1))

    # Check if both are in different trees
    index1 = findTree(dualVertex1, forest)
    index2 = findTree(dualVertex2, forest)
    if index1 != index2:
        # Merge the two trees
        for element in forest[index2]:
            forest[index1].append(element)
        forest.pop(index2)

        # Add the faces and vertices to the spanning tree

        spanningTree.append(sorteddualEdges[i])

    i = i + 1

sizeTree = len(spanningTree)
numUnfoldedEdges = 3 * numFaces - sizeTree
print("Der Baum hat " + str(sizeTree) + " Kanten.\nDas abgewickelte Netz wird " + str(
    numUnfoldedEdges) + " Kanten haben.")

# Make a tree of HalfEdges
halfEdgeTree = []
for edge in spanningTree:
    halfEdgeTree.append(mesh.halfedge_handle(edge, 0))
    halfEdgeTree.append(mesh.halfedge_handle(edge, 1))

# Unfolding
unfoldedMesh = om.TriMesh()
# Array to mark the folding edges
isFoldingEdge = np.zeros(numUnfoldedEdges, dtype=bool)
#Face onnection array
connections = np.empty(numFaces)

# Get the points of the triangle
startingTriangle = forest[0][0]

# Get all halfedges
firstHalfEdge = mesh.halfedge_handle(startingTriangle)
secondHalfEdge = mesh.next_halfedge_handle(firstHalfEdge)
thirdHalfEdge = mesh.next_halfedge_handle(secondHalfEdge)

edgelengths = [mesh.calc_edge_length(firstHalfEdge), mesh.calc_edge_length(secondHalfEdge),
               mesh.calc_edge_length(thirdHalfEdge)]

# The orientation of the first triangle is arbitrary
firstUnrolled = np.array([0, 0, 0])
secondUnrolled = np.array([edgelengths[0], 0, 0])

# Compute third point from lengths
[thirdUnrolled0, thirdUnrolled1] = getThirdPoint(firstUnrolled, secondUnrolled, edgelengths[0], edgelengths[1],
                                                 edgelengths[2])

if thirdUnrolled0[1] > 0:
    thirdUnrolled = thirdUnrolled0
else:
    thirdUnrolled = thirdUnrolled1

# Make triangle
# Add vertices
firstUnrolledVertex = unfoldedMesh.add_vertex(firstUnrolled)
secondUnrolledVertex = unfoldedMesh.add_vertex(secondUnrolled)
thirdUnrolledVertex = unfoldedMesh.add_vertex(thirdUnrolled)
# Create the face
f = unfoldedMesh.add_face(firstUnrolledVertex, secondUnrolledVertex, thirdUnrolledVertex)
# om.write_mesh('unfolding' + str(startingTriangle.idx()) + ".off", unfoldedMesh)
connections[startingTriangle.idx()] = f.idx()

# Now check the neighbours
# Find the unrolled half-edges inside the face
firstUnrolledHalfEdge = unfoldedMesh.opposite_halfedge_handle(unfoldedMesh.halfedge_handle(firstUnrolledVertex))
secondUnrolledHalfEdge = unfoldedMesh.opposite_halfedge_handle(unfoldedMesh.next_halfedge_handle(firstUnrolledHalfEdge))
thirdUnrolledHalfEdge = unfoldedMesh.opposite_halfedge_handle(unfoldedMesh.next_halfedge_handle(secondUnrolledHalfEdge))

if firstHalfEdge in halfEdgeTree:
    # Get the neighbouring face
    print("first")
    neighbourFace = mesh.face_handle(mesh.opposite_halfedge_handle(firstHalfEdge))
    unrollTree(neighbourFace, firstHalfEdge, firstUnrolledHalfEdge, secondUnrolledVertex, halfEdgeTree, mesh,
               unfoldedMesh, isFoldingEdge)
if secondHalfEdge in halfEdgeTree:
    print("sec")
    neighbourFace = mesh.face_handle(mesh.opposite_halfedge_handle(secondHalfEdge))
    unrollTree(neighbourFace, secondHalfEdge, secondUnrolledHalfEdge, thirdUnrolledVertex, halfEdgeTree, mesh,
               unfoldedMesh, isFoldingEdge)
if thirdHalfEdge in halfEdgeTree:
    print("third")
    neighbourFace = mesh.face_handle(mesh.opposite_halfedge_handle(thirdHalfEdge))
    unrollTree(neighbourFace, thirdHalfEdge, thirdUnrolledHalfEdge, firstUnrolledVertex, halfEdgeTree, mesh,
               unfoldedMesh, isFoldingEdge)

# print("Original mesh:")
# for face in mesh.faces():
#     print("Face "+ str(face.idx()))
#     for edge in mesh.fe(face):
#         print("Edge " + str(edge.idx()) + " has length " + str(mesh.calc_edge_length(edge)))
#
# print("Unfolded mesh:")
# for face in unfoldedMesh.faces():
#     print("Face "+ str(face.idx()))
#     for edge in unfoldedMesh.fe(face):
#         print("Edge " + str(edge.idx()) + " has length " + str(unfoldedMesh.calc_edge_length(edge)))

om.write_mesh('unfolding.off', unfoldedMesh)
writeSVG('unfolding.svg', unfoldedMesh, isFoldingEdge)
