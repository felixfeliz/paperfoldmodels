import numpy as np
import openmesh as om

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

def unrollTree(face, lastHalfEdge, unrolledLastHalfEdge, oppositeUnrolledVertex, halfEdgeTree, mesh, unfoldedMesh):
    print("Processing face " + str(face.idx()))
    #First unroll the current face
    #Get the two known unrolled Vertices
    unrolledFromVertex = unfoldedMesh.to_vertex_handle(unrolledLastHalfEdge)
    unrolledToVertex = unfoldedMesh.from_vertex_handle(unrolledLastHalfEdge)

    #Get the edge lengths
    lastHalfEdgeInFace = mesh.opposite_halfedge_handle(lastHalfEdge)
    secondHalfEdgeInFace = mesh.next_halfedge_handle(lastHalfEdgeInFace)
    thirdHalfEdgeInFace = mesh.next_halfedge_handle(secondHalfEdgeInFace)
    edgelengths=[mesh.calc_edge_length(lastHalfEdgeInFace), mesh.calc_edge_length(secondHalfEdgeInFace),mesh.calc_edge_length(thirdHalfEdgeInFace)]

    #Find the third unrolled Point
    #print([unfoldedMesh.point(unrolledFromVertex),unfoldedMesh.point(unrolledToVertex), edgelengths ])
    [newUnrolledVertex0, newUnrolledVertex1] = getThirdPoint(unfoldedMesh.point(unrolledFromVertex), unfoldedMesh.point(unrolledToVertex), edgelengths[0],
                               edgelengths[1], edgelengths[2])
    #Check which one is on the opposite side of the edge from  oppositeUnrolledVertex
    lastUnrolledEdgeVector = unfoldedMesh.point(unrolledToVertex) - unfoldedMesh.point(unrolledFromVertex)
    newUnrolledVertexVector0 = newUnrolledVertex0- unfoldedMesh.point(unrolledFromVertex)

    newUnrolledVertexVector1 = newUnrolledVertex1 - unfoldedMesh.point(unrolledFromVertex)
    oppositeVector = unfoldedMesh.point(oppositeUnrolledVertex) - unfoldedMesh.point(unrolledFromVertex)

    if  np.dot(np.cross(lastUnrolledEdgeVector, newUnrolledVertexVector0), np.cross(lastUnrolledEdgeVector, oppositeVector)) <= 0:
        newUnrolledVertex = unfoldedMesh.add_vertex(newUnrolledVertex0)
    else:
        newUnrolledVertex = unfoldedMesh.add_vertex(newUnrolledVertex1)

    #Make the face
    unfoldedMesh.add_face(unrolledFromVertex, unrolledToVertex, newUnrolledVertex)

    print("Unfolded face " + str(face.idx()) + ".")
    om.write_mesh('unfolding' + str(face.idx()) + ".off", unfoldedMesh)

    #Now find the neighbours
    # Get the unrolled halfEdge
    secondUnrolledHalfEdge = unfoldedMesh.next_halfedge_handle(unfoldedMesh.opposite_halfedge_handle(unrolledLastHalfEdge))
    thirdUnrolledHalfEdge = unfoldedMesh.next_halfedge_handle(secondUnrolledHalfEdge)

    #Check the two other half edges
    if secondHalfEdgeInFace in halfEdgeTree:
        #Get the face
        neighbourFace = mesh.face_handle(mesh.opposite_halfedge_handle(secondHalfEdgeInFace))
        unrollTree(neighbourFace, secondHalfEdgeInFace, secondUnrolledHalfEdge,unrolledFromVertex,halfEdgeTree,mesh,unfoldedMesh)
    if thirdHalfEdgeInFace in halfEdgeTree:
        #Get the face
        neighbourFace = mesh.face_handle(mesh.opposite_halfedge_handle(thirdHalfEdgeInFace))
        unrollTree(neighbourFace, thirdHalfEdgeInFace, thirdUnrolledHalfEdge,unrolledToVertex,halfEdgeTree,mesh,unfoldedMesh)


def findLeafIndex(forest, spanningTree, mesh):
    for i in range(len(forest)):
        outerEdges = []

        for edge in spanningTree:
            dualVertex1 = mesh.face_handle(mesh.halfedge_handle(edge, 0))
            dualVertex2 = mesh.face_handle(mesh.halfedge_handle(edge, 1))
            if (dualVertex1 in forest[i]) != (dualVertex2 in forest[i]):
                outerEdges.append(edge)

        if len(outerEdges) == 1:
            #Find the parent
            dualVertex1 = mesh.face_handle(mesh.halfedge_handle(outerEdges[0], 0))
            dualVertex2 = mesh.face_handle(mesh.halfedge_handle(outerEdges[0], 1))
            if dualVertex1 in forest[i]:
                parent = dualVertex1
            else:
                parent = dualVertex2
            return [True, i, outerEdges[0], parent]
    return [False, -1, -1, -1]

def getThirdPoint(v0,v1,l01,l12,l20):
    v2rotx = (l01**2 + l20*2 - l12**2) / (2 * l01)
    v2roty0 = np.sqrt((l01 + l20 + l12) * (l01 + l20 - l12) * (l01 - l20 + l12) * (-l01 + l20 +l12)) / (2*l01)
    v2roty1 = - v2roty0

    theta = np.arctan2(v1[1] - v0[1], v1[0] - v0[0])
    v2trans0 = np.array([v2rotx * np.cos(theta) - v2roty0 *np.sin(theta), v2rotx * np.sin(theta) + v2roty0*np.cos(theta),0])
    v2trans1 = np.array([v2rotx * np.cos(theta) - v2roty1 *np.sin(theta), v2rotx * np.sin(theta) + v2roty1*np.cos(theta),0])

    return[v2trans0 + v1, v2trans1 + v1]

#FILENAME = 'reduced.off'
FILENAME = 'models/icosahedron.obj'
#FILENAME = 'original.off'

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
for edge in mesh.edges():
    edgelength = mesh.calc_edge_length(edge)
    #Schneide kürzere Kanten um Umfang zu minimieren
    minPerimeterEdgeWeights[edge.idx()] = 1.0 - (edgelength - minLength) / (maxLength - minLength)


#Kruzkal's algorithm

#Sort the weights and find the order of the indices
sortedIndices = np.argsort(minPerimeterEdgeWeights)
#Make ordered list of (dual) edges
sorteddualEdges = []
for i in range(numEdges):
    sorteddualEdges.append(mesh.edge_handle(sortedIndices[i]))

#Initialize the forest with all dual-vertices (i.e. faces)
spanningTree = []
forest = []
for face in mesh.faces():
    forest.append([face])
#Perform the step as long as there is more than one tree in the forest
i = 0
while len(forest) > 1:
    # Get the two faces adjacent to the ith edge
    dualVertex1 = mesh.face_handle(mesh.halfedge_handle(sorteddualEdges[i], 0))
    dualVertex2 = mesh.face_handle(mesh.halfedge_handle(sorteddualEdges[i], 1))

    #Check if both are in different trees
    index1 = findTree(dualVertex1, forest)
    index2 = findTree(dualVertex2, forest)
    if index1 != index2:
        #Merge the two trees
        for element in forest[index2]:
            forest[index1].append(element)
        forest.pop(index2)

        #Add the faces and vertices to the spanning tree

        spanningTree.append(sorteddualEdges[i])

    i = i+1

sizeTree = 0
for edge in spanningTree:
    sizeTree = sizeTree + 1
print("Der Baum hat " + str(sizeTree) + " Kanten.")

#Make a tree of HalfEdges
halfEdgeTree = []
for edge in spanningTree:
    halfEdgeTree.append(mesh.halfedge_handle(edge,0))
    halfEdgeTree.append(mesh.halfedge_handle(edge, 1))

#Unfolding
unfoldedMesh = om.TriMesh()

#Get the points of the triangle
startingTriangle = forest[0][0]

#Get all halfedges
firstHalfEdge = mesh.halfedge_handle(startingTriangle)
secondHalfEdge = mesh.next_halfedge_handle(firstHalfEdge)
thirdHalfEdge = mesh.next_halfedge_handle(secondHalfEdge)

edgelengths = [mesh.calc_edge_length(firstHalfEdge), mesh.calc_edge_length(secondHalfEdge), mesh.calc_edge_length(thirdHalfEdge)]

#The orientation of the first triangle is arbitrary
firstUnrolled = np.array([0,0,0])
secondUnrolled = np.array([edgelengths[0],0,0])

#Compute third point from lengths
[thirdUnrolled0, thirdUnrolled1] = getThirdPoint(firstUnrolled, secondUnrolled, edgelengths[0], edgelengths[1], edgelengths[2])
#Make triangle
#Add vertices
firstUnrolledVertex = unfoldedMesh.add_vertex(firstUnrolled)
secondUnrolledVertex = unfoldedMesh.add_vertex(secondUnrolled)
thirdUnrolledVertex = unfoldedMesh.add_vertex(thirdUnrolled0)
#Create the face
f = unfoldedMesh.add_face(firstUnrolledVertex,secondUnrolledVertex,thirdUnrolledVertex)
om.write_mesh('unfolding' + str(startingTriangle.idx()) + ".off", unfoldedMesh)

#Now check the neighbours
#Find the unrolled half-edges inside the face
firstUnrolledHalfEdge = unfoldedMesh.opposite_halfedge_handle(unfoldedMesh.halfedge_handle(firstUnrolledVertex))
secondUnrolledHalfEdge = unfoldedMesh.opposite_halfedge_handle(unfoldedMesh.next_halfedge_handle(firstUnrolledHalfEdge))
thirdUnrolledHalfEdge = unfoldedMesh.opposite_halfedge_handle(unfoldedMesh.next_halfedge_handle(secondUnrolledHalfEdge))

if firstHalfEdge in halfEdgeTree:
    #Get the neighbouring face
    neighbourFace = mesh.face_handle(mesh.opposite_halfedge_handle(firstHalfEdge))
    unrollTree(neighbourFace, firstHalfEdge, firstUnrolledHalfEdge,secondUnrolledVertex, halfEdgeTree, mesh, unfoldedMesh)
if secondHalfEdge in halfEdgeTree:
    neighbourFace = mesh.face_handle(mesh.opposite_halfedge_handle(secondHalfEdge))
    unrollTree(neighbourFace, secondHalfEdge, secondUnrolledHalfEdge, thirdUnrolledVertex, halfEdgeTree, mesh, unfoldedMesh)
if thirdHalfEdge in halfEdgeTree:
    neighbourFace = mesh.face_handle(mesh.opposite_halfedge_handle(thirdHalfEdge))
    unrollTree(neighbourFace, thirdHalfEdge, thirdUnrolledHalfEdge, firstUnrolledVertex, halfEdgeTree, mesh, unfoldedMesh)


#Connection dict
#connections = {vertices[0].idx() : vh0, vertices[1].idx() : vh1, vertices[2].idx() : vh2}

#
# neighbours = []
# for edge in spanningTree:
#     dualVertex1 = mesh.face_handle(mesh.halfedge_handle(edge, 0))
#     dualVertex2 = mesh.face_handle(mesh.halfedge_handle(edge, 1))
#     if (startingTriangle == dualVertex1):
#         neighbours.append([dualVertex2,edge])
#     elif (startingTriangle == dualVertex2):
#         neighbours.append([dualVertex1, edge])
# for element in neighbours:
#     # Find the half edge which is in the face
#     halfEdge = mesh.halfedge_handle(element[1], 0)
#     if mesh.face_handle(halfEdge) != startingTriangle:
#         halfEdge = mesh.opposite_halfedge_handle(halfEdge)
#
#     firstVertex = mesh.from_vertex_handle(halfEdge)
#     secondVertex = mesh.to_vertex_handle(halfEdge)
#     #Find the third vertex
#     for vertex in vertices:
#         if vertex != firstVertex and vertex != secondVertex:
#             thirdVertex = vertex
#     unrollTree(element[0], element[1], connections[firstVertex.idx()], connections[secondVertex.idx()], connections[thirdVertex.idx()], spanningTree, mesh, unfoldedMesh)

om.write_mesh('unfolding.off',unfoldedMesh)


