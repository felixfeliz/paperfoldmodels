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

def unrollTree(vertex, last, tree):
    neighbours = []
    for edge in tree:
        dualVertex1 = mesh.face_handle(mesh.halfedge_handle(edge, 0))
        dualVertex2 = mesh.face_handle(mesh.halfedge_handle(edge, 1))
        if (vertex == dualVertex1) and (dualVertex2 != last):
            neighbours.append(dualVertex2)
        elif (vertex == dualVertex2) and (dualVertex1 != last):
            neighbours.append(dualVertex1)
    for face in neighbours:
        print(face.idx())
    for face in neighbours:
        unrollTree(face, vertex, tree)

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

#Go through the tree
unfoldedMesh = om.TriMesh()


#Get the points of the triangle
trianglePoints = np.empty([3,3])
i = 0
for vh in mesh.fv(forest[0][0]):
    trianglePoints[i,:] = mesh.point(vh)
    i = i+1
l01 = np.linalg.norm(trianglePoints[1,:] - trianglePoints[0,:])
l12 = np.linalg.norm(trianglePoints[2,:] - trianglePoints[1,:])
l20 = np.linalg.norm(trianglePoints[0,:] - trianglePoints[2,:])
#The orientation of the first triangle is arbitrary
vh0 = unfoldedMesh.add_vertex([0, 0, 0])
vh1 = unfoldedMesh.add_vertex([l01, 0, 0])
#Compute third point from lengths


#print(trianglePoints)
#print(forest[0][0].idx())
#unrollTree(forest[0][0], forest[0][0], spanningTree)
