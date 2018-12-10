import numpy as np
import openmesh as om
import networkx as nx
from unfoldfunctions import *
import getopt
import sys




printNumbers = True

FILENAME = 'models/icosahedron.obj'
#FILENAME = 'original.off'
#FILENAME = 'reduced.obj'
#FILENAME = 'models/reducedTeddy.obj'
#FILENAME = 'models/bunny.stl'
#FILENAME = 'models/tree.obj'
#FILENAME = 'models/polyhedron.obj'
#FILENAME = 'models/kndC.obj'

# for opt in sys.argv[1:]:
#     print(opt)
#     if opt == "-n":
#         printNumbers = True

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

    dihedralAngle = mesh.calc_dihedral_angle(edge)
    if dihedralAngle >= 0:
        dihedralAngleWeights[edge.idx()] = dihedralAngle
    else:
        dihedralAngleWeights[edge.idx()] = 10 - dihedralAngle
            #np.abs(mesh.calc_dihedral_angle(edge))

    edgeVector = mesh.calc_edge_vector(edge)
    flatAngleWeights[edge.idx()] = np.abs(np.dot(cutVector,edgeVector)) / np.linalg.norm(edgeVector)

weights = 3*minPerimeterEdgeWeights+ \
          1*flatAngleWeights + \
          0*dihedralAngleWeights
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
glueNumber = np.empty(numUnfoldedEdges, dtype = int)
#Face onnection array
connections = np.empty(numFaces, dtype=int)

# Get the points of the triangle
startingTriangle = forest[0][0]

unfold(unfoldedMesh, mesh, startingTriangle, halfEdgeTree, isFoldingEdge, connections, [], glueNumber)


#Resolve the intersection



#Find all intersections
epsilon= 1E-12
faceIntersections = []
isInterSected = np.zeros(numUnfoldedEdges, dtype=bool)
for face1 in unfoldedMesh.faces():
    for face2 in unfoldedMesh.faces():
        if face2.idx() < face1.idx():
            #Get the triangle faces
            triangle1 = []
            triangle2 = []
            for halfedge in unfoldedMesh.fh(face1):
                triangle1.append(unfoldedMesh.point(unfoldedMesh.from_vertex_handle(halfedge)))
            for halfedge in unfoldedMesh.fh(face2):
                triangle2.append(unfoldedMesh.point(unfoldedMesh.from_vertex_handle(halfedge)))
            if tri_intersect2(triangle1, triangle2, epsilon):
                #print("Intersection: " + str(face1.idx()) + " and " + str(face2.idx()))
                faceIntersections.append([connections[face1.idx()], connections[face2.idx()]])
                for edge in unfoldedMesh.fe(face1):
                    isInterSected[edge.idx()] = True
                for edge in unfoldedMesh.fe(face2):
                    isInterSected[edge.idx()] = True



print("Wir haben " + str(len(faceIntersections)) + " Überschneidungen.")

#print("Minimum: " + str(minIntersections) + " for i=" +str(minIntersectionIndex))


#Find the paths
#Make the tree a real graph
spanningGraph = nx.Graph()
for edge in spanningTree:
    face1 = mesh.face_handle(mesh.halfedge_handle(edge, 0))
    face2 = mesh.face_handle(mesh.halfedge_handle(edge, 1))
    spanningGraph.add_edge(face1.idx(), face2.idx())

print("Berechne Pfade")
paths = []
for intersection in faceIntersections:
    #Find the path
    paths.append(nx.algorithms.shortest_paths.shortest_path(spanningGraph, intersection[0], intersection[1]))
    #print(intersection[0])

#Save paths wrt edges
edgepaths = []
for path in paths:
    edgepath = []
    for i in range(len(path)-1):
        #Find edge between ith and next face
        for he in mesh.fh(mesh.face_handle(path[i])):
            if mesh.face_handle(mesh.opposite_halfedge_handle(he)) == mesh.face_handle(path[i+1]):
                edgepath.append(mesh.edge_handle(he).idx())
    edgepaths.append(edgepath)


allEdgesInPaths = list(set().union(*edgepaths))
#See how often all Edges arise in the paths
numEdgesInPaths = []
for edge in allEdgesInPaths:
    num = 0
    for path in edgepaths:
        if edge in path:
            num = num+1
    numEdgesInPaths.append(num)





S = []
C = []
gamma = 1
cutWeights = (1-gamma) * (1-weights) + gamma
print("Set covering algorithm")
while len(C) != len(paths):
    #Make the new weights of the spanningTree
    spanningWeights = np.empty(len(allEdgesInPaths))
    for i in range(len(allEdgesInPaths)):
        #Check how often the ith edge is in the path in C
        currentEdge = allEdgesInPaths[i]
        numInC = 0
        for path in C:
            if currentEdge in path:
                numInC = numInC + 1
        if (numEdgesInPaths[i] - numInC) > 0:
            spanningWeights[i] = cutWeights[edge]/(numEdgesInPaths[i] - numInC)
        else:
            spanningWeights[i] = 1000
    #Sort the edges in the paths
    minimal = np.argmin(spanningWeights)
    S.append(allEdgesInPaths[minimal])
    #Find all paths that the minimal edge is in and append to C
    for path in edgepaths:
        if allEdgesInPaths[minimal] in path and not path in C:
            C.append(path)
print("Wir müssen " + str(len(S)) + " Kanten schneiden um Überlappungen zu vermeiden.")

#Make halfedge list of cut edges
cutHalfEdges =[]
for ind in S:
    cutHalfEdges.append(mesh.halfedge_handle(mesh.edge_handle(ind), 0))
    cutHalfEdges.append(mesh.halfedge_handle(mesh.edge_handle(ind), 1))


print("Abwickeln")
#Find one face in each connected component
#Make forest
components = nx.Graph()
for edge in spanningTree:
    if not edge.idx() in S:
        face1 = mesh.face_handle(mesh.halfedge_handle(edge, 0))
        face2 = mesh.face_handle(mesh.halfedge_handle(edge, 1))
        components.add_edge(face1.idx(), face2.idx())

connectedComponents = nx.algorithms.components.connected_components(components)


count = 0
unfoldedComponents = []
for c in connectedComponents:
    startingTriangleInd = list(c)[0]
    startingTriangle = mesh.face_handle(startingTriangleInd)
    simpleUnfolded = om.TriMesh()
    simpleIsFoldingEdge = np.zeros(numUnfoldedEdges, dtype=bool)
    simpleconnections = np.empty(numFaces, dtype=int)
    simpleGlueNumber = np.empty(numUnfoldedEdges, dtype=int)

    unfold(simpleUnfolded, mesh, startingTriangle, halfEdgeTree, simpleIsFoldingEdge, simpleconnections, cutHalfEdges, simpleGlueNumber)
    #writeSVG("unfolding" + str(count) + ".svg", simpleUnfolded, simpleIsFoldingEdge, np.zeros(numUnfoldedEdges,dtype=bool), simpleGlueNumber)
    unfoldedComponents.append([simpleUnfolded, simpleIsFoldingEdge, simpleGlueNumber])
    count = count+1


#Compute maxSize of the components
maxSize = 0
for component in unfoldedComponents:
    # Get the bounding box
    firstpoint = component[0].point(mesh.vertex_handle(0))
    xmin = firstpoint[0]
    xmax = firstpoint[0]
    ymin = firstpoint[1]
    ymax = firstpoint[1]
    for vertex in component[0].vertices():
        coordinates = component[0].point(vertex)
        if (coordinates[0] < xmin):
            xmin = coordinates[0]
        if (coordinates[0] > xmax):
            xmax = coordinates[0]
        if (coordinates[1] < ymin):
            ymin = coordinates[1]
        if (coordinates[1] > ymax):
            ymax = coordinates[1]
    boxSize = np.maximum(np.abs(xmax - xmin), np.abs(ymax - ymin))
    if boxSize > maxSize:
        maxSize = boxSize


#Write
writeSVG('unfolding.svg', unfoldedMesh, isFoldingEdge, isInterSected, glueNumber, -1, printNumbers)

for i in range(len(unfoldedComponents)):
    writeSVG("unfolding" + str(i) + ".svg", unfoldedComponents[i][0], unfoldedComponents[i][1], np.zeros(numUnfoldedEdges,dtype=bool), unfoldedComponents[i][2], maxSize, printNumbers)

print("Wir haben " + str(count) + " Komponenten geschrieben.")


