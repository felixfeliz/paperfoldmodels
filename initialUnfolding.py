import numpy as np
import openmesh as om

FILENAME = 'reduced.off'
FILENAME = 'models/icosahedron.obj'
#FILENAME = 'original.off'

# Import the mode
mesh = om.read_trimesh(FILENAME)

# Berechne die Anzahl der Flächen, Kanten und Ecken, sowie die Längen der längsten kürzesten Kante
numEdges = 0
numVertices = 0
numFaces = 0
minLength = 1000
maxLength = 0
for edge in mesh.edges():
    numEdges = numEdges + 1
    edgelength = mesh.calc_edge_length(edge)
    if edgelength < minLength:
        minLength = edgelength
    if edgelength > maxLength:
        maxLength = edgelength
for vertex in mesh.vertices():
    numVertices = numVertices + 1
for face in mesh.faces():
    numFaces = numFaces + 1

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

