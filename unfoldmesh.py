import numpy as np
import openmesh as om
import networkx as nx


#Compute the third point of a triangle when two points and all edge lengths are given
def getThirdPoint(v0, v1, l01, l12, l20):
    v2rotx = (l01 ** 2 + l20 ** 2 - l12 ** 2) / (2 * l01)
    v2roty0 = np.sqrt((l01 + l20 + l12) * (l01 + l20 - l12) * (l01 - l20 + l12) * (-l01 + l20 + l12)) / (2 * l01)

    v2roty1 = - v2roty0

    theta = np.arctan2(v1[1] - v0[1], v1[0] - v0[0])

    v2trans0 = np.array(
        [v2rotx * np.cos(theta) - v2roty0 * np.sin(theta), v2rotx * np.sin(theta) + v2roty0 * np.cos(theta), 0])
    v2trans1 = np.array(
        [v2rotx * np.cos(theta) - v2roty1 * np.sin(theta), v2rotx * np.sin(theta) + v2roty1 * np.cos(theta), 0])
    return [v2trans0 + v0, v2trans1 + v0]


#Check if two lines intersect
def lineIntersection(v1, v2, v3, v4, epsilon):
    d = (v4[1] - v3[1]) * (v2[0] - v1[0]) - (v4[0] - v3[0]) * (v2[1] - v1[1])
    u = (v4[0] - v3[0]) * (v1[1] - v3[1]) - (v4[1] - v3[1]) * (v1[0] - v3[0])
    v = (v2[0] - v1[0]) * (v1[1] - v3[1]) - (v2[1] - v1[1]) * (v1[0] - v3[0])
    if d < 0:
        u, v, d = -u, -v, -d
    return ((0 + epsilon) <= u <= (d - epsilon)) and ((0 + epsilon) <= v <= (d - epsilon))

#Check if a point lies inside a triangle
def pointInTriangle(A, B, C, P, epsilon):
    v0 = [C[0] - A[0], C[1] - A[1]]
    v1 = [B[0] - A[0], B[1] - A[1]]
    v2 = [P[0] - A[0], P[1] - A[1]]
    cross = lambda u, v: u[0] * v[1] - u[1] * v[0]
    u = cross(v2, v0)
    v = cross(v1, v2)
    d = cross(v1, v0)
    if d < 0:
        u, v, d = -u, -v, -d
    return u >= (0 + epsilon) and v >= (0 + epsilon) and (u + v) <= (d - epsilon)


#Check if two triangles intersect
def triangleIntersection(t1, t2, epsilon):
    if lineIntersection(t1[0], t1[1], t2[0], t2[1], epsilon): return True
    if lineIntersection(t1[0], t1[1], t2[0], t2[2], epsilon): return True
    if lineIntersection(t1[0], t1[1], t2[1], t2[2], epsilon): return True
    if lineIntersection(t1[0], t1[2], t2[0], t2[1], epsilon): return True
    if lineIntersection(t1[0], t1[2], t2[0], t2[2], epsilon): return True
    if lineIntersection(t1[0], t1[2], t2[1], t2[2], epsilon): return True
    if lineIntersection(t1[1], t1[2], t2[0], t2[1], epsilon): return True
    if lineIntersection(t1[1], t1[2], t2[0], t2[2], epsilon): return True
    if lineIntersection(t1[1], t1[2], t2[1], t2[2], epsilon): return True
    inTri = True
    inTri = inTri and pointInTriangle(t1[0], t1[1], t1[2], t2[0], epsilon)
    inTri = inTri and pointInTriangle(t1[0], t1[1], t1[2], t2[1], epsilon)
    inTri = inTri and pointInTriangle(t1[0], t1[1], t1[2], t2[2], epsilon)
    if inTri == True: return True
    inTri = True
    inTri = inTri and pointInTriangle(t2[0], t2[1], t2[2], t1[0], epsilon)
    inTri = inTri and pointInTriangle(t2[0], t2[1], t2[2], t1[1], epsilon)
    inTri = inTri and pointInTriangle(t2[0], t2[1], t2[2], t1[2], epsilon)
    if inTri == True: return True
    return False


# Funktionen für die Visualisierung und Output
def addVisualisationData(mesh, unfoldedMesh, originalHalfedges, unfoldedHalfedges, glueNumber, foldingDirection):
    for i in range(3):
        # Faltungsrichtung
        if mesh.calc_dihedral_angle(originalHalfedges[i]) < 0:
            foldingDirection[unfoldedMesh.edge_handle(unfoldedHalfedges[i]).idx()] = -1
        else:
            foldingDirection[unfoldedMesh.edge_handle(unfoldedHalfedges[i]).idx()] = 1

        # Information, welche Kanten zusammengehören
        glueNumber[unfoldedMesh.edge_handle(unfoldedHalfedges[i]).idx()] = mesh.edge_handle(originalHalfedges[i]).idx()

# Funktion, die einen aufspannenden Baum abwickelt
def unfoldSpanningTree(mesh, spanningTree):
    unfoldedMesh = om.TriMesh()  # Das abgewickelte Netz

    numFaces = mesh.n_faces()
    sizeTree = spanningTree.number_of_edges()
    numUnfoldedEdges = 3 * numFaces - sizeTree

    isFoldingEdge = np.zeros(numUnfoldedEdges, dtype=bool)  # Gibt an, ob eine Kante gefaltet oder geschnitten wird
    glueNumber = np.empty(numUnfoldedEdges, dtype=int)  # Speichert, mit welcher Kante zusammegeklebt wird
    foldingDirection = np.empty(numUnfoldedEdges, dtype=int)  # Talfaltung oder Bergfaltung

    connections = np.empty(numFaces, dtype=int)  # Speichert, welches ursprüngliche Dreieck zum abgewickelten gehört

    # Wähle das erste Dreieck beliebig
    startingNode = list(spanningTree.nodes())[0]
    startingTriangle = mesh.face_handle(startingNode)

    # Wir wickeln das erste Dreieck ab

    # Alle Halbkanten des ersten Dreiecks
    firstHalfEdge = mesh.halfedge_handle(startingTriangle)
    secondHalfEdge = mesh.next_halfedge_handle(firstHalfEdge)
    thirdHalfEdge = mesh.next_halfedge_handle(secondHalfEdge)
    originalHalfEdges = [firstHalfEdge, secondHalfEdge, thirdHalfEdge]

    # Berechne die Längen der Kanten, hierdurch wird die Form des Dreiecks bestimmt (Kongruenz)
    edgelengths = [mesh.calc_edge_length(firstHalfEdge), mesh.calc_edge_length(secondHalfEdge),
                   mesh.calc_edge_length(thirdHalfEdge)]

    # Die beiden ersten Punkte
    firstUnfoldedPoint = np.array([0, 0, 0])
    secondUnfoldedPoint = np.array([edgelengths[0], 0, 0])

    # Wir berechnen den dritten Punkt des Dreiecks aus den ersten beiden. Es gibt zwei möglichkeiten
    [thirdUnfolded0, thirdUnfolded1] = getThirdPoint(firstUnfoldedPoint, secondUnfoldedPoint, edgelengths[0],
                                                     edgelengths[1],
                                                     edgelengths[2])
    if thirdUnfolded0[1] > 0:
        thirdUnfoldedPoint = thirdUnfolded0
    else:
        thirdUnfoldedPoint = thirdUnfolded1

    # Füge die neuen Ecken zum abgewickelten Netz hinzu
    # firstUnfoldedVertex = unfoldedMesh.add_vertex(secondUnfoldedPoint)
    # secondUnfoldedVertex = unfoldedMesh.add_vertex(thirdUnfoldedPoint)
    # thirdUnfoldedVertex = unfoldedMesh.add_vertex(firstUnfoldedPoint)

    firstUnfoldedVertex = unfoldedMesh.add_vertex(firstUnfoldedPoint)
    secondUnfoldedVertex = unfoldedMesh.add_vertex(secondUnfoldedPoint)
    thirdUnfoldedVertex = unfoldedMesh.add_vertex(thirdUnfoldedPoint)

    # Erzeuge die Seite
    unfoldedFace = unfoldedMesh.add_face(firstUnfoldedVertex, secondUnfoldedVertex, thirdUnfoldedVertex)

    # Speichere Eigenschaften der Fläche und Kanten
    # Die Halbkanten im abgewickelten Gitter
    firstUnfoldedHalfEdge = unfoldedMesh.next_halfedge_handle(unfoldedMesh.opposite_halfedge_handle(unfoldedMesh.halfedge_handle(firstUnfoldedVertex)))
    secondUnfoldedHalfEdge = unfoldedMesh.next_halfedge_handle(firstUnfoldedHalfEdge)
    thirdUnfoldedHalfEdge = unfoldedMesh.next_halfedge_handle(secondUnfoldedHalfEdge)

    unfoldedHalfEdges = [firstUnfoldedHalfEdge, secondUnfoldedHalfEdge, thirdUnfoldedHalfEdge]

    # Zugehöriges Dreieck im 3D-Netz
    connections[unfoldedFace.idx()] = startingTriangle.idx()
    # Faltungsrichtung und Klebenummer
    addVisualisationData(mesh, unfoldedMesh, originalHalfEdges, unfoldedHalfEdges, glueNumber, foldingDirection)

    halfEdgeConnections = {firstHalfEdge.idx(): firstUnfoldedHalfEdge.idx(),
                           secondHalfEdge.idx(): secondUnfoldedHalfEdge.idx(),
                           thirdHalfEdge.idx(): thirdUnfoldedHalfEdge.idx()}

    # Wir gehen durch den Baum
    for dualEdge in nx.dfs_edges(spanningTree, source=startingNode):
        foldingEdge = mesh.edge_handle(spanningTree[dualEdge[0]][dualEdge[1]]['idx'])
        # Finde die dazugehörige Halbkante im AusgangsDreieck
        foldingHalfEdge = mesh.halfedge_handle(foldingEdge, 0)
        if not (mesh.face_handle(foldingHalfEdge).idx() == dualEdge[0]):
            foldingHalfEdge = mesh.halfedge_handle(foldingEdge, 1)

        # Finde die dazugehörige abgewickelte Halbkante
        unfoldedLastHalfEdge = unfoldedMesh.halfedge_handle(halfEdgeConnections[foldingHalfEdge.idx()])

        # Finde den Punkt im Abgwickelten Dreieck, der nicht auf der Faltkante liegt
        oppositeUnfoldedVertex = unfoldedMesh.to_vertex_handle(unfoldedMesh.next_halfedge_handle(unfoldedLastHalfEdge))

        # Wir drehen die Halbkanten um, um im neuen Dreieck zu liegen
        foldingHalfEdge = mesh.opposite_halfedge_handle(foldingHalfEdge)
        unfoldedLastHalfEdge = unfoldedMesh.opposite_halfedge_handle(unfoldedLastHalfEdge)

        # Die beiden Ecken der Faltkante
        unfoldedFromVertex = unfoldedMesh.from_vertex_handle(unfoldedLastHalfEdge)
        unfoldedToVertex = unfoldedMesh.to_vertex_handle(unfoldedLastHalfEdge)

        # Berechne die Kantenlängen im neuen Dreieck
        secondHalfEdgeInFace = mesh.next_halfedge_handle(foldingHalfEdge)
        thirdHalfEdgeInFace = mesh.next_halfedge_handle(secondHalfEdgeInFace)

        originalHalfEdges = [foldingHalfEdge, secondHalfEdgeInFace, thirdHalfEdgeInFace]

        edgelengths = [mesh.calc_edge_length(foldingHalfEdge), mesh.calc_edge_length(secondHalfEdgeInFace),
                       mesh.calc_edge_length(thirdHalfEdgeInFace)]

        # Wir berechnen die beiden Möglichkeiten für den dritten Punkt im Dreieck
        [newUnfoldedVertex0, newUnfoldedVertex1] = getThirdPoint(unfoldedMesh.point(unfoldedFromVertex),
                                                                 unfoldedMesh.point(unfoldedToVertex), edgelengths[0],
                                                                 edgelengths[1], edgelengths[2])


        newUnfoldedVertex = unfoldedMesh.add_vertex(newUnfoldedVertex0)

        # Make the face
        newface = unfoldedMesh.add_face(unfoldedFromVertex, unfoldedToVertex, newUnfoldedVertex)

        secondUnfoldedHalfEdge = unfoldedMesh.next_halfedge_handle(unfoldedLastHalfEdge)
        thirdUnfoldedHalfEdge = unfoldedMesh.next_halfedge_handle(secondUnfoldedHalfEdge)
        unfoldedHalfEdges = [unfoldedLastHalfEdge, secondUnfoldedHalfEdge, thirdUnfoldedHalfEdge]

        # Speichern der Informationen über Kanten und Seite
        # Gestrichelte Linie in der Ausgabe
        unfoldedLastEdge = unfoldedMesh.edge_handle(unfoldedLastHalfEdge)
        isFoldingEdge[unfoldedLastEdge.idx()] = True

        # Klebenummer und Faltrichtung
        addVisualisationData(mesh, unfoldedMesh, originalHalfEdges, unfoldedHalfEdges, glueNumber, foldingDirection)

        # Zugehörige Seite
        connections[newface.idx()] = dualEdge[1]

        # Identifiziere die Halbkanten
        for i in range(3):
            halfEdgeConnections[originalHalfEdges[i].idx()] = unfoldedHalfEdges[i].idx()

    return [unfoldedMesh, isFoldingEdge, connections, glueNumber, foldingDirection]

def unfold(mesh):
    # Berechne die Anzahl der Flächen, Kanten und Ecken, sowie die Längen der längsten kürzesten Kante
    numEdges = mesh.n_edges()
    numVertices = mesh.n_vertices()
    numFaces = mesh.n_faces()

    # Erzeuge den dualen Graphen des Netzes und berechne die Gewichte
    dualGraph = nx.Graph()

    # Für die Gewichte: Berechne die längste und kürzeste Kante des Dreiecks
    minLength = 1000
    maxLength = 0
    for edge in mesh.edges():
        edgelength = mesh.calc_edge_length(edge)
        if edgelength < minLength:
            minLength = edgelength
        if edgelength > maxLength:
            maxLength = edgelength

    # Alle Kanten im Netz
    for edge in mesh.edges():
        # Die beiden Seiten, die an die Kante angrenzebn
        face1 = mesh.face_handle(mesh.halfedge_handle(edge, 0))
        face2 = mesh.face_handle(mesh.halfedge_handle(edge, 1))

        # Das Gewicht
        edgeweight = 1.0 - (mesh.calc_edge_length(edge) - minLength) / (maxLength - minLength)

        # Berechne die Mittelpunkte der Seiten (nur für die Visualisierung notwendig)
        center1 = (0, 0)
        for vertex in mesh.fv(face1):
            center1 = center1 + 0.3333333333333333 * np.array([mesh.point(vertex)[0], mesh.point(vertex)[2]])
        center2 = (0, 0)
        for vertex in mesh.fv(face2):
            center2 = center2 + 0.3333333333333333 * np.array([mesh.point(vertex)[0], mesh.point(vertex)[2]])

        # Füge die neuen Knoten und Kante zum dualen Graph hinzu
        dualGraph.add_node(face1.idx(), pos=center1)
        dualGraph.add_node(face2.idx(), pos=center2)
        dualGraph.add_edge(face1.idx(), face2.idx(), idx=edge.idx(), weight=edgeweight)

    # Berechne den minimalen aufspannenden Baum
    spanningTree = nx.minimum_spanning_tree(dualGraph)

    #Unfold the tree

    fullUnfolding = unfoldSpanningTree(mesh, spanningTree)
    [unfoldedMesh, isFoldingEdge, connections, glueNumber, foldingDirection] = fullUnfolding


    # Resolve the intersections
    # Find all intersections
    epsilon = 1E-12  # Genauigkeit
    faceIntersections = []
    for face1 in unfoldedMesh.faces():
        for face2 in unfoldedMesh.faces():
            if face2.idx() < face1.idx():  # Damit wir die Paare nicht doppelt durchgehen
                # Get the triangle faces
                triangle1 = []
                triangle2 = []
                for halfedge in unfoldedMesh.fh(face1):
                    triangle1.append(unfoldedMesh.point(unfoldedMesh.from_vertex_handle(halfedge)))
                for halfedge in unfoldedMesh.fh(face2):
                    triangle2.append(unfoldedMesh.point(unfoldedMesh.from_vertex_handle(halfedge)))
                if triangleIntersection(triangle1, triangle2, epsilon):
                    faceIntersections.append([connections[face1.idx()], connections[face2.idx()]])


    # Find the paths
    # Wir finden die minimale Anzahl von Schnitten, um jede Selbstüberschneidung aufzulösen

    # Suche alle Pfade zwischen sich überschneidenden Dreiecken
    paths = []
    for intersection in faceIntersections:
        paths.append(
            nx.algorithms.shortest_paths.shortest_path(spanningTree, source=intersection[0], target=intersection[1]))

    # Finde alle Kanten in allen Pfäden
    edgepaths = []
    for path in paths:
        edgepath = []
        for i in range(len(path) - 1):
            edgepath.append((path[i], path[i + 1]))
        edgepaths.append(edgepath)

    # Liste aller Kanten in allen Pfaden
    allEdgesInPaths = list(set().union(*edgepaths))

    # Zähle, wie oft jede Kante vorkommt
    numEdgesInPaths = []
    for edge in allEdgesInPaths:
        num = 0
        for path in edgepaths:
            if edge in path:
                num = num + 1
        numEdgesInPaths.append(num)

    S = []
    C = []

    while len(C) != len(paths):
        # Berechne die Gewichte um zu entscheiden, welche Kante wir schneiden
        cutWeights = np.empty(len(allEdgesInPaths))
        for i in range(len(allEdgesInPaths)):
            currentEdge = allEdgesInPaths[i]

            # Zähle, wievele der Pfade, in denen die Kante vorkommt, bereits geschnitten wurden
            numInC = 0
            for path in C:
                if currentEdge in path:
                    numInC = numInC + 1

            # Bestimme das Gewicht
            if (numEdgesInPaths[i] - numInC) > 0:
                cutWeights[i] = 1 / (numEdgesInPaths[i] - numInC)
            else:
                cutWeights[i] = 1000  # 1000 = unendlich
        # Finde die Kante mit dem kleinsten Gewicht
        minimalIndex = np.argmin(cutWeights)
        S.append(allEdgesInPaths[minimalIndex])
        # Finde alle Pfade, in denen die Kante vorkommt und füge sie zu C hinzu
        for path in edgepaths:
            if allEdgesInPaths[minimalIndex] in path and not path in C:
                C.append(path)

    # Nun entfernen wir die Schnittkanten aus dem minimalen aufspannenden Baum
    spanningTree.remove_edges_from(S)

    # Finde die Zusammenhangskomponenten
    connectedComponents = nx.algorithms.components.connected_components(spanningTree)
    connectedComponentList = list(connectedComponents)

    # Abwicklung der Komponenten
    unfoldings = []
    for component in connectedComponentList:
        unfoldings.append(unfoldSpanningTree(mesh, spanningTree.subgraph(component)))

    return fullUnfolding, unfoldings


def findBoundingBox(mesh):
    firstpoint = mesh.point(mesh.vertex_handle(0))
    xmin = firstpoint[0]
    xmax = firstpoint[0]
    ymin = firstpoint[1]
    ymax = firstpoint[1]
    for vertex in mesh.vertices():
        coordinates = mesh.point(vertex)
        if (coordinates[0] < xmin):
            xmin = coordinates[0]
        if (coordinates[0] > xmax):
            xmax = coordinates[0]
        if (coordinates[1] < ymin):
            ymin = coordinates[1]
        if (coordinates[1] > ymax):
            ymax = coordinates[1]
    boxSize = np.maximum(np.abs(xmax - xmin), np.abs(ymax - ymin))

    return [xmin, ymin, boxSize]


def writeSVG(filename, unfolding, size, printNumbers):
    mesh = unfolding[0]
    isFoldingEdge = unfolding[1]
    glueNumber = unfolding[3]
    foldingDirection = unfolding[4]

    # Berechne die bounding box
    [xmin, ymin, boxSize] = findBoundingBox(unfolding[0])

    if size > 0:
        boxSize = size

    strokewidth = 0.002 * boxSize
    dashLength = 0.008 * boxSize
    spaceLength = 0.02 * boxSize

    textDistance = 0.02 * boxSize
    textStrokewidth = 0.05 * strokewidth
    textLength = 0.001 * boxSize
    fontsize = 0.015 * boxSize

    frame = 0.05 * boxSize

    # Öffne Datei im Schreibmodus (write)
    file = open(filename, 'w')

    # Schreibe xml-header
    file.write("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n")

    # Die Papiergröße und die Skalierung
    file.write("<svg width=\"30.5cm\" height=\"30.5cm\" viewBox = \"" + str(xmin - frame) + " " + str(
        ymin - frame) + " " + str(boxSize + 2 * frame) + " " + str(
        boxSize + 2 * frame) + "\" version = \"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n")

    # Gehe über alle Kanten des GItters
    for edge in mesh.edges():
        # Die beiden Endpunkte
        he = mesh.halfedge_handle(edge, 0)
        vertex0 = mesh.point(mesh.from_vertex_handle(he))
        vertex1 = mesh.point(mesh.to_vertex_handle(he))

        # Schreibe eine Gerade zwischen den beiden Ecken
        file.write("<path d =\"M " + str(vertex0[0]) + "," + str(vertex0[1]) + " " + str(vertex1[0]) + "," + str(
            vertex1[1]) + "\" style=\"fill:none;stroke:")

        # Farbe je nach Faltrichtung
        if foldingDirection[edge.idx()] > 0:
            file.write("#ff0000")
        elif foldingDirection[edge.idx()] < 0:
            file.write("#0066ff")

        file.write(";stroke-width:" + str(
            strokewidth) + ";stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:4;stroke-dasharray:")

        # Gestrichelte Linien für Faltkanten
        if isFoldingEdge[edge.idx()]:
            file.write((str(dashLength) + ", " + str(spaceLength)))
        else:
            file.write("none")

        file.write(";stroke-dashoffset:0;stroke-opacity:1")
        file.write("\" />\n")

        # Die Nummer der Kante, mit der zusammengeklebt wird
        if not isFoldingEdge[edge.idx()]:
            # Find halfedge in the face
            halfEdge = mesh.halfedge_handle(edge, 0)
            if mesh.face_handle(halfEdge).idx() == -1:
                halfEdge = mesh.opposite_halfedge_handle(halfEdge)
            vector = mesh.calc_edge_vector(halfEdge)
            # normalize
            vector = vector / np.linalg.norm(vector)
            midPoint = 0.5 * (
                    mesh.point(mesh.from_vertex_handle(halfEdge)) + mesh.point(mesh.to_vertex_handle(halfEdge)))
            rotatedVector = np.array([-vector[1], vector[0], 0])
            angle = np.arctan2(vector[1], vector[0])
            position = midPoint + textDistance * rotatedVector
            rotation = 180 / np.pi * angle

            if (printNumbers):
                file.write("<text x=\"" + str(position[0]) + "\" y=\"" + str(position[1]) + "\" font-size=\"" + str(
                    fontsize) + "\" stroke-width=\"" + str(textStrokewidth) + "\" transform=\"rotate(" + str(
                    rotation) + "," + str(position[0]) + "," + str(position[1]) + ")\">" + str(
                    glueNumber[edge.idx()]) + "</text>\n")

    file.write("</svg>")
    file.close()