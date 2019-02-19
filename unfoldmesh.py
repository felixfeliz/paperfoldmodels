import numpy as np
import openmesh as om
import networkx as nx


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


def unfoldNeighborsInTree(face, lastHalfEdge, unrolledLastHalfEdge, oppositeUnrolledVertex, halfEdgeTree, mesh, unfoldedMesh,
                          isFoldingEdge, connections, cutEdges, glueNumber, foldingDirection):
    # print("Processing face " + str(face.idx()))

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
    newface = unfoldedMesh.add_face(unrolledFromVertex, unrolledToVertex, newUnrolledVertex)

    # Make the last edge a folding edge
    unrolledLastEdge = unfoldedMesh.edge_handle(unrolledLastHalfEdge)
    isFoldingEdge[unrolledLastEdge.idx()] = True



    # Save the connections
    connections[newface.idx()] = face.idx()

    # print("Unfolded face " + str(face.idx()) + ".")
    # om.write_mesh('unfolding' + str(face.idx()) + ".off", unfoldedMesh)

    # Now find the neighbours
    # Get the unrolled halfEdge
    secondUnrolledHalfEdge = unfoldedMesh.next_halfedge_handle(
        unfoldedMesh.opposite_halfedge_handle(unrolledLastHalfEdge))
    thirdUnrolledHalfEdge = unfoldedMesh.next_halfedge_handle(secondUnrolledHalfEdge)

    # Compute the angle to get the folding directions
    if mesh.calc_dihedral_angle(lastHalfEdge) < 0:
        foldingDirection[unrolledLastEdge.idx()] = -1
    else:
        foldingDirection[unrolledLastEdge.idx()] = 1
    if mesh.calc_dihedral_angle(secondHalfEdgeInFace) < 0:
        foldingDirection[unfoldedMesh.edge_handle(secondUnrolledHalfEdge).idx()] = -1
    else:
        foldingDirection[unfoldedMesh.edge_handle(secondUnrolledHalfEdge).idx()] = 1

    if mesh.calc_dihedral_angle(thirdHalfEdgeInFace) < 0:
        foldingDirection[unfoldedMesh.edge_handle(thirdUnrolledHalfEdge).idx()] = -1
    else:
        foldingDirection[unfoldedMesh.edge_handle(thirdUnrolledHalfEdge).idx()] = 1


    # Set glue numbers
    glueNumber[unfoldedMesh.edge_handle(unrolledLastHalfEdge).idx()] = mesh.edge_handle(lastHalfEdgeInFace).idx()
    glueNumber[unfoldedMesh.edge_handle(secondUnrolledHalfEdge).idx()] = mesh.edge_handle(secondHalfEdgeInFace).idx()
    glueNumber[unfoldedMesh.edge_handle(thirdUnrolledHalfEdge).idx()] = mesh.edge_handle(
        thirdHalfEdgeInFace).idx()

    # Check the two other half edges
    if secondHalfEdgeInFace in halfEdgeTree and not secondHalfEdgeInFace in cutEdges:
        # Get the face
        neighbourFace = mesh.face_handle(mesh.opposite_halfedge_handle(secondHalfEdgeInFace))
        unfoldNeighborsInTree(neighbourFace, secondHalfEdgeInFace, secondUnrolledHalfEdge, unrolledFromVertex, halfEdgeTree, mesh,
                              unfoldedMesh, isFoldingEdge, connections, cutEdges, glueNumber, foldingDirection)
    if thirdHalfEdgeInFace in halfEdgeTree and not thirdHalfEdgeInFace in cutEdges:
        # Get the face
        neighbourFace = mesh.face_handle(mesh.opposite_halfedge_handle(thirdHalfEdgeInFace))
        unfoldNeighborsInTree(neighbourFace, thirdHalfEdgeInFace, thirdUnrolledHalfEdge, unrolledToVertex, halfEdgeTree, mesh,
                              unfoldedMesh, isFoldingEdge, connections, cutEdges, glueNumber, foldingDirection)


def unfoldSpanningTree(unfoldedMesh, mesh, startingTriangle, halfEdgeTree, isFoldingEdge, connections, cutEdges, glueNumber, foldingDirection):
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
    firstUnrolledVertex = unfoldedMesh.add_vertex(secondUnrolled)
    secondUnrolledVertex = unfoldedMesh.add_vertex(thirdUnrolled)
    thirdUnrolledVertex = unfoldedMesh.add_vertex(firstUnrolled)
    # Create the face
    f = unfoldedMesh.add_face(firstUnrolledVertex, secondUnrolledVertex, thirdUnrolledVertex)
    # om.write_mesh('unfolding' + str(startingTriangle.idx()) + ".off", unfoldedMesh)


    connections[f.idx()] = startingTriangle.idx()

    # Now check the neighbours
    # Find the unrolled half-edges around the face
    firstUnrolledHalfEdge = unfoldedMesh.opposite_halfedge_handle(
        unfoldedMesh.halfedge_handle(firstUnrolledVertex))  # unfoldedMesh.opposite_halfedge_handle(
    # secondUnrolledHalfEdge = unfoldedMesh.opposite_halfedge_handle(
    #     unfoldedMesh.next_halfedge_handle(firstUnrolledHalfEdge))
    # thirdUnrolledHalfEdge = unfoldedMesh.opposite_halfedge_handle(
    #     unfoldedMesh.next_halfedge_handle(secondUnrolledHalfEdge))
    secondUnrolledHalfEdge = unfoldedMesh.next_halfedge_handle(firstUnrolledHalfEdge)
    thirdUnrolledHalfEdge = unfoldedMesh.next_halfedge_handle(secondUnrolledHalfEdge)

    if mesh.calc_dihedral_angle(firstHalfEdge) < 0:
        foldingDirection[unfoldedMesh.edge_handle(firstUnrolledHalfEdge).idx()] = -1
    else:
        foldingDirection[unfoldedMesh.edge_handle(firstUnrolledHalfEdge).idx()] = 1

    if mesh.calc_dihedral_angle(secondHalfEdge) < 0:
        foldingDirection[unfoldedMesh.edge_handle(secondUnrolledHalfEdge).idx()] = -1
    else:
        foldingDirection[unfoldedMesh.edge_handle(secondUnrolledHalfEdge).idx()] = 1

    if mesh.calc_dihedral_angle(thirdHalfEdge) < 0:
        foldingDirection[unfoldedMesh.edge_handle(thirdUnrolledHalfEdge).idx()] = -1
    else:
        foldingDirection[unfoldedMesh.edge_handle(thirdUnrolledHalfEdge).idx()] = 1

    # isFoldingEdge[unfoldedMesh.edge_handle(firstUnrolledHalfEdge).idx()] = True
    # print(unfoldedMesh.point(unfoldedMesh.from_vertex_handle(firstUnrolledHalfEdge)), unfoldedMesh.point(unfoldedMesh.to_vertex_handle(firstUnrolledHalfEdge)))
    # print(unfoldedMesh.point(unfoldedMesh.from_vertex_handle(secondUnrolledHalfEdge)),
    #       unfoldedMesh.point(unfoldedMesh.to_vertex_handle(secondUnrolledHalfEdge)))
    # print(unfoldedMesh.point(unfoldedMesh.from_vertex_handle(thirdUnrolledHalfEdge)),
    #       unfoldedMesh.point(unfoldedMesh.to_vertex_handle(thirdUnrolledHalfEdge)))

    # print(unfoldedMesh.edge_handle(firstUnrolledHalfEdge).idx(), unfoldedMesh.edge_handle(secondUnrolledHalfEdge).idx(), unfoldedMesh.edge_handle(thirdUnrolledHalfEdge).idx())
    # print(mesh.edge_handle(firstHalfEdge).idx(),mesh.edge_handle(secondHalfEdge).idx(),mesh.edge_handle(thirdHalfEdge).idx())
    # isFoldingEdge[unfoldedMesh.edge_handle(secondUnrolledHalfEdge).idx()] = True
    # isFoldingEdge[unfoldedMesh.edge_handle(thirdUnrolledHalfEdge).idx()] = True

    # Set glue numbers
    glueNumber[unfoldedMesh.edge_handle(firstUnrolledHalfEdge).idx()] = mesh.edge_handle(firstHalfEdge).idx()
    glueNumber[unfoldedMesh.edge_handle(secondUnrolledHalfEdge).idx()] = mesh.edge_handle(secondHalfEdge).idx()
    glueNumber[unfoldedMesh.edge_handle(thirdUnrolledHalfEdge).idx()] = mesh.edge_handle(thirdHalfEdge).idx()

    if firstHalfEdge in halfEdgeTree and not firstHalfEdge in cutEdges:
        # Get the neighbouring face
        # print("first")
        neighbourFace = mesh.face_handle(mesh.opposite_halfedge_handle(firstHalfEdge))
        unfoldNeighborsInTree(neighbourFace, firstHalfEdge, firstUnrolledHalfEdge, secondUnrolledVertex, halfEdgeTree, mesh,
                              unfoldedMesh, isFoldingEdge, connections, cutEdges, glueNumber, foldingDirection)

    if secondHalfEdge in halfEdgeTree and not secondHalfEdge in cutEdges:
        # print("sec")
        neighbourFace = mesh.face_handle(mesh.opposite_halfedge_handle(secondHalfEdge))
        unfoldNeighborsInTree(neighbourFace, secondHalfEdge, secondUnrolledHalfEdge, thirdUnrolledVertex, halfEdgeTree, mesh,
                              unfoldedMesh, isFoldingEdge, connections, cutEdges, glueNumber, foldingDirection)

    if thirdHalfEdge in halfEdgeTree and not thirdHalfEdge in cutEdges:
        # print("third")
        neighbourFace = mesh.face_handle(mesh.opposite_halfedge_handle(thirdHalfEdge))
        unfoldNeighborsInTree(neighbourFace, thirdHalfEdge, thirdUnrolledHalfEdge, firstUnrolledVertex, halfEdgeTree, mesh,
                              unfoldedMesh, isFoldingEdge, connections, cutEdges, glueNumber, foldingDirection)

    return None


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

def writeToBoth(string, file1, file2):
    file1.write(string)
    file2.write(string)

def writeToAll(string, file1, file2, file3):
    file1.write(string)
    file2.write(string)
    file3.write(string)

def writeSVG(filename, mesh, isFoldingEdge, isIntersected, glueNumber, foldingDirection, size, printNumbers):
    # Get the bounding box
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

    if size <= 0:
        boxSize = np.maximum(np.abs(xmax - xmin), np.abs(ymax - ymin))
    else:
        boxSize = size
    # xmin = np.minimum(xmin,ymin)
    # ymin = np.minimum(xmin,ymin)

    strokewidth = 0.002 * boxSize
    dashLength = 0.05 * boxSize
    spaceLength = 0.002 * boxSize

    textDistance = 0.02 * boxSize
    textStrokewidth = 0.05 * strokewidth
    textLength = 0.001 * boxSize
    fontsize = 0.015 * boxSize

    # dashLength = 0.2
    # spaceLength = 0.1
    frame = 0.05 * boxSize

    file = open(filename, 'w')
    fileCutting = open("cut" + filename, 'w')
    fileNumbers = open("numbers" + filename, 'w')

    writeToAll("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n", file, fileCutting, fileNumbers)
    writeToAll("<svg width=\"30.5cm\" height=\"30.5cm\" viewBox = \"" + str(xmin - frame) + " " + str(
        ymin - frame) + " " + str(boxSize + 2 * frame) + " " + str(
        boxSize + 2 * frame) + "\" version = \"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n", file, fileCutting, fileNumbers)

    # file.write("<rect x=\"" + str(xmin) + "\" y=\"" +str(ymin) + "\" width=\"100%\" height=\"100%\"/>\n")

    writeToAll("<g stroke = \"black\" stroke-width = \"" + str(strokewidth) + "\" >\n", file, fileCutting, fileNumbers)
    for edge in mesh.edges():
        # Get the two points
        he = mesh.halfedge_handle(edge, 0)
        vertex0 = mesh.point(mesh.from_vertex_handle(he))
        vertex1 = mesh.point(mesh.to_vertex_handle(he))

        # file.write("<line x1=\"" + str(vertex0[0]) + "\" y1=\"" + str(vertex0[1]) + "\" x2 = \"" + str(
        #    vertex1[0]) + "\" y2 = \"" + str(vertex1[1]) + "\"")

        writeToBoth(
            "<path d =\"M " + format(vertex0[0], '.10f') + "," + format(vertex0[1], '.10f') + " " + format(vertex1[0],
                                                                                                           '.10f') + "," + format(
                vertex1[1], '.10f') + "\" style=\"fill:none;stroke:", file, fileCutting)

        if isIntersected[edge.idx()]:
            writeToBoth("#ff0000", file, fileCutting)
        elif foldingDirection[edge.idx()] > 0:
            writeToBoth("#ff0000", file, fileCutting)
            #file.write("#00e64d")
        elif foldingDirection[edge.idx()] < 0:
            writeToBoth("#0066ff", file, fileCutting)
        else:
            writeToBoth("#000000", file, fileCutting)

        writeToBoth(";stroke-width:" + format(strokewidth,
                                             '.10f') + ";stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:4;stroke-dasharray:", file, fileCutting)

        if isFoldingEdge[edge.idx()]:
            writeToBoth(format(dashLength, '.10f') + ", " + format(spaceLength, '.10f'), file, fileCutting)
        else:
            writeToBoth("none", file, fileCutting)
            # file.write(" stroke-dasharray=\"" + str(dashLength) + " " + str(spaceLength) + "\"")
            # file.write(" stroke-dasharray=\"none\"")
            # file.write("style=\"stroke-dasharray:%.2f" %dashLength + ",%.2f" % spaceLength + ";stroke-dashoffset:0\"")
            # file.write(" style=\"fill:none;stroke:#000000;stroke-width:0.31999999;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:1;stroke-dasharray:1.28, 2.56;stroke-dashoffset:0;stroke-opacity:1\"")
        # else:
        #    file.write(" stroke-dasharray=\"none\"")
        # if isIntersected[edge.idx()]:
        #    file.write(" stroke=\"red\" ")
        writeToBoth(";stroke-dashoffset:0;stroke-opacity:1", file, fileCutting)
        writeToBoth("\" />\n", file, fileCutting)
        # write number if it is a cut edge
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

            if(printNumbers):
                writeToBoth("<text x=\"" + str(position[0]) + "\" y=\"" + str(position[1]) + "\" font-size=\"" + str(
                fontsize) + "\" stroke-width=\"" + str(textStrokewidth) + "\" transform=\"rotate(" + str(
                rotation) + "," + str(position[0]) + "," + str(position[1]) + ")\">" + str(
                glueNumber[edge.idx()]) + "</text>\n", file, fileNumbers)

    writeToAll("</g>\n", file, fileNumbers, fileCutting)
    writeToAll("</svg>", file, fileNumbers, fileCutting)
    file.close()


def line_intersect2(v1, v2, v3, v4, epsilon):
    '''
    judge if line (v1,v2) intersects with line(v3,v4)
    '''
    d = (v4[1] - v3[1]) * (v2[0] - v1[0]) - (v4[0] - v3[0]) * (v2[1] - v1[1])
    u = (v4[0] - v3[0]) * (v1[1] - v3[1]) - (v4[1] - v3[1]) * (v1[0] - v3[0])
    v = (v2[0] - v1[0]) * (v1[1] - v3[1]) - (v2[1] - v1[1]) * (v1[0] - v3[0])
    if d < 0:
        u, v, d = -u, -v, -d
    return ((0 + epsilon) <= u <= (d - epsilon)) and ((0 + epsilon) <= v <= (d - epsilon))


def point_in_triangle2(A, B, C, P, epsilon):
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


def tri_intersect2(t1, t2, epsilon):
    '''
    judge if two triangles in a plane intersect
    '''
    if line_intersect2(t1[0], t1[1], t2[0], t2[1], epsilon): return True
    if line_intersect2(t1[0], t1[1], t2[0], t2[2], epsilon): return True
    if line_intersect2(t1[0], t1[1], t2[1], t2[2], epsilon): return True
    if line_intersect2(t1[0], t1[2], t2[0], t2[1], epsilon): return True
    if line_intersect2(t1[0], t1[2], t2[0], t2[2], epsilon): return True
    if line_intersect2(t1[0], t1[2], t2[1], t2[2], epsilon): return True
    if line_intersect2(t1[1], t1[2], t2[0], t2[1], epsilon): return True
    if line_intersect2(t1[1], t1[2], t2[0], t2[2], epsilon): return True
    if line_intersect2(t1[1], t1[2], t2[1], t2[2], epsilon): return True
    inTri = True
    inTri = inTri and point_in_triangle2(t1[0], t1[1], t1[2], t2[0], epsilon)
    inTri = inTri and point_in_triangle2(t1[0], t1[1], t1[2], t2[1], epsilon)
    inTri = inTri and point_in_triangle2(t1[0], t1[1], t1[2], t2[2], epsilon)
    if inTri == True: return True
    inTri = True
    inTri = inTri and point_in_triangle2(t2[0], t2[1], t2[2], t1[0], epsilon)
    inTri = inTri and point_in_triangle2(t2[0], t2[1], t2[2], t1[1], epsilon)
    inTri = inTri and point_in_triangle2(t2[0], t2[1], t2[2], t1[2], epsilon)
    if inTri == True: return True
    return False

def unfold(mesh, weightFactors):
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

    # phi = -np.pi
    # theta = 0.9*np.pi
    # cutVector = np.array([np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta), np.cos(theta)])
    # print(cutVector)
    cutVector = np.array([1.0, 1.0, 1.0])

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
            # np.abs(mesh.calc_dihedral_angle(edge))

        edgeVector = mesh.calc_edge_vector(edge)
        flatAngleWeights[edge.idx()] = np.abs(np.dot(cutVector, edgeVector)) / np.linalg.norm(edgeVector)

    weights = weightFactors[0] * minPerimeterEdgeWeights + \
              weightFactors[1] * flatAngleWeights + \
              weightFactors[2] * dihedralAngleWeights
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
    glueNumber = np.empty(numUnfoldedEdges, dtype=int)
    foldingDirections = np.empty(numUnfoldedEdges, dtype=int)

    # Face onnection array
    connections = np.empty(numFaces, dtype=int)

    # Get the points of the triangle
    startingTriangle = forest[0][0]

    unfoldSpanningTree(unfoldedMesh, mesh, startingTriangle, halfEdgeTree, isFoldingEdge, connections, [], glueNumber, foldingDirections)

    # Resolve the intersection

    # Find all intersections
    epsilon = 1E-12
    faceIntersections = []
    isInterSected = np.zeros(numUnfoldedEdges, dtype=bool)
    for face1 in unfoldedMesh.faces():
        for face2 in unfoldedMesh.faces():
            if face2.idx() < face1.idx():
                # Get the triangle faces
                triangle1 = []
                triangle2 = []
                for halfedge in unfoldedMesh.fh(face1):
                    triangle1.append(unfoldedMesh.point(unfoldedMesh.from_vertex_handle(halfedge)))
                for halfedge in unfoldedMesh.fh(face2):
                    triangle2.append(unfoldedMesh.point(unfoldedMesh.from_vertex_handle(halfedge)))
                if tri_intersect2(triangle1, triangle2, epsilon):
                    # print("Intersection: " + str(face1.idx()) + " and " + str(face2.idx()))
                    faceIntersections.append([connections[face1.idx()], connections[face2.idx()]])
                    for edge in unfoldedMesh.fe(face1):
                        isInterSected[edge.idx()] = True
                    for edge in unfoldedMesh.fe(face2):
                        isInterSected[edge.idx()] = True

    print("Wir haben " + str(len(faceIntersections)) + " Überschneidungen.")

    # print("Minimum: " + str(minIntersections) + " for i=" +str(minIntersectionIndex))

    # Find the paths
    # Make the tree a real graph
    spanningGraph = nx.Graph()
    for edge in spanningTree:
        face1 = mesh.face_handle(mesh.halfedge_handle(edge, 0))
        face2 = mesh.face_handle(mesh.halfedge_handle(edge, 1))
        spanningGraph.add_edge(face1.idx(), face2.idx())

    print("Berechne Pfade")
    paths = []
    for intersection in faceIntersections:
        # Find the path
        paths.append(nx.algorithms.shortest_paths.shortest_path(spanningGraph, intersection[0], intersection[1]))
        # print(intersection[0])

    # Save paths wrt edges
    edgepaths = []
    for path in paths:
        edgepath = []
        for i in range(len(path) - 1):
            # Find edge between ith and next face
            for he in mesh.fh(mesh.face_handle(path[i])):
                if mesh.face_handle(mesh.opposite_halfedge_handle(he)) == mesh.face_handle(path[i + 1]):
                    edgepath.append(mesh.edge_handle(he).idx())
        edgepaths.append(edgepath)

    allEdgesInPaths = list(set().union(*edgepaths))
    # See how often all Edges arise in the paths
    numEdgesInPaths = []
    for edge in allEdgesInPaths:
        num = 0
        for path in edgepaths:
            if edge in path:
                num = num + 1
        numEdgesInPaths.append(num)

    S = []
    C = []
    gamma = 1
    cutWeights = (1 - gamma) * (1 - weights) + gamma
    print("Set covering algorithm")
    while len(C) != len(paths):
        # Make the new weights of the spanningTree
        spanningWeights = np.empty(len(allEdgesInPaths))
        for i in range(len(allEdgesInPaths)):
            # Check how often the ith edge is in the path in C
            currentEdge = allEdgesInPaths[i]
            numInC = 0
            for path in C:
                if currentEdge in path:
                    numInC = numInC + 1
            if (numEdgesInPaths[i] - numInC) > 0:
                spanningWeights[i] = cutWeights[edge] / (numEdgesInPaths[i] - numInC)
            else:
                spanningWeights[i] = 1000
        # Sort the edges in the paths
        minimal = np.argmin(spanningWeights)
        S.append(allEdgesInPaths[minimal])
        # Find all paths that the minimal edge is in and append to C
        for path in edgepaths:
            if allEdgesInPaths[minimal] in path and not path in C:
                C.append(path)
    print("Wir müssen " + str(len(S)) + " Kanten schneiden um Überlappungen zu vermeiden.")

    # Make halfedge list of cut edges
    cutHalfEdges = []
    for ind in S:
        cutHalfEdges.append(mesh.halfedge_handle(mesh.edge_handle(ind), 0))
        cutHalfEdges.append(mesh.halfedge_handle(mesh.edge_handle(ind), 1))

    print("Abwickeln")
    # Find one face in each connected component
    # Make forest
    components = nx.Graph()
    for edge in spanningTree:
        if not edge.idx() in S:
            face1 = mesh.face_handle(mesh.halfedge_handle(edge, 0))
            face2 = mesh.face_handle(mesh.halfedge_handle(edge, 1))
            components.add_edge(face1.idx(), face2.idx())

    #Find single triangles
    for face in mesh.faces():
        allCut = True
        for edge in mesh.fe(face):
            if not edge.idx() in S:
                components.add_edge(face.idx(), face.idx())



    connectedComponents = nx.algorithms.components.connected_components(components)

    unfoldedComponents = []
    for c in connectedComponents:
        startingTriangleInd = list(c)[0]
        startingTriangle = mesh.face_handle(startingTriangleInd)
        simpleUnfolded = om.TriMesh()
        simpleIsFoldingEdge = np.zeros(numUnfoldedEdges, dtype=bool)
        simpleconnections = np.empty(numFaces, dtype=int)
        simpleGlueNumber = np.empty(numUnfoldedEdges, dtype=int)
        simpleFoldingDirection = np.empty(numUnfoldedEdges, dtype=int)

        unfoldSpanningTree(simpleUnfolded, mesh, startingTriangle, halfEdgeTree, simpleIsFoldingEdge, simpleconnections,
                           cutHalfEdges, simpleGlueNumber, simpleFoldingDirection)
        # writeSVG("unfolding" + str(count) + ".svg", simpleUnfolded, simpleIsFoldingEdge, np.zeros(numUnfoldedEdges,dtype=bool), simpleGlueNumber)
        unfoldedComponents.append([simpleUnfolded, simpleIsFoldingEdge, simpleGlueNumber, simpleFoldingDirection])

    return [unfoldedMesh, isFoldingEdge, isInterSected, glueNumber, foldingDirections], unfoldedComponents