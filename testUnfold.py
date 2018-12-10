import numpy as np
import openmesh as om
from unfoldMesh import unfold, writeSVG
#import sys



def main():
    printNumbers = True

    #FILENAME = 'models/icosahedron.obj'
    #FILENAME = 'original.off'
    #FILENAME = 'reduced.obj'
    #FILENAME = 'models/reducedTeddy.obj'
    FILENAME = 'models/bunny.stl'
    #FILENAME = 'models/tree.obj'
    #FILENAME = 'models/polyhedron.obj'
    #FILENAME = 'models/kndC.obj'

    # for opt in sys.argv[1:]:
    #     print(opt)
    #     if opt == "-n":
    #         printNumbers = True

    # Import the mode
    mesh = om.read_trimesh(FILENAME)

    fullUnfolded, unfoldedComponents = unfold(mesh,[3.0,1.0,0.0])


    #Compute maxSize of the components
    maxSize = 0.0
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
    writeSVG('unfolding.svg', fullUnfolded[0], fullUnfolded[1], fullUnfolded[2], fullUnfolded[3], -1, printNumbers)

    for i in range(len(unfoldedComponents)):
        writeSVG("unfolding" + str(i) + ".svg", unfoldedComponents[i][0], unfoldedComponents[i][1], np.zeros(fullUnfolded[0].n_edges(), dtype=bool), unfoldedComponents[i][2], maxSize, printNumbers)

    print("Wir haben " + str(len(unfoldedComponents)) + " Komponenten geschrieben.")

if __name__=='__main__':
    main()


