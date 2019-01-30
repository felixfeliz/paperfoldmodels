import numpy as np
import openmesh as om
from unfoldmesh import unfold, writeSVG
from argparse import ArgumentParser



def main():

    #FILENAME = 'models/icosahedron.obj'
    #FILENAME = 'original.off'
    #FILENAME = 'reduced.obj'
    #FILENAME = 'models/reducedTeddy.obj'
    #FILENAME = 'models/bunny.stl'
    #FILENAME = 'models/tree.obj'
    #FILENAME = 'models/polyhedron.obj'
    #FILENAME = 'models/kndC.obj'

    parser = ArgumentParser()
    parser.add_argument("-f", "--file", dest="filename", default="models/icosahedron.obj",
                        help="path to the model", metavar="FILE")
    parser.add_argument("-n", "--numbers",
                        action="store_true", dest="printNumbers", default=False,
                        help="Print numbers on the cut edges")
    args = parser.parse_args()


    printNumbers = args.printNumbers

    # Import the mode
    mesh = om.read_trimesh(args.filename)

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
    writeSVG('unfolding.svg', fullUnfolded[0], fullUnfolded[1], fullUnfolded[2], fullUnfolded[3], fullUnfolded[4], -1, printNumbers)

    for i in range(len(unfoldedComponents)):
        writeSVG("unfolding" + str(i) + ".svg", unfoldedComponents[i][0], unfoldedComponents[i][1], np.zeros(fullUnfolded[0].n_edges(), dtype=bool), unfoldedComponents[i][2], unfoldedComponents[i][3], maxSize, printNumbers)

    print("Wir haben " + str(len(unfoldedComponents)) + " Komponenten geschrieben.")

if __name__=='__main__':
    main()


