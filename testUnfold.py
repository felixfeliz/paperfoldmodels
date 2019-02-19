import openmesh as om
from unfoldmesh import unfold, writeSVG, findBoundingBox
from argparse import ArgumentParser



def main():
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

    fullUnfolded, unfoldedComponents = unfold(mesh)


    #Compute maxSize of the components
    # Alle Komponenten müssen auf die selbe Größe skaliert werden, die vom größten Bauteil
    maxSize = 0
    for unfolding in unfoldedComponents:
        [xmin, ymin, boxSize] = findBoundingBox(unfolding[0])
        if boxSize > maxSize:
            maxSize = boxSize


    #Write SVG
    if printNumbers:
        basefilename = "unfoldingNumbers"
    else:
        basefilename = "unfolding"

    for i in range(len(unfoldedComponents)):
        writeSVG(basefilename + str(i) + ".svg", unfoldedComponents[i], maxSize, printNumbers)

    print("We wrote " + str(len(unfoldedComponents)) + " connected components.")

if __name__=='__main__':
    main()


