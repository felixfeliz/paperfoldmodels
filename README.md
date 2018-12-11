# paperfoldmodels
A python module that unfolds triangulated surfaces to a two-dimensional net. The result can be used to create papercraft models. It is mostly based on the algorithm presented in this [report](https://geom.ivd.kit.edu/downloads/proj-paper-models_cut_out_sheets.pdf) of Straub and Prautzsch.

The python bindings of [openmesh](http://www.openmesh.org) are used and thus the original mesh can be loaded from any format supported by it.
The output consists of several SVG files, each containing an intersection-free component of the unfolding.

## Dependencies
numpy

openmesh

networkx

## Example
The usage is shown in testUnfold.py

Unfolding of an icosahedron:
![Icosahedron](icosahedron.svg)

## Method
The algorithm consists of three steps:

1. Find a minimum spanning tree of the dual graph of the mesh.
2. Unfold the dual graph.
3. Remove self-intersections by adding additional cuts along edges. 
