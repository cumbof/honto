#!/usr/bin/env python

__authors__ = ("Paolo Giulio Franciosa (paolo.franciosa@uniroma1.it)"
               "Nicola Apollonio (nicola.apollonio@cnr.it)"
               "Daniele Santoni (daniele.santoni@iasi.cnr.it)"
               "Fabio Cumbo (fabio.cumbo@gmail.com)")

__version__ = "0.1.0"
__date__ = "May 5, 2022"

import os

def inputGraphFromFile(file_colors, file_network, graph, colors, colorGroups, weight_threshold, isolated):
    """
    Read a graph from files

    :param file_colors:         Path to the file with a mapping between node names and color.
                                Default color value is "X"
    :param graph:               NetworkX graph, initially empty
    :param colors:              Set of colors in graph, initially empty
    :param colorGroups:         Set of color groups in graph, initially empty
    :param file_network:        Path to the file with couples of edges, one for each line. The representing graph is not oriented.
                                Nodes could not be defined in "file_colors". In this case, nodes are considered with color "X".
                                There could be a weight for each edge. In this case, edges are considered if their weight is greater 
                                then or equals to "weight_threshold"
    :param weight_threshold:    Threshold on edges weight. Edges are not considered if their weight is lower than this threshold.
    :param isolated:            Also insert isolated nodes if true
    """
    
    # Add edgens to the NetworkX graph
    addEdges(graph, file_network, weight_threshold)
    # Read colors from input file
    readColors(graph, file_colors, isolated)
    # Populate colors
    populateColors(graph, colors, colorGroups)

def readColors(graph, filepath, isolated):
    """
    Read nodes with colors from file
    Set node attribute as color and color group
    Nodes with names "xxx_y" are allowed
    Extensions are trimmed out from edges while reading input file
    If a node has an extension, verify whether it already exists with a color:
        - assign the color defined in the input file if color is empty or equals to the color defined in the input file
        - assign color "X" if color is different from the one defined in the input file
    
    :param graph:       NetworkX graph
    :param filepath:    Input file with colors
    :param isolated:    Also insert isolated nodes if true
    """

    if filepath and os.path.exists(filepath):
        with open(filepath) as f_in:
            for line in f_in:
                line = line.strip()
                if line:
                    fields = line.split()
                    node = fields[0][:fields[0].index("_")].upper() if "_" in fields[0] else fields[0].upper()
                    if len(fields) == 1:
                        color = "X"
                        colorGroup = "X"
                    else:
                        color = fields[1]
                        if color in "RS":
                            color = "X"
                        if color in "JAKLB":
                            colorGroup = "JAKLB"
                        elif color in "DYVTMNZWUO":
                            colorGroup = "DYVTMNZWUO"
                        elif color in "CGEFHIPQ":
                            colorGroup = "CGEFHIPQ"
                        elif color in "RSX":
                            colorGroup = "X"  
                        else:
                            colorGroup = fields[1]
                    
                    if node not in graph.nodes() and isolated:
                        graph.add_node(node, color="", colorGroup="colGroup")
                    
                    if node in graph.nodes():
                        oldColor = graph.nodes[node]["color"]
                        if oldColor and oldColor != colore:
                            color = "X"
                            colorGroup = "X"  
                        graph.add_node(node, color=color)
                        graph.add_node(node, colorGroup=colorGroup)

    # Set color to "X" for all nodes in graph with no color
    for node in graph.nodes:
        if graph.nodes[node]["color"] == "":
            graph.add_node(node, color="X", colorGroup="X")

def populateColors(graph, colors, colorGroups):
    """
    Iterate over nodes and add color and group of colors to "colors" and "colorGroups"

    :param graph:           NetworkX graph
    :param colors:          List of colors
    :param colorGroups:     List of color groups
    """

    for node in graph.nodes:
        colors.add(graph.nodes[node]["color"])
        colorGroups.add(graph.nodes[node]["colorGroup"])

def addEdges(graph, filepath, weight_threshold):
    """
    Read edges from file and add them to the NetworkX graph (it also adds nodes).
    The edge is added in both the adjacency lists if not already present.
    A simmetric couple of edges can be defined in the input file for each couple of adjacent nodes
    
    :param graph:               NetworkX graph
    :param filepath:            Input file
    :param weight_threshold:    Threshold on edges weight. Edges are not considered if their weight is lower than this threshold
    """

    with open(filepath) as f_in:
        for line in f_in:
            line = line.strip()
            if line:
                fields = line.split()
                node1 = fields[0][:fields[0].index("_")].upper() if "_" in fields[0] else fields[0].upper()
                node2 = fields[1][:fields[1].index("_")].upper() if "_" in fields[1] else fields[1].upper()
                if len(fields) == 3:
                    weight = int(fields[2])
                    if weight >= weight_threshold:
                        graph.add_node(node1, color="", colorGroup="")
                        graph.add_node(node2, color="", colorGroup="")
                        graph.add_edge(node1, node2)
                        graph.add_edge(node2, node1)
                else:
                    graph.add_node(node1, color="", colorGroup="")
                    graph.add_node(node2, color="", colorGroup="")
                    graph.add_edge(node1, node2)
                    graph.add_edge(node2, node1)

def dumpGraph(graph, outnet_filepath, outcolors_filepath):
    """
    Dump a NetworkX graph to file

    :param graph:               NetworkX graph
    :param outnet_filepath:     Path to the output file with the network
    :param outcolors_filepath:  Path to the output file with colos
    """

    # Dump network
    with open(outnet_filepath, "w+") as f_net:
        for edge in graph.edges():
            f_net.write("{}\t{}\n".format(edge[0], edge[1]))

    # Dump node colors
    with open(outcolors_filepath, "w+") as f_colors:
        for node, dolor in graph.nodes.data("color"):
            f_colors.write("{}\t{}\n".format(node, color))
