#!/usr/bin/env python

__authors__ = ("Paolo Giulio Franciosa (paolo.franciosa@uniroma1.it)",
               "Nicola Apollonio (nicola.apollonio@cnr.it)",
               "Daniele Santoni (daniele.santoni@iasi.cnr.it)",
               "Fabio Cumbo (fabio.cumbo@gmail.com)")

__version__ = "0.1.0"
__date__ = "May 10, 2022"

import sys

# Define tool name
TOOL_ID="honto"

# Control current Python version
# It requires Python 3 or higher
if sys.version_info[0] < 3:
    raise Exception("{} requires Python 3, your current Python version is {}.{}.{}"
                    .format(TOOL_ID, sys.version_info[0], sys.version_info[1], sys.version_info[2]))

import os, math, time, tqdm
import honto.utils as utils
import argparse as ap
import pandas as pd
import networkx as nx
import seaborn as sb
import multiprocessing as mp
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm, LogNorm, Normalize
from functools import partial

def read_params():
    p = ap.ArgumentParser( description = ( "A novel method for assessing and measuring homophily in networks" ),
                           formatter_class = ap.ArgumentDefaultsHelpFormatter )
    # Inputs
    p.add_argument( '--input_edges', 
                    type = str,
                    help = "Path to the file with the list of edges" )
    p.add_argument( '--input_nodes', 
                    type = str,
                    help = "Path to the file with nodes and colors" )
    p.add_argument( '--weight_threshold', 
                    type = int,
                    default = 700,
                    help = "Threshold for considering edges based in their weight" )
    p.add_argument( '--isolated',
                    action = 'store_true',
                    default = False,
                    help = "Insert isolated nodes" )
    # Rescaling z-scores arguments
    p.add_argument( '--log_transform',
                    action = 'store_true',
                    default = False,
                    help = "Log-transform z-scores" )
    p.add_argument( '--scale_factor',
                    type = float,
                    help = "Rescale z-scores with this constant before log-transforming values" )
    p.add_argument( '--scale_from_one',
                    action = 'store_true',
                    default = False,
                    help = "Set z-scores to 1 if lower than 1 before log-transforming values" )
    # Matplotlib arguments
    p.add_argument( '--cmap', 
                    type = str,
                    default = "PiYG",
                    help = "Heatmap colormap" )
    p.add_argument( '--vmin', 
                    type = float,
                    default = 2.2957,
                    help = "Min value to anchor the colormap" )
    p.add_argument( '--vmax', 
                    type = float,
                    default = 4.3957,
                    help = "Max value to anchor the colormap" )
    p.add_argument( '--center', 
                    type = float,
                    default = 2.9957,
                    help = "The value at which to center the colormap when plotting divergant data" )
    p.add_argument( '--cbar', 
                    action = 'store_true',
                    default = False,
                    help = "Whether to draw a colorbar" )
    # General purpose arguments
    p.add_argument( '--nproc', 
                    type = int,
                    default = 1,
                    help = "Make the computation of the z-scores parallel for singletons" )
    p.add_argument( '--overwrite',
                    action = 'store_true',
                    default = False,
                    help = "Overwrite results if already exist" )
    p.add_argument( '--verbose',
                    action = 'store_true',
                    default = False,
                    help = "Print results in real time" )
    p.add_argument( '-v', 
                    '--version', 
                    action = 'version',
                    version = '{} v{} ({})'.format( TOOL_ID, __version__, __date__ ),
                    help = "Print current {} version and exit".format(TOOL_ID) )
    return p.parse_args()

def plot_heatmap(input_edges, z_score_edges, cmap="PiYG", vmin=None, vmax=None, center=None, cbar=False):
    scores = pd.DataFrame(z_score_edges).T.fillna(0)
    fig, ax = plt.subplots(figsize=(5,5))
    sb.heatmap(scores, square=True, cmap=cmap, vmin=vmin, vmax=vmax, center=center, cbar=cbar, ax=ax)
    plt.tick_params(axis="both", which="major", labelsize=12, labelleft=True, labeltop=True, 
                    left=False, right=False, labelbottom=False, bottom=False, top=False)
    plt.savefig("{}.pdf".format(os.path.splitext(input_edges)[0]))

def transform_z_score(zscores, scale_value=None, from_one=False):
    """
    Rescale z-scores

    :param zscores:         Z-scores
    :param scale_value:     Rescale z-scare with this value: z-score - scale_value + 1
    :param from_one:        Set z-scores to 1 if lower than 1
    """

    if scale_value is None and not from_one:
        auto_scale = min([zscores[color1][color2] for color1 in zscores for color2 in zscores[color1]])

    for color1 in zscores:
        for color2 in zscores[color1]:
            value = zscores[color1][color2]
            if scale_value is not None:
                value += scale_value
            elif from_one:
                value = 1 if value < 1 else value
            else:
                value += 1 - auto_scale
            zscores[color1][color2] = math.log(value)
    return zscores

def compute_print_z_score_singletons(input_edges, n_nodes, memo_falling_frac, 
                                     supergraph_degrees, supergraph_degree_histo, supergraph, 
                                     col_list, excluded_colors, nproc=1, verbose=False):
    pairs_p_3 = None
    with open("{}_zscores_singletons.txt".format(os.path.splitext(input_edges)[0]), "w+") as f:
        f.write("# {}\n".format(os.path.basename(input_edges)))
        f.write("# color\tn\tm\tz_score\n")
        for color in col_list:
            if color not in excluded_colors:
                induced_subgraph = induced_color_subgraph(supergraph, "color", color)
                n_sub, m_sub, z_score_singletons, pairs_p_3 = subgraph_profile(induced_subgraph, n_nodes, memo_falling_frac, 
                                                                               supergraph_degrees, supergraph_degree_histo, supergraph, pairs_p_3, 
                                                                               nproc=nproc, verbose=verbose)
                f.write("{}\t{}\t{}\t{:.4f}\n".format(color, n_sub, m_sub, z_score_singletons))
    
def subgraph_profile(graph, n_nodes, memo_falling_frac, 
                     supergraph_degrees, supergraph_degree_histo, supergraph, pairs_p_3, nproc=1, verbose=False):
    """
    Compute singleton for subgraph z_scores
    """

    n_sub = graph.number_of_nodes()
    m_sub = graph.number_of_edges()
    z_score_singletons, pairs_p_3 = z_score_number_singletons(graph, n_nodes, memo_falling_frac, 
                                                              supergraph_degrees, supergraph_degree_histo, supergraph, pairs_p_3, nproc=nproc, verbose=verbose)
    
    return n_sub, m_sub, z_score_singletons, pairs_p_3

def print_z_score_edges(input_edges, col_list, z_score_edges):
    with open("{}_zscores_edges.txt".format(os.path.splitext(input_edges)[0]), "w+") as f:
        f.write("# {}\n".format(os.path.basename(input_edges)))
        f.write("# color1\tcolor2\tz_score\n")
        for color1 in col_list:
            for color2 in col_list:
                f.write("{}\t{}\t{:.4f}\n".format(color1, color2, z_score_edges[color1][color2]))

def z_score_number_edges(memo_falling_frac, n_edges, n_nodes, col_list, 
                         node_distribution, memo_falling_power, edge_distribution, pi_3):
    """
    Compute the z_score of the number of bichromatic and monochromatic edges
    """

    expected_edges = expected_number_edges(memo_falling_frac, n_edges, n_nodes, col_list, node_distribution, memo_falling_power)
    variance_edges = variance_number_edges(n_edges, n_nodes, memo_falling_power, node_distribution, col_list, pi_3)

    z_scores = dict({color: dict() for color in col_list})
    min_value = None
    for color1 in col_list:
        for color2 in col_list:
            sigma = math.sqrt(variance_edges[color1][color2])
            if sigma == 0:
                z_scores[color1][color2] = 0
            else:
                z_scores[color1][color2] = (edge_distribution[color1][color2] - expected_edges[color1][color2]) / sigma
            
            if min_value is None:
                min_value = z_scores[color1][color2]
            else:
                if z_scores[color1][color2] < min_value:
                    min_value = z_scores[color1][color2]
        
    return z_scores, min_value

def variance_number_edges(n_edges, n_nodes, memo_falling_power, node_distribution, col_list, pi_3):
    """
    Compute the expected variance of the number of bichromatic and monochromatic edges
    """

    n2 = falling_power(n_nodes, 2, memo_falling_power)
    n3 = falling_power(n_nodes, 3, memo_falling_power)
    n4 = falling_power(n_nodes, 4, memo_falling_power)

    variance = dict({color: dict() for color in col_list})
    for i, color1 in enumerate(col_list):
        for j in range(i, len(col_list)):
            color2 = col_list[j]
            if color1 == color2:
                a = node_distribution[color1]
                a2 = falling_power(a, 2, memo_falling_power)
                a3 = falling_power(a, 3, memo_falling_power)
                a4 = falling_power(a, 4, memo_falling_power)

                exp_m = n_edges * a2 / n2
                variance[color1][color1] = exp_m - exp_m*exp_m + a4*(n_edges*n_edges-n_edges)/n4 + 2*(n_nodes-a)*a3*pi_3/n4
            else:
                a = node_distribution[color1]
                b = node_distribution[color2]

                exp_m = 2 * n_edges * a * b / n2
                term = 4 * falling_power(a, 2, memo_falling_power) * falling_power(b, 2, memo_falling_power) / n4
                variance[color1][color2] = exp_m - exp_m*exp_m + 2 * (((a * falling_power(b, 2, memo_falling_power) + falling_power(a, 2, memo_falling_power) * b) / n3 - term) * pi_3 + term * n_edges * (n_edges-1) / 2)
                variance[color2][color1] = variance[color1][color2]
    
    return variance

def expected_number_edges(memo_falling_frac, n_edges, n_nodes, col_list, node_distribution, memo_falling_power):
    """
    Compute the expected value of bichromatic and monochromatic edges
    """

    expected = dict({color: dict() for color in col_list})
    for i, color1 in enumerate(col_list):
        for j in range(i, len(col_list)):
            color2 = col_list[j]
            if color1 == color2:
                a = node_distribution[color1]
                expected[color1][color1] = n_edges * falling_frac(a, n_nodes, 2, memo_falling_frac)
            else:
                expected[color1][color2] = 2 * n_edges * node_distribution[color1] * node_distribution[color2] / falling_power(n_nodes, 2, memo_falling_power)
                expected[color2][color1] = expected[color1][color2]
    
    return expected

def z_score_number_singletons(graph, n_nodes, memo_falling_frac, 
                              supergraph_degrees, supergraph_degree_histo, supergraph, pairs_p_3, nproc=1, verbose=False):
    """
    Compute the z_score of the number of isolated nodes
    """

    v_fast, pairs_p_3 = variance_number_singletons_fast(graph, n_nodes, memo_falling_frac, 
                                                        supergraph_degrees, supergraph_degree_histo, supergraph, pairs_p_3, nproc=nproc, verbose=verbose)
    sigma = math.sqrt(v_fast)
    if sigma == 0:
        return 0, pairs_p_3

    z = (nx.degree_histogram(graph)[0] - expected_number_singletons(graph, n_nodes, supergraph_degrees, memo_falling_frac)) / sigma

    return z, pairs_p_3

def variance_number_singletons_fast(graph, n_nodes, memo_falling_frac, 
                                    supergraph_degrees, supergraph_degree_histo, supergraph, pairs_p_3, nproc=1, verbose=False):
    graph_nodes = graph.number_of_nodes()

    summation = 0
    summation += all_pairs_contribution_rough(n_nodes, graph_nodes, memo_falling_frac, supergraph_degree_histo)
    summation -= edge_contribution_rough(n_nodes, graph_nodes, memo_falling_frac, supergraph, supergraph_degrees)
    summation -= vertices_contribution_rough(n_nodes, graph_nodes, memo_falling_frac, supergraph, supergraph_degrees)
    
    correction, pairs_p_3 = p_3_correction(n_nodes, graph_nodes, memo_falling_frac, 
                                           supergraph_degrees, supergraph, pairs_p_3, nproc=nproc, verbose=verbose)
    summation += correction
    
    expected = expected_number_singletons(graph, n_nodes, supergraph_degrees, memo_falling_frac)
    var = expected * (1 - expected) + falling_frac(graph_nodes, n_nodes, 2, memo_falling_frac) * summation

    return var, pairs_p_3

def p_3_correction(n, a, memo_falling_frac, supergraph_degrees, supergraph, pairs_p_3, nproc=1, verbose=False):
    correction = 0
    pairs_p_3 = all_p3_pairs(pairs_p_3, supergraph, nproc=nproc, verbose=verbose)
    for u, v in pairs_p_3:
        correction -= falling_frac(n-a, n-2, supergraph_degrees(u)+supergraph_degrees(v), memo_falling_frac)
        correction += falling_frac(n-a, n-2, neighbor_union_size(u, v, supergraph_degrees, supergraph), memo_falling_frac)
    
    return 2 * correction, pairs_p_3

def vertices_contribution_rough(n, a, memo_falling_frac, supergraph, supergraph_degrees):
    result = 0
    for node in supergraph.nodes():
        degree = supergraph_degrees(node)
        result += falling_frac(n-a, n-2, 2*degree, memo_falling_frac)

    return result

def edge_contribution_rough(n, a, memo_falling_frac, supergraph, supergraph_degrees):
    result = 0
    for edge in supergraph.edges():
        degree1 = supergraph_degrees[edge[0]]
        degree2 = supergraph_degrees[edge[1]]
        result += falling_frac(n-a, n-2, degree1+degree2, memo_falling_frac)
    
    return 2 * result

def all_pairs_contribution_rough(n, a, memo_falling_frac, supergraph_degree_histo):
    # Consider every couple as they have disjunct neighbors
    # Also consider (u, u)
    result = 0
    for degree1 in range(len(supergraph_degree_histo)):
        if supergraph_degree_histo[degree1] > 0:
            for degree2 in range(len(supergraph_degree_histo)):
                result += supergraph_degree_histo[degree1] * supergraph_degree_histo[degree2] * falling_frac(n-a, n-2, degree1+degree2, memo_falling_frac)
    
    return result

def all_p3_pairs_par(node, supergraph):
    pairs_p_3_local = set()
    
    for edge1 in supergraph.edges(node):
        for edge2 in supergraph.edges(node):
            u, v = edge1[1], edge2[1]
            if u >= v:
                continue
            if not supergraph.has_edge(u, v):
                # No effect if already present
                pairs_p_3_local.add((u, v))
    
    return pairs_p_3_local

def all_p3_pairs(pairs_p_3, supergraph, nproc=1, verbose=False):
    """
    Build the set of all edges (u,v) such that
        - u < v
        - u, v have a common neighbor
        - u, v are not adjacent
    """
    
    if pairs_p_3 is not None:
        return pairs_p_3
    
    pairs_p_3 = set()        

    nodes = list(nx.nodes(supergraph))
    if verbose:
        print("Processing nodes")

    all_p3_pairs_par_partial = partial(all_p3_pairs_par, supergraph=supergraph)

    with mp.Pool(processes=nproc) as pool:
        # Run jobs
        jobs = tqdm.tqdm(pool.imap(all_p3_pairs_par_partial, nodes), total=len(nodes), disable=(not verbose))
        pairs_p_3 = set().union(*jobs)
        
    return pairs_p_3

def neighbor_union_size(u, v, supergraph_degrees, supergraph):
    b = supergraph_degrees[u]
    for n in supergraph.neighbors(v):
        if not supergraph.has_edge(u, n):
            b += 1
    
    return b

def compute_falling_frac(num, den, b):
    """
    Return falling_power(num, b)/falling_power(den, b) avoiding generating too large intermediate values
    Assumes num <= den
    """

    if b > num:
        return 0

    p = 1
    while b > 0:
        p *= num / den
        num -= 1
        den -= 1
        b -= 1
    
    return p

def falling_frac(num, den, b, memo_falling_frac):
    if (num, den, b) not in memo_falling_frac:
        memo_falling_frac[(num, den, b)] = compute_falling_frac(num, den, b)
    
    return memo_falling_frac[(num, den, b)]

def compute_falling_power(n, k):
    if k > n:
        return 0

    prod = 1
    for x in range(n, n-k, -1):
        prod *= x
    
    return prod

def falling_power(x, y, memo_falling_power):
    if (x, y) not in memo_falling_power.keys():
        memo_falling_power[(x, y)] = compute_falling_power(x, y)
    
    return memo_falling_power[(x, y)]

def expected_number_singletons(graph, n_nodes, supergraph_degrees, memo_falling_frac):
    """
    Compute the expected value of isolated nodes in the induced subgraph
    """

    graph_nodes = graph.number_of_nodes()
    expected = 0
    for _, degree in supergraph_degrees:
        expected += falling_frac(n_nodes-graph_nodes, n_nodes-1, degree, memo_falling_frac)
    expected *= graph_nodes/n_nodes

    return expected

def induced_color_subgraph(supergraph, attribute, value):
    """
    Return the subgraph induced by nodes where "attribute" value is "value"
    Returned subgraph may contain isolated nodes
    """
    
    return supergraph.subgraph([node for node in supergraph.nodes() if supergraph.nodes[node][attribute]==value]).copy()

def compute_graph_properties(supergraph, colors, excluded_colors, memo_falling_power, verbose=False):
    if verbose:
        print("\nGraph properties:")

    supergraph_degrees = nx.degree(supergraph)
    n_nodes = supergraph.number_of_nodes()
    n_edges = supergraph.number_of_edges()

    if verbose:
        print("\tNodes: {}".format(n_nodes))
        print("\tEdges: {}".format(n_edges))

    supergraph_degree_histo = nx.degree_histogram(supergraph)
    
    # Compute sum of squared degree
    sum_of_squares = 0
    for degree in range(len(supergraph_degree_histo)):
        sum_of_squares += supergraph_degree_histo[degree] * degree * (degree-1)
    sum_degrees_squared = sum_of_squares

    if verbose:
        print("\tSum of squared degrees: {}".format(sum_of_squares))
    
    col_list = [col for col in colors if col not in excluded_colors]
    col_list.sort()

    node_distribution = dict()
    edge_distribution = dict()
    for color1 in col_list:
        edge_distribution[color1] = {}
        node_distribution[color1] = 0
        for color2 in col_list:
            edge_distribution[color1][color2] = 0
    
    for node in supergraph.nodes():
        color = supergraph.nodes[node]["color"]
        if color not in excluded_colors:
            node_distribution[color] += 1
    
    for edge in supergraph.edges():
        color1 = supergraph.nodes[edge[0]]["color"]
        color2 = supergraph.nodes[edge[1]]["color"]
        if color1 not in excluded_colors and color2 not in excluded_colors:
            edge_distribution[color1][color2] += 1
            if color1 != color2:
                edge_distribution[color2][color1] += 1
    
    # Compute P3
    pi_3 = 0
    for _, degree in supergraph_degrees:
        pi_3 += falling_power(degree, 2, memo_falling_power) / 2

    if verbose:
        print("\tPi3: {}".format(pi_3))
    
    return supergraph_degrees, n_nodes, n_edges, supergraph_degree_histo, col_list, node_distribution, edge_distribution, pi_3, sum_degrees_squared, memo_falling_power

def main():
    t0 = time.time()

    # Load command line parameters
    args = read_params()
    if args.verbose:
        print('{} v{} ({})'.format(TOOL_ID, __version__, __date__))
        print("\n--input_edges {}".format(args.input_edges))
        if args.input_nodes:
            print("--input_nodes {}".format(args.input_nodes))

    # input_edges is required
    if not args.input_edges or not os.path.exists(args.input_edges):
        raise FileNotFoundError("Unable to locate input file: {}".format(args.input_edges))
    # input_nodes is optional
    if args.input_nodes and not os.path.exists(args.input_nodes):
        raise FileNotFoundError("Unable to locate input file: {}".format(args.input_nodes))

    # Reassign nproc if the provided one is greater than current cpu_count
    nproc = args.nproc
    if nproc < 1:
        raise ValueError("--nproc must be greater than 1")
    elif nproc > os.cpu_count():
        nproc = os.cpu_count()
    if args.verbose:
        print("--nproc {}".format(nproc))

    # Check whether --cmap is a valid colormap
    if args.cmap not in plt.colormaps():
        raise ValueError("--cmap is not a valid Matplotlib colormap")
    
    # Check whether the output already exists
    if not args.overwrite:
        out_files = [
            "{}_zscores_edges.txt".format(os.path.splitext(args.input_edges)[0]),
            "{}_zscores_singletons.txt".format(os.path.splitext(args.input_edges)[0]),
            "{}.pdf".format(os.path.splitext(args.input_edges)[0])
        ]
        for out_path in out_files:
            if os.path.exists(out_path):
                raise OSError("Output file already exists: {}".format(out_path))
    
    # Initialize the supergraph
    supergraph = nx.Graph()
    # Assign the input file name as the supergraph name
    supergraph.graph["name"] = os.path.splitext(os.path.basename(args.input_edges))[0]
    colors = set()
    colorGroups = set()
    utils.inputGraphFromFile(args.input_nodes,          # Input file with nodes and colors
                             args.input_edges,          # Input file with edges
                             supergraph,                # NetworkX supergraph
                             colors,                    # Set of colors
                             colorGroups,               # Set of color groups
                             args.weight_threshold,     # Weight threshold
                             isolated=args.isolated)    # Add isolated nodes

    edge_zscore_t0 = time.time()
    if args.verbose:
        print("\nComputing edge z-scores")
    excluded_colors = ""
    memo_falling_power = dict()
    memo_falling_frac = dict()
    supergraph_degrees, n_nodes, n_edges, supergraph_degree_histo, col_list, \
        node_distribution, edge_distribution, pi_3, sum_degrees_squared, memo_falling_power = compute_graph_properties(supergraph, 
                                                                                                                       colors, 
                                                                                                                       excluded_colors, 
                                                                                                                       memo_falling_power,
                                                                                                                       verbose=args.verbose)

    z_score_edges, min_value = z_score_number_edges(memo_falling_frac, n_edges, n_nodes, col_list, 
                                                    node_distribution, memo_falling_power, edge_distribution, pi_3)

    edge_zscore_t1 = time.time()
    if args.verbose:
        print("\nTotal elapsed time for computing edge z-score: {}s\n".format(int(edge_zscore_t1 - edge_zscore_t0)))
    
    print_z_score_edges(args.input_edges, col_list, z_score_edges)

    # Z-score log transformation for plotting heatmap
    if args.log_transform:
        z_score_edges = transform_z_score(z_score_edges, scale_value=args.scale_factor, from_one=args.scale_from_one)
    # Plot heatmap
    plot_heatmap(args.input_edges, z_score_edges, 
                 cmap=args.cmap, vmin=args.vmin, vmax=args.vmax, center=args.center, cbar=args.cbar)
    
    singleton_zscore_t0 = time.time()
    if args.verbose:
        print("Computing singleton z-scores")
    compute_print_z_score_singletons(args.input_edges, 
                                     n_nodes, 
                                     memo_falling_frac, 
                                     supergraph_degrees, 
                                     supergraph_degree_histo, 
                                     supergraph, 
                                     col_list, 
                                     excluded_colors, 
                                     nproc=nproc,
                                     verbose=args.verbose)
    
    singleton_zscore_t1 = time.time()
    if args.verbose:
        print("Total elapsed time for computing singleton z-scores: {}s\n".format(int(singleton_zscore_t1 - singleton_zscore_t0)))

    t1 = time.time()
    if args.verbose:
        print("Total elapsed time {}s".format(int(t1 - t0)))

if __name__ == "__main__":
    main()
