import argparse
import json
import os

import gseapy as gp

import matplotlib
matplotlib.use('Agg')  # Use a non-interactive matplotlib backend
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D

import networkx as nx
from pyvis.network import Network

import numpy as np
import pandas as pd
import re
import random

from gensim.models import Word2Vec
from sklearn.cluster import SpectralClustering


TITLE = ("""
  ██████╗ ██████╗  █████╗ ██████╗ ██╗  ██╗    ██████╗  █████╗  ██████╗██╗  ██╗
 ██╔════╝ ██╔══██╗██╔══██╗██╔══██╗██║  ██║    ██╔══██╗██╔══██╗██╔════╝██║ ██╔╝
 ██║  ███╗██████╔╝███████║██████╔╝███████║    ██████╔╝███████║██║     █████╔╝ 
 ██║   ██║██╔══██╗██╔══██║██╔═══╝ ██╔══██║    ██╔═══╝ ██╔══██║██║     ██╔═██╗ 
 ╚██████╔╝██║  ██║██║  ██║██║     ██║  ██║    ██║     ██║  ██║╚██████╗██║  ██╗
  ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝  ╚═╝    ╚═╝     ╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝""")
VIZ = ("""
    ┌───┐          ┌───┐                                              
    │███├──────────┤███│                                    ╔═══╗          
    └─┬─┘          └┬──┘                       .            ║███║  
      ╰──────╮╭─────╯           .. ............;;.          ╚═┃═╝              
            ┌┴┴─┐               ..::::::::::::;;;;.         ┌─┸─┐          
       ╭────┤   ├──╮          . . ::::::::::::;;:'          │   ├──╮       
       │    └┬──┘  │  ┌───┐                   :'            └╥──┘  │  ┌───┐
     ┌─┴─┐   │     ╰──┤░░░│                          ┌───┐   ║     ╰──┤░░░│
     │▚▚▚├───╯        └───┘                          │▚▚▚╞═══╝        └───┘
     └───┘                                           └───┘                 
""")

DESCR = ("""
GraphPack is a Python tool engineered to facilitate the compression and visualization
of large-scale networks, such as protein-protein interaction networks or metabolic
pathways. It offers a user-friendly interface that enables the application of diverse
graph compression algorithms and the visualization of the compressed networks.""")

LONG_DESCR = ("""
GraphPack includes a variety of graph compression algorithms:
- Louvain Clustering algorithm optimizes modularity to uncover community structures,
- Greedy algorithm iteratively merges nodes to minimize graph size while preserving
  structural integrity,
- Label Propagation detects communities using labels distributed across the network,
- Asynchronous Fluid Communities algorithm identifies fluid communities in the network.
- Spectral Clustering utilizes eigenvalues of the graph Laplacian for clustering,
- Hierarchical Clustering builds a hierarchy of clusters,
- Node2Vec embeds nodes through random walks for clustering purposes,
- DeepWalk learns node representations via truncated random walks,
- Clique Percolation Method (CPM) detects overlapping communities,
- Non-negative Matrix Factorization (NMF) factorizes the adjacency matrix for
  community detection.

To ensure continuity, GraphPack maps the nodes of the compressed graph to information
regarding the original nodes they represent. This mapping is maintained in an additional
annotation file, preserving the relationship between the new and old nodes.

The tool comes with example input data, including a small network, and comprehensive
documentation to guide users on how to utilize the tool effectively. GraphPack also
supports Gene Set Enrichment Analysis (GSEA) to identify enriched biological pathways
within the compressed networks, offering deeper biological insights into the network
structure.

GraphPack provides robust visualization options for both the original and compressed
networks, facilitating easy comparison and in-depth analysis.

Networks can be input in various formats, such as adjacency lists, edge lists, or
networkx graphs. Users can select a compression algorithm and specify any necessary
parameters, such as compression strength. The tool then outputs a compressed graph
along with a mapping file that details the relationship between the new and original
nodes. 
""")

# ANSI escape codes for red and bold text
RED_BOLD = "\033[1;31m"
ORANGE_BOLD = "\033[1;33m"
GREEN_BOLD = "\033[1;32m"
BOLD = "\033[1m"
RESET = "\033[0m"

FIG_NUM = 0
FIG_EXT = "pdf"
SEED = 123
random.seed(SEED)

MAX_K = 50  # Maximum number of clusters for community detection, to guarantee a nice visualization
MAX_N_NODES_INTERACTIVE = 5000  # Maximum number of nodes to enable interactive visualization
MAX_N_NODES = 5000  # 10000 # Maximum number of nodes to enable graph visualization
MAX_N_NODES_NMF = 1000  # Maximum number of nodes for NMF n_components = 'auto' community detection

METHODS = ['louvain', 'greedy', 'label_propagation', 'asyn_fluidc', 'spectral', 'hclust', 'node2vec', 'deepwalk',
           'cpm', 'nmf']

EXTENSIONS = ['edgelist', 'txt', 'csv', 'tsv', 'json', 'gpickle', 'gml',
              'graphml', 'net', 'pajek', 'gexf', 'yaml', 'yml']

COLORS = ['red', 'gold', 'mediumseagreen', 'deepskyblue', 'blue', 'darkmagenta', 'sandybrown', 'yellowgreen',
          'turquoise', 'dodgerblue', 'blueviolet', 'mediumvioletred', 'coral', 'olivedrab', 'mediumspringgreen',
          'steelblue', 'mediumslateblue', 'orchid', 'brown', 'darkorange', 'forestgreen', 'teal', 'lightsteelblue',
          'mediumorchid', 'maroon', 'orange', 'seagreen', 'cyan', 'royalblue', 'violet']

# Extended color palette
COLORS_EXT = COLORS + [
    'crimson', 'darkgoldenrod', 'lightcoral', 'springgreen', 'midnightblue', 'limegreen', 'peru', 'greenyellow',
    'cadetblue', 'slateblue', 'indigo', 'hotpink', 'firebrick', 'darkolivegreen', 'aquamarine',
    'cornflowerblue', 'darkseagreen', 'darkorchid', 'chocolate', 'orangered', 'darkgreen', 'darkcyan', 'skyblue',
    'plum', 'darkred', 'darkorange', 'lime', 'darkturquoise', 'mediumblue', 'lavender'
]

CONFIG_PATH = os.path.join(os.path.dirname(__file__), 'config')


class CustomArgumentParser(argparse.ArgumentParser):
    """
    Custom ArgumentParser class to modify the help message and usage format.

    Args:
        argparse.ArgumentParser: The ArgumentParser class to inherit from.

    Returns:
        CustomArgumentParser: A custom ArgumentParser class with modified help message and usage format.
    """

    def format_help(self):
        formatter = self._get_formatter()

        # Customize the usage message
        formatter.add_text(self.format_usage().replace("usage:", f"{BOLD}Usage:{RESET}"))

        # Add description
        formatter.add_text(self.description)

        # Customize the options message
        for action_group in self._action_groups:
            formatter.start_section(action_group.title)
            formatter.add_text(action_group.description)
            formatter.add_arguments(action_group._group_actions)
            formatter.end_section()

        # Add epilog
        formatter.add_text(self.epilog)

        return formatter.format_help()


def read_graph(file_path):
    """
    Reads a graph from a file. Supports edgelist, JSON, and various other formats.

    Args:
        file_path (str): The path to the file.

    Returns:
        nx.Graph: A NetworkX graph created from the file data.

    Raises:
        IOError: If the file cannot be read.
        ValueError: If the file extension is not supported.
        json.JSONDecodeError: If the JSON file is not valid.
        KeyError: If the JSON data does not contain the expected 'edges' key.

    Examples:
        >>> import networkx as nx
        >>> from graphpack.utils import read_graph

        >>> # Example 1: Reading an edgelist file
        >>> G = read_graph('path/to/edgelist.txt')
        >>> isinstance(G, nx.Graph)
        True
        >>> len(G.nodes) > 0  # Ensure the graph has nodes
        True

        >>> # Example 2: Reading a JSON file
        >>> G = read_graph('path/to/graph.json')
        >>> isinstance(G, nx.Graph)
        True
        >>> len(G.edges) > 0  # Ensure the graph has edges
        True

        >>> # Example 3: Reading a GML file
        >>> G = read_graph('path/to/graph.gml')
        >>> isinstance(G, nx.Graph)
        True
        >>> len(G.nodes) > 0  # Ensure the graph has nodes
        True
    """
    _, file_extension = os.path.splitext(file_path)

    if file_extension == '.json':
        return read_graph_from_json(file_path)
    elif file_extension in ['.edgelist', '.txt', '.csv', '.tsv']:
        return read_graph_from_edge_list(file_path, file_extension)
    elif file_extension == 'gpickle':
        return nx.read_gpickle(file_path)
    elif file_extension == 'gml':
        return nx.read_gml(file_path)
    elif file_extension == 'graphml':
        return nx.read_graphml(file_path)
    elif file_extension in ['net', 'pajek']:
        return nx.read_pajek(file_path)
    elif file_extension == 'gexf':
        return nx.read_gexf(file_path)
    elif file_extension in ['yaml', 'yml']:
        return nx.read_yaml(file_path)
    else:
        raise ValueError(f"Unsupported file extension: {file_extension}")


def read_graph_from_edge_list(file_path, file_extension):
    """
    Reads a graph from an edgelist, CSV, or TSV file.

    Args:
        file_path (str): The path to the file.
        file_extension (str): The extension of the file.

    Returns:
        nx.Graph: A NetworkX graph created from the file data.

    Raises:
        IOError: If the file cannot be read.

    Examples:
        >>> import networkx as nx
        >>> from graphpack.utils import read_graph_from_edge_list

        >>> # Example 1: Reading an edgelist file
        >>> # Contents of 'graph.edgelist':
        >>> # A B
        >>> # B C
        >>> G = read_graph_from_edge_list('path/to/graph.edgelist', '.edgelist')
        >>> isinstance(G, nx.Graph)
        True
        >>> len(G.nodes) == 3  # Ensure the graph has 3 nodes
        True

        >>> # Example 2: Reading a CSV file, unweighted graph
        >>> # Contents of 'graph.csv':
        >>> # source,target
        >>> # A,B
        >>> # B,C
        >>> G = read_graph_from_edge_list('path/to/graph.csv', '.csv')
        >>> isinstance(G, nx.Graph)
        True
        >>> len(G.edges) == 2  # Ensure the graph has 2 edges
        True

        >>> # Example 3: Reading a TSV file, weighted graph
        >>> # Contents of 'graph.tsv':
        >>> # source\ttarget\tweight
        >>> # A\tB\t1.0
        >>> # B\tC\t2.0
        >>> G = read_graph_from_edge_list('path/to/graph.tsv', '.tsv')
        >>> isinstance(G, nx.Graph)
        True
        >>> len(G.edges) == 2  # Ensure the graph has 2 edges
        True
        >>> G['A']['B']['weight'] == 1.0  # Ensure the graph has edge weights
        True
    """
    if file_extension in ['.edgelist', '.txt']:
        return nx.read_edgelist(file_path)

    elif file_extension in ['.csv', '.tsv']:
        delimiter = ',' if file_extension == '.csv' else '\t'
        df = pd.read_csv(file_path, delimiter=delimiter)

        if 'weight' in df.columns:
            edges = df[['source', 'target', 'weight']].to_records(index=False).tolist()
        else:
            edges = df[['source', 'target']].to_records(index=False).tolist()

        graph = nx.Graph()
        for edge in edges:
            if len(edge) == 3:
                graph.add_edge(edge[0], edge[1], weight=edge[2])
            else:
                graph.add_edge(edge[0], edge[1])

        return graph

    else:
        raise IOError(f"Unsupported file extension: {file_extension}")


def read_graph_from_json(file_path):
    """
    Reads a graph from a JSON file.

    Args:
        file_path (str): The path to the JSON file.

    Returns:
        nx.Graph: A NetworkX graph created from the JSON data.

    Raises:
        IOError: If the file cannot be read.
        json.JSONDecodeError: If the file is not valid JSON.
        KeyError: If the JSON data does not contain the expected 'edges' key.

    Examples:
        >>> import networkx as nx
        >>> from graphpack.utils import read_graph_from_json

        >>> # Example 1: Reading a JSON file with weighted edges
        >>> # Contents of 'graph.json':
        >>> # {
        >>> #     "nodes": ["A", "B", "C"],
        >>> #     "edges": [
        >>> #         {"source": "A", "target": "B", "weight": 1.0},
        >>> #         {"source": "B", "target": "C", "weight": 2.0}
        >>> #     ]
        >>> # }
        >>> G = read_graph_from_json('path/to/graph.json')
        >>> isinstance(G, nx.Graph)
        True
        >>> len(G.nodes) == 3  # Ensure the graph has nodes
        True
        >>> len(G.edges) == 2  # Ensure the graph has edges
        True
        >>> G['A']['B']['weight'] == 1.0  # Ensure the graph has weights
        True
    """
    with open(file_path, 'r') as f:
        data = json.load(f)
    graph = nx.Graph()
    graph.add_nodes_from(data["nodes"])

    for edge in data['edges']:
        graph.add_edge(edge['source'], edge['target'], weight=edge.get('weight', 1))

    return graph


def save_graph(graph, file_path, save_data=False):
    """
    Saves a NetworkX graph to a file in the format specified by the file extension.

    Args:
        graph (nx.Graph): The NetworkX graph to save.
        file_path (str): The path to the file to save the graph to.
        save_data (bool): Whether to save edge data (e.g., weights) to the file. Defaults to False.

    Raises:
        ValueError: If the file extension is not supported.

    Examples:
        >>> import networkx as nx
        >>> from graphpack.utils import save_graph

        >>> # Example 1: Save a graph to a JSON file
        >>> G = nx.Graph()
        >>> G.add_edge('A', 'B', weight=1.0)
        >>> G.add_edge('B', 'C', weight=2.0)
        >>> save_graph(G, 'path/to/graph.json', save_data=True)

        >>> # Example 2: Save a graph to an edgelist file
        >>> save_graph(G, 'path/to/graph.edgelist')

        >>> # Example 3: Save a graph to a CSV file
        >>> save_graph(G, 'path/to/graph.csv')
    """
    _, file_extension = os.path.splitext(file_path)

    if file_extension == '.json':
        data = {
            'nodes': list(graph.nodes()),
            'edges': [
                {
                    'source': u,
                    'target': v,
                    'weight': data.get('weight', 1)
                } for u, v, data in graph.edges(data=True)
            ]
        }
        with open(file_path, 'w') as f:
            json.dump(data, f, indent=4)

    elif file_extension in ['.edgelist', '.txt']:
        nx.write_edgelist(graph, file_path, data=save_data)

    elif file_extension in ['.csv', '.tsv']:
        delimiter = ',' if file_extension == '.csv' else '\t'
        if graph.number_of_edges() > 0:
            edges = [(u, v, data.get('weight', 1)) for u, v, data in graph.edges(data=save_data)]
            df = pd.DataFrame(edges, columns=['source', 'target', 'weight'])

        else:
            # If there are no edges with data, save only source and target nodes
            edges = [(u, v) for u, v in graph.edges()]
            df = pd.DataFrame(edges, columns=['source', 'target'])

        df.to_csv(file_path, sep=delimiter, index=False)

    elif file_extension == 'gpickle':
        nx.write_gpickle(graph, file_path)
    elif file_extension == 'gml':
        nx.write_gml(graph, file_path)
    elif file_extension == 'graphml':
        nx.write_graphml(graph, file_path)
    elif file_extension in ['net', 'pajek']:
        nx.write_pajek(graph, file_path)
    elif file_extension == 'gexf':
        nx.write_gexf(graph, file_path)
    elif file_extension in ['yaml', 'yml']:
        nx.write_yaml(graph, file_path)
    else:
        raise ValueError(f"Unsupported file extension: {file_extension}")


def generate_figure_filepath(title=None, file_path=None, fig_ext='png'):
    """
    Generate a file path for a figure based on title and file_path parameters.

    Args:
        title (str or None): Title of the figure. If None, a default name 'fig' will be used.
        file_path (str or None): Base directory where the figure file should be saved. If None, current directory will be used.
        fig_ext (str): Extension of the figure file (e.g., 'png', 'jpg').

    Returns:
        file_path_ext (str): Generated file path including the extension.
    """

    if title is not None:
        name = (title.lower().replace(" ", "_").replace("(", "").replace(")", "").replace(",", "").replace("- ", ""))
    else:
        name = "fig"

    if file_path is None:
        file_path = ""
    else:
        os.makedirs(file_path, exist_ok=True)  # Ensure the directory exists

    file_path_ext = os.path.join(file_path, name + '.' + fig_ext)
    return file_path_ext


def assign_community_colors(graph, compressed_graph, decompression_mapping, labels=None):
    """
    Assigns colors to nodes in the original graph based on community labels of the compressed graph.

    Args:
        graph (nx.Graph): The original graph.
        compressed_graph (nx.Graph): The compressed graph.
        decompression_mapping (dict): Mapping from nodes in compressed graph to nodes in original graph.
        labels (dict): Community labels for nodes in the compressed graph. Defaults to None.

    Returns:
        partition_colors (list): List of colors assigned to communities in the compressed graph.
        node_colors (list): List of colors assigned to nodes in the original graph.
        color_map (dict): Mapping from community labels to colors.
    """
    # Extract unique communities from the compressed graph
    communities = list(labels.values()) if labels else [community_id for community_id in compressed_graph.nodes]

    # Create a color map from community IDs/labels to colors
    colors = [mcolors.CSS4_COLORS[name] for name in COLORS_EXT]
    color_map = {community: colors[i % len(colors)] for i, community in enumerate(communities)}

    # Assign colors to nodes in the original graph based on the decompression mapping and labels
    if labels is not None:
        node_colors = [color_map[labels[decompression_mapping[node]]] for node in graph.nodes()]
        partition_colors = [color_map[community] for community in communities]
    else:
        node_colors = [color_map[decompression_mapping[node]] for node in graph.nodes()]
        partition_colors = [color_map[community] for community in communities]

    return partition_colors, node_colors, color_map


def draw_interactive_graph_pyvis(graph, labels=None, title=None, file_path=None, node_color=None, node_sizes=None,
                                 color_map=None):
    """
    Draws an interactive graph using config and saves it to an HTML file.

    Args:
        graph (nx.Graph): The graph to be drawn.
        labels (dict, optional): Node labels. Defaults to None.
        title (str, optional): Title of the graph plot. Defaults to None.
        file_path (str, optional): File path to save the plot as an HTML file. If None, the plot is saved in the current directory with a default name. Defaults to None.
        node_color (list, optional): Node colors. Defaults to None.
        node_sizes (list, optional): Node sizes. Defaults to None.
        color_map (dict, optional): Color mapping for legend. Defaults to None.

    Returns:
        None

    Examples:
        >>> import networkx as nx
        >>> from graphpack.utils import draw_interactive_graph_pyvis

        >>> # Example 1: Draw a basic graph
        >>> G = nx.Graph()
        >>> G.add_edge('A', 'B')
        >>> G.add_edge('B', 'C')
        >>> draw_interactive_graph_pyvis(G, title='Basic Graph', file_path='simple_graph')

        >>> # Example 2: Draw a graph with custom node labels
        >>> labels = {'A': 'Node A', 'B': 'Node B', 'C': 'Node C'}
        >>> draw_interactive_graph_pyvis(G, labels=labels, title='Graph with Labels', file_path='simple_graph')

        >>> # Example 3: Draw a graph with custom node colors
        >>> node_colors = ['red', 'green', 'blue']
        >>> draw_interactive_graph_pyvis(G, node_color=node_colors, title='Graph with Colors', file_path='simple_graph')

        >>> # Example 4: Draw a graph with custom node sizes
        >>> node_sizes = [10, 20, 30]
        >>> draw_interactive_graph_pyvis(G, node_sizes=node_sizes, title='Graph with Sizes', file_path='simple_graph')
    """
    # Create a config Network
    net = Network(height="800px", width="100%", font_color="gray", filter_menu=True, select_menu=True,
                  cdn_resources='remote')
    node_color = node_color if node_color else ['dodgerblue'] * len(graph.nodes())
    color_map = {v: k for k, v in color_map.items()} if color_map else None

    # Normalize node sizes to range from 0 to 1
    if node_sizes:
        min_node_size = 10
        max_node_size = 20
        min_size, max_size = min(node_sizes), max(node_sizes)
        if min_size != max_size:
            node_sizes = [((size - min_size) / (max_size - min_size) * (max_node_size - 1) + min_node_size) for size in
                          node_sizes]
        else:
            node_sizes = [10] * len(graph.nodes())
    else:
        node_sizes = [10] * len(graph.nodes())

    # Add nodes and edges from the NetworkX graph
    for i, node in enumerate(graph.nodes()):
        color = node_color[i]
        label = labels[node] if labels else str(node)
        size = node_sizes[i]
        title_attr = f"{color_map[color]}" if color_map and color in color_map else str(
            node)  # f"{node} - {color_map[color]}" if color_map and color in color_map else str(node)

        net.add_node(node, label=label, color=color, size=size, title=title_attr)

    # Extract the edges and add them to the config Network
    for edge in graph.edges():
        width = graph[edge[0]][edge[1]].get('weight', 1) / 2  # Default to 1 if weight is None or not present
        net.add_edge(edge[0], edge[1], width=width, color="lightgray")

    # Set buttons (uncomment ONLY to modify the options in GUI mode)
    # net.show_buttons()

    # Define the common part of the options, as specified in config/const_options file
    with open(f"{CONFIG_PATH}/const_options", "r") as file:
        common_options = file.read()

    # Define the conditional part of the physics options
    physics_option_file = f"{CONFIG_PATH}/physics_options_small" if len(
        net.nodes) < 30 else f"{CONFIG_PATH}/physics_options_big"
    with open(physics_option_file, "r") as file:
        physics_options = file.read()

    # Concatenate the common options with the conditional physics options
    options = common_options + physics_options

    # Set the options (comment if using GUI mode for options setup)
    net.set_options(options)

    # Save the plot to an image file if file_path is specified
    file_path_ext = generate_figure_filepath(title, file_path, 'html')
    html = net.generate_html(file_path_ext)

    with open(file_path_ext, "w+") as out:
        out.write(html)

    def add_physics_switch(filepath):
        with open(filepath, 'r', encoding="utf-8") as file:
            content = file.read()

        # Define the pattern to match the existing <style> section
        pattern = r'<style\s+type="text/css">.*?</style>'

        # Define the replacement content for the <style> section
        with open(f"{CONFIG_PATH}/style_options", "r") as file:
            replacement = file.read()

        # Perform the replacement in the content
        new_content = re.sub(pattern, replacement, content, flags=re.DOTALL)

        # Write the modified content back to the file
        with open(filepath, 'w', encoding="utf-8") as file:
            file.write(new_content)

    def add_physics_stop_to_html(filepath):
        with open(filepath, 'r', encoding="utf-8") as file:
            content = file.read()

        # Search for the stabilizationIterationsDone event and insert the network.setOptions line
        pattern = r'(network.once\("stabilizationIterationsDone", function\(\) {)'
        replacement = r'\1\n\t\t\t\t\t\t  // Disable the physics after stabilization is done.\n\t\t\t\t\t\t  network.setOptions({ physics: false });\n\t\t\t\t\t\t  updateCheckboxState();'

        new_content = re.sub(pattern, replacement, content, flags=re.DOTALL)

        # Write the modified content back to the file
        with open(filepath, 'w', encoding="utf-8") as file:
            file.write(new_content)

    def add_unselected_nodes_hiding(filepath):
        with open(filepath, 'r', encoding="utf-8") as file:
            content = file.read()

            # Function to add JavaScript code for handling node click and reset visibility
            with open(f"{CONFIG_PATH}/additional_js", "r") as file:
                js_code = file.read()

            # Insert the JavaScript code before the closing </body> tag
            new_content = content.replace('</body>', js_code + '</body>')

            # Write the modified content back to the file
            with open(filepath, 'w', encoding="utf-8") as file:
                file.write(new_content)

    # Set options such as layout, physics, etc. in the HTML file
    add_physics_switch(file_path_ext)
    add_physics_stop_to_html(file_path_ext)
    add_unselected_nodes_hiding(file_path_ext)


def draw_graph(graph, labels=None, edge_thickness=None, node_color=None, color_map=None, title=None, file_path=None,
               is_interactive=False, plot_disconnected_components=False, separate_communities=False, **kwargs):
    """
    Draws a graph using NetworkX and Matplotlib and saves it to an image file.

    Args:
        graph (nx.Graph): The graph to be drawn.
        labels (dict, optional): Node labels. Defaults to None.
        edge_thickness (list, optional): List of edge thicknesses. Defaults to None.
        node_color (list, optional): Node colors. Defaults to None.
        color_map (dict, optional): Color mapping for legend. Defaults to None.
        title (str, optional): Title of the graph plot. Defaults to None.
        file_path (str, optional): File path to save the plot as an image. Defaults to None.
        is_interactive (bool, optional): Whether to draw an interactive plot using config. Defaults to False.
        plot_disconnected_components (bool, optional): How to handle disconnected graphs.
            - True: Draw the entire graph including disconnected nodes.
            - False: Draw only the largest connected component.
        separate_communities (bool, optional): Whether to separate communities in the layout. Defaults to False.
        **kwargs: Additional keyword arguments to be passed to `plt.figure`.

    Returns:
        None

    Examples:
        >>> import networkx as nx
        >>> import matplotlib.pyplot as plt
        >>> from graphpack.utils import draw_graph, draw_interactive_graph_pyvis

        >>> # Example 1: Simple graph with default settings
        >>> G = nx.Graph()
        >>> G.add_edges_from([(1, 2), (2, 3), (3, 1)])
        >>> draw_graph(G, title="Simple Graph", file_path="simple_graph")

        >>> # Example 2: Graph with specified node colors and edge thickness
        >>> node_colors = ['blue', 'green', 'red']
        >>> edge_thickness = [2.0, 1.5, 2.5]
        >>> draw_graph(G, node_color=node_colors, edge_thickness=edge_thickness, title="Graph with Custom Colors and Thickness", file_path="simple_graph")

        >>> # Example 3: Interactive graph using config
        >>> draw_graph(G, is_interactive=True, title="Interactive Graph", file_path="simple_graph")

        >>> # Example 4: Graph with communities and color mapping
        >>> communities = {1: 'blue', 2: 'blue', 3: 'red'}
        >>> draw_graph(G, node_color=[communities[node] for node in G.nodes()], color_map={'Blue Nodes': 'blue', 'Red Nodes': 'red'}, title="Graph with Communities", file_path="simple_graph")

        >>> # Example 5: Handling disconnected components
        >>> G.add_edges_from([(4, 5), (5, 6)])
        >>> draw_graph(G, plot_disconnected_components=True, title="Graph with Disconnected Components", file_path="simple_graph")

        >>> # Final example
        >>> G = nx.Graph()
        >>> G.add_edges_from([(1, 2), (2, 3), (3, 1), (4, 5), (5, 6), (6, 7), (7, 8), (8, 6)])
        >>> communities = {1: 'blue', 2: 'blue', 3: 'yellow', 4: 'green', 5: 'green', 6: 'red', 7: 'red', 8: 'red'}
        >>> edge_thickness = [2.0, 1.5, 2.5, 1.0, 2.0, 1.5, 2.5, 5.0]
        >>> color_map = {'Blue Nodes': 'blue', 'Red Nodes': 'red', 'Green Nodes': 'green', 'Yellow Nodes': 'yellow'}
        >>> draw_graph(G, labels=None, edge_thickness=edge_thickness, node_color=[communities[node] for node in G.nodes()], color_map=color_map, title="Example Graph with All Options Set", file_path="simple_graph", is_interactive=True, plot_disconnected_components=True, separate_communities=True)
    """
    # Check if the graph is connected
    if not nx.is_connected(graph) and plot_disconnected_components is False:
        print(
            f"{RED_BOLD}The graph is not connected. Only the largest connected component will be displayed, thus some nodes will not be present in the graph.{RESET}")
        print(f"{BOLD}To display all nodes, consider disabling the connected_component flag.{RESET}")

    # Handle connected components
    if not plot_disconnected_components:
        # Use only the largest connected component
        largest_cc = max(nx.connected_components(graph), key=len)
        subgraph = graph.subgraph(largest_cc).copy()
    else:
        # Use the entire graph including disconnected nodes
        subgraph = graph

    # If number of nodes exceeds MAX_N_NODES_INTERACTIVE, remove all leaf nodes
    if subgraph.number_of_nodes() > MAX_N_NODES_INTERACTIVE:
        print(f"\n{ORANGE_BOLD}Leaf nodes and nodes with degree <= 2 will be removed to enable visualization.{RESET}")
        print(f"{BOLD}To display all nodes, consider increasing the MAX_N_NODES_INTERACTIVE constant.{RESET}")
        leaves = [node for node, degree in dict(subgraph.degree()).items() if degree <= 2]
        subgraph.remove_nodes_from(leaves)

        print(f"{ORANGE_BOLD}Removed {len(leaves)} leaves or nodes with degree 2.{RESET}\n")

    # Update node_color list to match remaining nodes
    if node_color:
        remaining_nodes = set(subgraph.nodes)
        node_color = [color for node, color in zip(graph.nodes, node_color) if node in remaining_nodes]

    # Get the edges with their data (attributes)
    edges = subgraph.edges(data=True)
    num_edges = len(edges)
    num_nodes = len(subgraph.nodes)

    # Define a base thickness scaling factor inversely related to the number of edges
    scaling_factor = 1 / (4 * (num_edges ** 0.5)) if num_edges > 50 else 5

    # Check if edge_thickness is not provided
    if edge_thickness is None:
        # Check if any edge has a 'weight' attribute
        if any('weight' in data for _, _, data in edges):
            # Extract all weights from the edges
            weights = [data['weight'] for _, _, data in edges]
            # If there's variation in the weights, normalize the thickness
            if min(weights) != max(weights):
                # Normalize weights to range [0.1, 1] for thickness and apply scaling factor
                edge_thickness = [scaling_factor * (0.1 + (weight - min(weights)) / (max(weights) - min(weights))) for
                                  weight in weights]
            else:
                # If all weights are the same, set a default thickness and apply scaling factor
                edge_thickness = [scaling_factor for _ in weights]
        else:
            # If no weights are specified, set all thicknesses to a default value and apply scaling factor
            edge_thickness = [scaling_factor for _ in edges]

    # Ensure edge_thickness is a list, even if a single value was given
    if isinstance(edge_thickness, (int, float)):
        edge_thickness = [edge_thickness * scaling_factor] * len(edges)

    # Assign as subgraph edges weights the edge thicknesses
    for idx, (u, v, data) in enumerate(subgraph.edges(data=True)):
        subgraph[u][v]['weight'] = edge_thickness[idx]

    # Initialize edge colors
    edge_colors = ['darkgray' for _ in range(num_edges)]

    # Calculate node size
    base_size = 5000  # base size for nodes
    min_node_size = 25.0 if num_nodes > 50 else 60.0
    max_node_size = max(min_node_size, base_size / num_nodes)

    # Calculate the label size
    label_size = max(1, 10 - num_nodes // 100)

    # Extract node sizes from the 'size' attribute, with a default minimum size
    node_sizes_raw = [subgraph.nodes[node].get('size', 1) for node in subgraph.nodes()]

    # Normalize node sizes to range from 1 to max_node_size
    min_size, max_size = min(node_sizes_raw), max(node_sizes_raw)
    if min_size != max_size:
        node_sizes = [(((size - min_size) / (max_size - min_size) * (max_node_size - 1) + min_node_size) * 5) for size
                      in node_sizes_raw]
    else:
        node_sizes = [max_node_size for _ in node_sizes_raw]

    # Assign as subgraph nodes weights the node sizes
    for idx, node in enumerate(subgraph.nodes()):
        subgraph.nodes[node]['size'] = node_sizes[idx]

    # Layout the subgraph
    if color_map and separate_communities:
        # Precompute clusters
        unique_colors = list(set(node_color))
        unique_clusters = list(set(color_map.values()))
        if len(unique_clusters) > MAX_K:
            print(
                f"{ORANGE_BOLD}Warning: Number of unique clusters ({len(unique_clusters)}) is greater than the  ({MAX_K}).")

        # Determine grid size
        num_clusters = len(unique_colors)
        grid_size = int(np.ceil(np.sqrt(num_clusters)))

        # Generate grid positions for clusters
        color_positions = {}
        for i, color in enumerate(unique_colors):
            row = i // grid_size
            col = i % grid_size
            color_positions[color] = np.array([col, row])

        # Apply the custom layout
        pos = nx.spring_layout(subgraph, scale=0.5)
        for node, color in zip(subgraph.nodes(), node_color):
            pos[node] += color_positions[color]  # Shift position based on color

        # Adjust edge colors for inter-community edges
        for i, (u, v) in enumerate(subgraph.edges()):
            if node_color[list(subgraph.nodes()).index(u)] != node_color[list(subgraph.nodes()).index(v)]:
                edge_colors[i] = '#E4E4E4'  # Lighter color for inter-community edges

    else:
        pos = nx.spring_layout(subgraph, scale=1, k=1 / num_nodes)

    # Filter labels for nodes present in the graph
    if labels:
        labels = {node: labels.get(node, node) for node in subgraph.nodes()}
        # Check node colors consistency

    if node_color:
        if len(node_color) != len(subgraph.nodes):
            raise ValueError(
                f"Number of node colors ({len(node_color)}) does not match number of nodes in the graph ({len(subgraph.nodes)}).")

    # If interactive, draw the graph using config
    if is_interactive:
        draw_interactive_graph_pyvis(subgraph,
                                     labels=labels,
                                     title=title,
                                     file_path=file_path,
                                     node_color=node_color,
                                     node_sizes=node_sizes,
                                     color_map=color_map)

    # Create the figure
    plt.figure(**kwargs)

    # Set the title if provided
    if title is not None:
        plt.title(title)

    nx.draw(subgraph, pos, with_labels=True, labels=labels,
            font_size=label_size, node_size=node_sizes, width=edge_thickness,
            edge_color=edge_colors, node_color=node_color)

    # Add legend for color mapping if provided
    if color_map is not None:
        color_handles = [Line2D([0], [0], marker='o', color='w', label=label,
                                markersize=10, markerfacecolor=color)
                         for label, color in color_map.items()]

        # Disable the legend if there are too many communities
        if len(color_handles) <= MAX_K:
            plt.legend(handles=color_handles, title='Communities', loc='upper left', bbox_to_anchor=(1, 1))

    # Save the plot to an image file if file_path is specified
    file_path_ext = generate_figure_filepath(title, file_path, FIG_EXT)

    # Set figure size and dpi (adjust as needed)
    plt.gcf().set_size_inches((20, 14))
    plt.gcf().set_dpi(200)

    # Save figure as PDF with adjusted bounding box
    plt.savefig(file_path_ext, bbox_inches='tight')
    plt.close()


def generate_random_walks(graph, num_walks, walk_length):
    """
    Generates random walks from each node in the graph.

    Args:
        graph (nx.Graph): The input graph.
        num_walks (int): Number of random walks to perform per node.
        walk_length (int): Length of each random walk.

    Returns:
        list of list of str: A list of random walks.
    """
    walks = []
    nodes = list(graph.nodes())

    for _ in range(num_walks):
        random.shuffle(nodes)
        for node in nodes:
            walks.append(random_walk(graph, node, walk_length))

    return walks


def random_walk(graph, start_node, walk_length):
    """
    Performs a random walk starting from the given node.

    Args:
        graph (nx.Graph): The input graph.
        start_node (str or int): The starting node for the random walk.
        walk_length (int): Length of the random walk.

    Returns:
        list of str or int: A random walk.
    """
    walk = [start_node]

    for _ in range(walk_length - 1):
        cur = walk[-1]
        neighbors = list(graph.neighbors(cur))
        if len(neighbors) > 0:
            next_node = random.choice(neighbors)
            walk.append(next_node)
        else:
            break

    return walk


def deepwalk_embedding(graph, num_walks, walk_length):
    """
    Generates node embeddings using the DeepWalk algorithm.

    This function performs random walks on the input graph and learns node embeddings using the Word2Vec algorithm.

    Args:
        graph (nx.Graph): The input graph.
        num_walks (int): Number of random walks to perform per node.
        walk_length (int): Length of each random walk.

    Returns:
        gensim.models.word2vec.Word2Vec: A Word2Vec model trained on the generated random walks.
    """
    # Step 1: Generate random walks
    walks = generate_random_walks(graph, num_walks, walk_length)

    # Convert each walk to a list of strings (Word2Vec expects strings)
    walks = [[str(node) for node in walk] for walk in walks]

    # Step 2: Train Word2Vec model
    model = Word2Vec(walks, vector_size=64, window=10, min_count=0, sg=1, workers=4)

    return model


def cluster_graph_embeddings(graph, model, k, seed=SEED):
    """
    Generate embeddings for the graph nodes and perform spectral clustering.

    Args:
        graph (nx.Graph): The input graph.
        model (gensim.models.Word2Vec): The trained embedding model.
        k (int): The number of clusters.
        seed (int, optional): The random seed for reproducibility. Default to 123.

    Returns:
        dict: A dictionary where keys are node IDs and values are cluster labels.

    Examples:
        >>> import networkx as nx
        >>> from gensim.models import Word2Vec
        >>> from sklearn.cluster import SpectralClustering
        >>> import numpy as np
        >>> from graphpack.utils import cluster_graph_embeddings

        >>> # Example 1: Simple graph clustering
        >>> G = nx.karate_club_graph()
        >>> model = Word2Vec(sentences=[[str(node) for node in G.neighbors(n)] for n in G.nodes()], vector_size=16, window=5, min_count=1, sg=1)
        >>> partition = cluster_graph_embeddings(G, model, k=2)
        >>> print(partition)
        {0: 1, 1: 0, 2: 1, 3: 0, 4: 1, 5: 1, 6: 1, 7: 0, 8: 0, 9: 0, 10: 1, 11: 0, 12: 1, 13: 1, 14: 0, 15: 0, 16: 0, 17: 0, 18: 1, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 1, 26: 1, 27: 0, 28: 1, 29: 1, 30: 1, 31: 0, 32: 1, 33: 0}

        >>> # Example 2: Clustering with different number of clusters
        >>> partition = cluster_graph_embeddings(G, model, k=4)
        >>> print(partition)
        {0: 0, 1: 1, 2: 0, 3: 1, 4: 2, 5: 0, 6: 2, 7: 3, 8: 2, 9: 2, 10: 0, 11: 1, 12: 0, 13: 0, 14: 3, 15: 3, 16: 1, 17: 1, 18: 0, 19: 3, 20: 3, 21: 1, 22: 3, 23: 1, 24: 2, 25: 2, 26: 3, 27: 3, 28: 2, 29: 0, 30: 0, 31: 3, 32: 2, 33: 1}
    """
    # Generate embeddings for the graph nodes
    embeddings = np.array([model.wv[str(node)] for node in graph.nodes()])

    # Perform spectral clustering on the embeddings
    sc = SpectralClustering(n_clusters=k, random_state=seed)
    labels = sc.fit_predict(embeddings)

    # Create a partition dictionary with node IDs as keys and cluster labels as values
    partition = {node: int(labels[i]) for i, node in enumerate(graph.nodes())}
    return partition


def perform_gsea(gene_list, organism='human', gene_sets='KEGG_2019_Human', k=5):
    """
    Performs Gene Set Enrichment Analysis (GSEA) on a list of nodes.

    Args:
        gene_list (list): A list of genes to perform GSEA on.
        gene_sets (str): Name of the gene sets to use for enrichment analysis. Default is 'KEGG_2019_Human'.
        k (int): Number of top enriched terms to return. Default is 5.

    Returns:
        list: A list of labels representing the top k enriched gene sets.

    Raises:
        ValueError: If the gene list is empty.

    Examples:
        >>> from graphpack.utils import perform_gsea

        >>> # Example 1: Basic GSEA with default settings
        >>> gene_list = ['TP53', 'BRCA1', 'EGFR', 'MYC', 'MTOR']
        >>> perform_gsea(gene_list)
        ['Breast cancer', 'Central carbon metabolism in cancer', 'MicroRNAs in cancer', 'Colorectal cancer', 'PI3K-Akt signaling pathway']

        >>> # Example 2: GSEA with a different organism and gene set
        >>> perform_gsea(gene_list, organism='mouse', gene_sets='KEGG_2019_Mouse')
        ['Breast cancer', 'MicroRNAs in cancer', 'Central carbon metabolism in cancer', 'PI3K-Akt signaling pathway', 'ErbB signaling pathway']

        >>> # Example 3: GSEA with a custom number of top enriched terms
        >>> perform_gsea(gene_list, k=3)
        ['Breast cancer', 'Central carbon metabolism in cancer', 'MicroRNAs in cancer']

        >>> # Example 4: Handling gene list with fewer than 4 genes
        >>> perform_gsea(['TP53'])
        ['TP53']
        >>> perform_gsea(['TP53', 'BRCA1'])
        [['TP53', 'BRCA1']]
    """
    # Check if the gene list is empty
    if not gene_list:
        raise ValueError("The gene list is empty. Please provide a list of genes for analysis.")

    # Handle gene lists with fewer than 4 genes
    if len(gene_list) <= 3:
        return [gene_list] if len(gene_list) > 1 else gene_list

    # Perform GSEA
    results = gp.enrichr(gene_list=gene_list, organism=organism, gene_sets=gene_sets)

    # Get the top k enriched terms
    if not results.res2d.empty:
        top_enriched_terms = results.res2d.iloc[:k]['Term'].tolist()
        return top_enriched_terms
    else:
        return ["No Enriched Gene Set"]
