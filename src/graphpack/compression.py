from collections import defaultdict

import msgpack
import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities, label_propagation_communities, asyn_fluidc, \
    k_clique_communities, louvain_communities
from node2vec import Node2Vec
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import NMF
from sklearn.impute import SimpleImputer
from tqdm import tqdm

from graphpack.utils import *


def compute_adjacency_matrix(graph):
    """
    Computes the adjacency matrix of a graph.

    Args:
        graph: The input graph.

    Returns:
        np.ndarray: The adjacency matrix of the input graph.
    """
    adj_matrix = nx.to_numpy_array(graph)

    if np.isnan(adj_matrix).any():
        # Impute missing values
        imputer = SimpleImputer(strategy='mean')  # You can also use 'median', 'most_frequent', or 'constant'
        adj_matrix = imputer.fit_transform(adj_matrix)

    return adj_matrix


def detect_communities(graph, method='louvain', seed=123, resolution=1.25, k=3):
    """
    Detects communities in a graph using various partitioning methods.

    Args:
        graph (nx.Graph): The graph to detect communities in.
        method (str, optional): The partitioning method to use. Options: 'louvain' (default), 'greedy', 'label_propagation', 'asyn_fluidc', 'spectral', 'hclust', 'node2vec', 'deepwalk', 'cpm', 'nmf'.
        seed (int, optional): Random seed for reproducibility. Defaults to 123.
        resolution (float, optional): Resolution parameter for Louvain and greedy methods. Defaults to 1.25.
        k (int, optional): Number of communities/clusters/components for 'asyn_fluidc', 'spectral', 'hclust', 'node2vec', 'deepwalk', 'nmf' methods. For 'cpm' is the size of the smallest clique. Defaults to 3.

    Returns:
        dict: A dictionary mapping nodes to community IDs.

    Raises:
        ValueError: If the specified method is not supported or if the graph is not connected.

    Example:
        >>> import networkx as nx
        >>> from graphpack.compression import detect_communities

        >>> # Example: Louvain method
        >>> G = nx.Graph()
        >>> G.add_edges_from([(1, 2), (2, 3), (3, 1), (4, 5), (5, 6), (6, 7), (7, 8), (8, 6)])
        >>> partition = detect_communities(G, method='louvain')
        >>> print(partition)
        {1: 1, 2: 1, 3: 1, 4: 0, 5: 0, 6: 2, 7: 2, 8: 2}
    """
    partition = {}
    if method == 'louvain':
        communities = list(louvain_communities(graph, resolution=resolution, seed=seed))
        partition = {node: cid for cid, community in enumerate(communities) for node in community}

    elif method == 'greedy':
        communities = list(greedy_modularity_communities(graph,
                                                         resolution=resolution,
                                                         cutoff=min(k, graph.number_of_nodes()),
                                                         best_n=min(MAX_K, graph.number_of_nodes())))
        partition = {node: cid for cid, community in enumerate(communities) for node in community}

    elif method == 'label_propagation':
        communities = list(label_propagation_communities(graph))
        partition = {node: cid for cid, community in enumerate(communities) for node in community}

    elif method == 'asyn_fluidc':
        if nx.is_connected(graph) is False:
            raise ValueError("Graph must be connected for community detection.")

        communities = list(asyn_fluidc(graph, k=k, seed=seed))
        partition = {node: cid for cid, community in enumerate(communities) for node in community}

    elif method == 'spectral':
        adj_matrix = compute_adjacency_matrix(graph)
        sc = SpectralClustering(n_clusters=k, affinity='precomputed', random_state=seed)
        labels = sc.fit_predict(adj_matrix)
        partition = {node: int(labels[i]) for i, node in enumerate(graph.nodes())}

    elif method == 'hclust':
        adj_matrix = compute_adjacency_matrix(graph)
        hc = AgglomerativeClustering(n_clusters=k, metric='precomputed', linkage='average')
        labels = hc.fit_predict(adj_matrix)
        partition = {node: int(labels[i]) for i, node in enumerate(graph.nodes())}

    elif method == 'node2vec':
        node2vec = Node2Vec(graph, dimensions=64, walk_length=30, num_walks=200, workers=4, seed=seed)
        model = node2vec.fit(window=10, min_count=1, batch_words=4)
        partition = cluster_graph_embeddings(graph, model, k, seed)

    elif method == 'deepwalk':
        model = deepwalk_embedding(graph, walk_length=80, num_walks=10)
        partition = cluster_graph_embeddings(graph, model, k, seed)

    elif method == 'cpm':
        for i in range(3):
            communities = list(k_clique_communities(graph, k=k))
            partition = {node: cid for cid, community in enumerate(communities) for node in community}

            # Include nodes not in any k-clique community
            unassigned_nodes = set(graph.nodes) - set(partition.keys())
            for node in unassigned_nodes:
                partition[node] = len(communities)

    elif method == 'nmf':
        adj_matrix = compute_adjacency_matrix(graph)
        model = NMF(n_components=k, random_state=seed)
        w = model.fit_transform(adj_matrix)
        h = model.components_
        labels = np.argmax(w, axis=1)
        partition = {node: int(labels[i]) for i, node in enumerate(graph.nodes())}

    else:
        raise ValueError(f"Unsupported detection method: {method}")

    return partition


def compress_graph_partition_based(graph, partition, is_weighted=True, is_gene_network=True):
    """
    Compresses a network based on detected communities.

    Args:
        graph (nx.Graph): The original network graph.
        partition (dict): A dictionary mapping nodes to community IDs.
        is_weighted (bool, optional): Whether to consider edge weights. Defaults to True.
        is_gene_network (bool, optional): Whether the network is a gene network to perform GSEA. Defaults to True.

    Returns:
        tuple: A tuple containing:
            - community_graph (nx.Graph): The compressed network graph where nodes represent communities.
            - compression_mapping (dict): A mapping of original nodes to their corresponding community nodes.
            - decompression_mapping (dict): A mapping of community nodes to the original nodes they represent.

    Raises:
        ValueError: If the partition is invalid or missing nodes.

    Example:
        >>> import networkx as nx
        >>> import random
        >>> random.seed(123)
        >>> from networkx.algorithms.community import louvain_communities
        >>> from graphpack.compression import compress_graph_partition_based

        >>> # Example 1: Compress an unweighted graph based on detected communities
        >>> G = nx.Graph()
        >>> G.add_edges_from([(1, 2), (2, 3), (3, 1), (4, 5), (5, 6), (6, 7), (7, 8), (8, 6)])
        >>> communities = louvain_communities(graph)
        >>> partition = {node: cid for cid, community in enumerate(communities) for node in community}
        >>> compressed_graph, _, _ = compress_graph_partition_based(G, partition, is_weighted=False)
        >>> print(compressed_graph.nodes(data=True))
        [(1, {'nodes': [1, 2, 3], 'label': [[1, 2, 3]], 'size': 3}), (0, {'nodes': [4, 5], 'label': [[4, 5]], 'size': 2}), (2, {'nodes': [6, 7, 8], 'label': [[6, 7, 8]], 'size': 3})]
        >>> print(compressed_graph.edges(data=True))
        [(0, 2, {'weight': 1})]

        >>> # Example 2: Compress a weighted gene network based on detected communities
        >>> weights = {edge: random.random() for edge in G.edges()}
        >>> nx.set_edge_attributes(G, weights, 'weight')
        >>> compressed_graph, compression_mapping, decompression_mapping = compress_graph_partition_based(G, partition)
        >>> print(compressed_graph.edges(data=True))
        [(0, 2, {'weight': 0.9011988779516946})]
        >>> print(compression_mapping)
        {1: [1, 2, 3], 0: [4, 5], 2: [6, 7, 8]}
        >>> print(decompression_mapping)
        {1: 1, 2: 1, 3: 1, 4: 0, 5: 0, 6: 2, 7: 2, 8: 2}
    """
    # Total number of iterations for progress bar
    n_iterations = len(set(partition.values())) * 2 if is_gene_network else len(set(partition.values()))

    pbar = tqdm(desc="📉 Compressing network (adding nodes)", total=n_iterations)

    try:
        # Create a dictionary to hold the list of nodes for each community
        community_nodes = {}
        community_labels = {}
        for node, community_id in partition.items():
            if community_id not in community_nodes:
                community_nodes[community_id] = []
            community_nodes[community_id].append(node)

        # Perform GSEA on the communities to assign biologically meaningful labels
        if is_gene_network:
            for community_id, nodes in community_nodes.items():
                try:
                    community_labels[community_id] = perform_gsea(nodes)
                except Exception as e:
                    print(f"Error performing GSEA: {e}")
                    community_labels[community_id] = nodes

                pbar.update(1)  # Update progress bar for each community

        # Create a new graph where nodes are communities represented by their node lists
        community_graph = nx.Graph()

        # Add community nodes (each community as a list of nodes)
        for community_id, nodes in community_nodes.items():
            size = len(nodes)  # Calculate the size based on the number of original nodes
            community_graph.add_node(community_id,
                                     nodes=nodes,
                                     label=community_labels[community_id] if is_gene_network else None,
                                     size=size)  # Assign size to each community node

            pbar.update(1)  # Update progress bar for each community

        pbar.close()

        n_iterations = len(graph.edges)
        pbar = tqdm(desc="📉 Compressing network (computing new edge weights)", total=n_iterations)

        # Add edges between communities with appropriate weights
        edge_weights = defaultdict(int)

        for node1, node2, data in graph.edges(data=True):
            community1 = partition[node1]
            community2 = partition[node2]
            weight = data['weight'] if is_weighted else 1  # Weight from the original graph, or 1 if not weighted

            if community1 != community2:
                # Use a tuple (min, max) to ensure each edge is uniquely identified
                edge = (min(community1, community2), max(community1, community2))
                edge_weights[edge] += weight

            pbar.update(1)  # Update progress bar for each graph edge

        pbar.close()

        n_iterations = len(edge_weights.items())
        pbar = tqdm(desc="📉 Compressing network (adding edges)", total=n_iterations)

        # Add the edges to the community graph
        for (community1, community2), weight in edge_weights.items():
            community_graph.add_edge(community1, community2, weight=weight)
            pbar.update(1)  # Update progress bar for each edge weight

        # Create compression mapping: {community_id: [nodes]}
        compression_mapping = {community_id: nodes for community_id, nodes in community_graph.nodes(data='nodes')}

        # Create decompression mapping: {old_node: new_node}
        decompression_mapping = {old_node: new_node for new_node, old_nodes in compression_mapping.items() for old_node
                                 in old_nodes}

        return community_graph, compression_mapping, decompression_mapping

    except KeyError as e:
        raise ValueError(f"Invalid partition or missing nodes: {e}")

    finally:
        pbar.close()


def compress_graph(graph, method='louvain', seed=123, resolution=1.25, k=3, is_weighted=True, is_gene_network=True):
    """
    Compresses a graph using the specified method and generates a new graph where nodes represent communities.

    Args:
        graph (nx.Graph): The original network graph.
        method (str, optional): The compression method to use. Options: 'louvain' (default), 'greedy', 'label_propagation', 'asyn_fluidc', 'spectral', 'hclust', 'node2vec', 'deepwalk', 'cpm', 'nmf'.
        seed (int, optional): Random seed for reproducibility. Defaults to 123.
        resolution (float, optional): Resolution parameter for Louvain and greedy methods. Defaults to 1.25.
        k (int, optional): Number of communities/clusters/components for 'asyn_fluidc', 'spectral', 'hclust', 'node2vec', 'deepwalk' methods. For 'cpm' is the size of the smallest clique. Defaults to 3.
        is_weighted (bool, optional): Whether to consider edge weights. Defaults to True.
        is_gene_network (bool, optional): Whether the network is a gene network to perform GSEA. Defaults to True.

    Returns:
        tuple: A tuple containing:
            - compressed_graph (nx.Graph): The compressed network graph where nodes represent communities.
            - compression_mapping (dict): A mapping of original nodes to their corresponding community nodes.
            - decompression_mapping (dict): A mapping of community nodes to the original nodes they represent.

    Raises:
        ValueError: If the specified method is not supported or if the graph is not connected and connectivity is required.

    Example:
        >>> import networkx as nx
        >>> from graphpack.compression import compress_graph, detect_communities, compress_graph_partition_based

        >>> # Example: Compress a weighted graph (gene network) using the Louvain method
        >>> G = nx.Graph()
        >>> edges = [("HIF1A", "EGFR", 0.934), ("HIF1A", "JAK2", 0.784), ("HIF1A", "IGF1R", 0.752), ("EGFR", "IGF1R", 0.989), ("JAK2", "IGF1R", 0.981)]
        >>> G.add_weighted_edges_from(edges)
        >>> compressed_graph, compression_mapping, decompression_mapping = compress_graph(G, method='louvain', is_weighted=True, is_gene_network=True)
        >>> print(compressed_graph.edges(data=True))
        [(1, 0, {'weight': 2.525})]
        >>> print(compression_mapping)
        {1: ['HIF1A', 'EGFR'], 0: ['JAK2', 'IGF1R']}
        >>> print(decompression_mapping)
        {'HIF1A': 1, 'EGFR': 1, 'JAK2': 0, 'IGF1R': 0}
    """
    if method in ['louvain', 'greedy', 'label_propagation', 'asyn_fluidc', 'spectral',
                  'hclust', 'node2vec', 'deepwalk', 'cpm', 'nmf']:
        # Detect communities
        try:
            with tqdm(desc="🕵️Community detection", total=1) as pbar:
                partition = detect_communities(graph, method=method, resolution=resolution, seed=seed, k=k)
                pbar.update(1)

        except ValueError as e:
            raise ValueError(f"Error detecting communities: {e}")

        # Compress network based on communities
        return compress_graph_partition_based(graph,
                                              partition,
                                              is_weighted=is_weighted,
                                              is_gene_network=is_gene_network)

    else:
        raise ValueError(f"Unsupported compression method: {method}")


def compute_and_save_edges_to_remove(graph, community_graph, output_folder, edges_to_remove_file='edges_to_remove.txt'):
    """
    Computes and saves edges to remove from the original network based on the compressed network.

    Args:
        graph (nx.Graph): The original network graph.
        community_graph (nx.Graph): The compressed network graph where nodes represent communities.
        output_folder (str): The folder path to save the output file.
        edges_to_remove_file (str): The name of the file to save the edges to remove (default is 'edges_to_remove.txt').

    Returns:
        None

    Raises:
        ValueError: If the output folder does not exist or is not a directory.
    """
    global FIG_NUM
    try:
        # Ensure that the output folder exists
        if not os.path.exists(output_folder) or not os.path.isdir(output_folder):
            raise ValueError("Output folder does not exist or is not a directory.")

        edges_to_remove = []

        # Reconstruct original network
        reconstructed_network = nx.Graph()

        # draw_graph(reconstructed_network, title=f"{FIG_NUM}_step_0")
        FIG_NUM += 1

        # Add nodes and edges from compressed network
        for community_id, nodes in community_graph.nodes(data='nodes'):
            reconstructed_network.add_nodes_from(nodes)
            # draw_graph(reconstructed_network, title=f"{FIG_NUM}_adding_nodes_from_community_{community_id}")
            FIG_NUM += 1

            # Add edges among all the nodes in the same supernode of the compressed network
            for node1 in nodes:
                for node2 in nodes:
                    if node1 != node2:
                        reconstructed_network.add_edge(node1, node2)

            # draw_graph(reconstructed_network, title=f"{FIG_NUM}_adding_edges_within_community_{community_id}")
            FIG_NUM += 1

        # Add edges between all the nodes in a supernode and all the nodes in a connected supernode
        for community_id, nodes in community_graph.nodes(data='nodes'):
            for neighbor_id in community_graph.neighbors(community_id):
                # print(f"Community {community_id} has neighbor {neighbor_id}")
                neighbor_nodes = community_graph.nodes[neighbor_id]['nodes']

                for node1 in nodes:
                    for node2 in neighbor_nodes:
                        reconstructed_network.add_edge(node1, node2)

                # draw_graph(reconstructed_network, title=f"{FIG_NUM}_adding_edges_between_community_{community_id}_and_{neighbor_id}")
                FIG_NUM += 1

        # Iterate through original network to identify edges to remove and add
        for edge in graph.edges(data=True):
            source = edge[0]
            target = edge[1]

            # If edge exists in reconstructed network, remove it
            if reconstructed_network.has_edge(source, target):
                reconstructed_network.remove_edge(source, target)

        # draw_graph(reconstructed_network, title=f"{FIG_NUM}_step_1")
        FIG_NUM = 100

        # Add remaining edges in the reconstructed network to the list of edges to remove
        for edge in reconstructed_network.edges(data=True):
            edges_to_remove.append((edge[0], edge[1]))

        # Save edge information to a file in edge list format
        with open(os.path.join(output_folder, edges_to_remove_file), 'w') as edge_list_file:
            for edge in edges_to_remove:
                edge_list_file.write(f"{edge[0]}\t{edge[1]}\n")

    except Exception as e:
        raise ValueError(f"Error computing and saving edges to remove: {e}")


def save_network_files(graph,
                       community_graph,
                       compression_mapping,
                       decompression_mapping,
                       output_folder,
                       compression_mapping_filename='compression_mapping',
                       decompression_mapping_filename='decompression_mapping',
                       labels_mapping=None,
                       save_data=False,
                       file_format='txt'):
    """
    Saves network files and mappings to the specified output folder.

    Args:
        graph (nx.Graph): The original network graph.
        community_graph (nx.Graph): The compressed network graph where nodes represent communities.
        compression_mapping (dict): A mapping of original nodes to their corresponding community nodes.
        decompression_mapping (dict): A mapping of community nodes to the original nodes they represent.
        output_folder (str): The folder path to save the output files.
        compression_mapping_filename (str): The base name for the compression mapping files (default is 'compression_mapping').
        decompression_mapping_filename (str): The base name for the decompression mapping files (default is 'decompression_mapping').
        labels_mapping (dict, optional): A mapping of node labels (default is None).
        save_data (bool, optional): Whether to save edge data (default is False).
        file_format (str, optional): The file format to save the network files (default is 'txt').

    Returns:
        None

    Raises:
        ValueError: If the output folder does not exist or cannot be created.
    """
    try:
        # Ensure that the output folder exists
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # Save original and compressed networks as edge lists
        save_graph(graph,
                   file_path=os.path.join(output_folder, 'original_network.txt'),
                   save_data=save_data)
        save_graph(community_graph,
                   file_path=os.path.join(output_folder, f'compressed_network.{file_format}'),
                   save_data=save_data)

        # Save compression and decompression mappings as MessagePack files
        with open(os.path.join(output_folder, f"{compression_mapping_filename}.msgpack"), 'wb') as file:
            file.write(msgpack.packb(compression_mapping))

        with open(os.path.join(output_folder, f"{decompression_mapping_filename}.msgpack"), 'wb') as file:
            file.write(msgpack.packb(decompression_mapping))

        # Save the compression_mapping dictionary to a JSON file
        with open(os.path.join(output_folder, f"{compression_mapping_filename}.json"), 'w') as json_file:
            json.dump(compression_mapping, json_file, indent=4)

        # Save the decompression_mapping to a JSON file
        with open(os.path.join(output_folder, f"{decompression_mapping_filename}.json"), 'w') as json_file:
            json.dump(decompression_mapping, json_file, indent=4)

        # Save the labels_mapping to a JSON file, if given
        if labels_mapping:
            with open(os.path.join(output_folder, 'labels_mapping.json'), 'w') as json_file:
                json.dump(labels_mapping, json_file, indent=4)

    except Exception as e:
        raise ValueError(f"Error saving network files: {e}")


def reconstruct_original_network(output_folder,
                                 compressed_network_file='compressed_network.txt',
                                 compression_mapping_file='compression_mapping.msgpack',
                                 edges_to_remove_file='edges_to_remove.txt',
                                 output_file='reconstructed_network.txt'):
    """
    Reconstructs the original network from compressed data and removes specified edges.

    Args:
        output_folder (str): The folder path to read input files from and save the output file.
        compressed_network_file (str): The file name of the compressed network (default is 'compressed_network.txt').
        compression_mapping_file (str): The file name of the compression mapping (default is 'compression_mapping.msgpack').
        edges_to_remove_file (str): The file name of the edges to remove (default is 'edges_to_remove.txt').
        output_file (str): The name of the file to save the reconstructed network (default is 'reconstructed_network.txt').

    Returns:
        None

    Raises:
        ValueError: If any of the required files are missing in the specified folder.
    """
    global FIG_NUM

    try:
        # Ensure that all required files exist
        for file_name in [compressed_network_file, compression_mapping_file, edges_to_remove_file]:
            if not os.path.exists(os.path.join(output_folder, file_name)):
                raise ValueError(f"File '{file_name}' not found in the output folder.")

        # Load compressed network from edge list file
        with open(os.path.join(output_folder, compressed_network_file), 'rb') as file:
            community_graph = nx.read_edgelist(file, delimiter='\t')

        # Load compression mapping from MessagePack file
        with open(os.path.join(output_folder, compression_mapping_file), 'rb') as file:
            compression_mapping = msgpack.unpack(file, raw=False, strict_map_key=False)

            # Convert all keys in compression_mapping to strings for consistency
            compression_mapping = {str(k): v for k, v in compression_mapping.items()}

        # Load edges to remove from the edge list file
        edges_to_remove = set()
        with open(os.path.join(output_folder, edges_to_remove_file), 'r') as edge_list_file:
            for line in edge_list_file:
                source, target = line.strip().split('\t')
                edges_to_remove.add((source, target))

        # Reconstruct original network
        reconstructed_network = nx.Graph()
        # draw_graph(reconstructed_network, title=f"{FIG_NUM}_reconstruction_step_0")
        FIG_NUM += 1

        # Add nodes and edges from compressed network
        for community_id, nodes in compression_mapping.items():
            # print(f"Nodes in community {community_id}: {nodes}")
            reconstructed_network.add_nodes_from(nodes)
            # draw_graph(reconstructed_network, title=f"{FIG_NUM}_reconstruction_adding_nodes_from_community_{community_id}")
            FIG_NUM += 1

            # Add edges among all the nodes in the same supernode of the compressed network
            for node1 in nodes:
                for node2 in nodes:
                    if node1 != node2:
                        reconstructed_network.add_edge(node1, node2)

            # draw_graph(reconstructed_network, title=f"{FIG_NUM}_reconstruction_adding_edges_within_community_{community_id}")
            FIG_NUM += 1

        # Add edges between all the nodes in a supernode and all the nodes in a connected supernode
        for community_id, nodes in compression_mapping.items():
            if community_id in community_graph.nodes:
                for neighbor_id in community_graph.neighbors(community_id):
                    neighbor_nodes = compression_mapping[str(neighbor_id)]

                    for node1 in nodes:
                        for node2 in neighbor_nodes:
                            reconstructed_network.add_edge(node1, node2)

                    # draw_graph(reconstructed_network, title=f"{FIG_NUM}_adding_edges_between_community_{community_id}_and_{neighbor_id}")
                    FIG_NUM += 1
            else:
                print(f"{ORANGE_BOLD}Warning: Community {community_id} not found in the compressed network.{RESET}")

        # Remove edges specified in edges_to_remove
        for edge_data in edges_to_remove:
            # Check if the edge exists before attempting to remove it
            if reconstructed_network.has_edge(edge_data[0], edge_data[1]):
                reconstructed_network.remove_edge(edge_data[0], edge_data[1])
            else:
                # If the edge doesn't exist, raise a warning
                print(f"{ORANGE_BOLD}Warning: Edge {edge_data[0]}-{edge_data[1]} does not exist in the graph.{RESET}")

        # draw_graph(reconstructed_network, title=f"{FIG_NUM}_reconstruction_step_1")
        FIG_NUM += 1

        # Save the reconstructed original network
        nx.write_edgelist(reconstructed_network, os.path.join(output_folder, output_file), delimiter='\t', data=False)

        return reconstructed_network

    except Exception as e:
        raise ValueError(f"Error reconstructing original network: {e}")


def calculate_compression_efficacy(output_folder,
                                   is_lossless=True,
                                   original_network_file='original_network.txt',
                                   compressed_network_file='compressed_network.txt',
                                   compression_mapping_file='compression_mapping.msgpack',
                                   decompression_mapping_file='decompression_mapping.msgpack',
                                   edges_to_remove_file='edges_to_remove.txt'):
    """
    Calculates and prints the efficacy of the network compression.

    Args:
        output_folder (str): The folder path containing the network files.
        is_lossless (bool, optional): Flag indicating whether the compression is lossless (default: True).
        original_network_file (str, optional): The file name of the original network (default: 'original_network.txt').
        compressed_network_file (str, optional): The file name of the compressed network (default: 'compressed_network.txt').
        compression_mapping_file (str, optional): The file name of the compression mapping (default: 'compression_mapping.msgpack').
        decompression_mapping_file (str, optional): The file name of the decompression mapping (default: 'decompression_mapping.msgpack').
        edges_to_remove_file (str, optional): The file name containing the edges to remove (default: 'edges_to_remove.txt').

    Returns:
        None

    Raises:
        ValueError: If any required file is missing or cannot be read.
    """
    try:
        input_file_path = os.path.join(output_folder, original_network_file)
        compressed_file_path = os.path.join(output_folder, compressed_network_file)
        edges_to_remove_file_path = os.path.join(output_folder, edges_to_remove_file)
        compression_mapping_file_path = os.path.join(output_folder, compression_mapping_file)
        decompression_mapping_file_path = os.path.join(output_folder, decompression_mapping_file)

        # Calculate sizes of files
        original_size = os.path.getsize(input_file_path)
        compressed_size = os.path.getsize(compressed_file_path)
        compression_mapping_size = os.path.getsize(compression_mapping_file_path)
        decompression_mapping_size = os.path.getsize(decompression_mapping_file_path)

        # Calculate the lossy compression size
        lossy_compression_size = compressed_size + min(compression_mapping_size, decompression_mapping_size)

        # Calculate compression percentages
        percentage_compression = ((original_size - compressed_size) / original_size) * 100
        percentage_lossy_compression = ((original_size - lossy_compression_size) / original_size) * 100

        # Print results
        print(f"\nOriginal size: {original_size} bytes")
        print(f"Compressed size: {compressed_size} bytes")
        print(f"Percentage compression: {percentage_compression:.2f}%")

        print(f"\nLossy compression size: {lossy_compression_size} bytes")
        print(f"Percentage lossy compression: {percentage_lossy_compression:.2f}%")

        # Calculate lossless compression size and percentage
        if is_lossless:
            edges_to_remove_size = os.path.getsize(edges_to_remove_file_path)
            lossless_compression_size = compressed_size + edges_to_remove_size + min(compression_mapping_size,
                                                                                     decompression_mapping_size)
            percentage_lossless_compression = ((original_size - lossless_compression_size) / original_size) * 100

            print(f"\nLossless compression size: {lossless_compression_size} bytes")
            print(f"Percentage lossless compression: {percentage_lossless_compression:.2f}%")

    except Exception as e:
        raise ValueError(f"Error calculating compression efficacy: {e}")


def remove_unnecessary_files(output_folder_method):
    """
    Removes unnecessary files from the output folder.

    Args:
        output_folder_method (str): The folder path containing the output files.

    Returns:
        None
    """
    try:
        # Remove unnecessary files
        for file_name in ['original_network.txt', 'compression_mapping.msgpack', 'decompression_mapping.msgpack', ]:
            file_path = os.path.join(output_folder_method, file_name)
            if os.path.exists(file_path):
                os.remove(file_path)

    except Exception as e:
        raise ValueError(f"Error removing unnecessary files: {e}")


def perform_compression(input, output='data/output', output_format='txt', method='louvain', resolution=1.25, k=3,
                        seed=123,
                        is_weighted=False, is_gene_network=False, is_lossless=False, plot=False, is_interactive=False,
                        plot_disconnected=False, separate_communities=False, title='', verbosity=2):
    """
    Suggested pipeline for the GraphPack tool.

    This function performs graph compression using various community detection methods. It reads an input graph file,
    compresses the graph, and saves the compressed version along with optional visualizations and statistical summaries.

    Args:
        input (str): Path to the input graph file. This is a required argument.
        output (str): Path to the output folder where results will be saved. Default is 'data/output'.
        output_format (str): File format to save the network files. Options: 'edgelist', 'txt' (default), 'csv', 'tsv', 'json', 'gpickle', 'gml', 'graphml', 'net', 'pajek', 'gexf', 'yaml', 'yml'.
        method (str): Community detection method to use for graph compression. Options: 'louvain' (default), 'greedy', 'label_propagation', 'asyn_fluidc', 'spectral', 'hclust', 'node2vec', 'deepwalk', 'cpm', 'nmf'.
        resolution (float): Resolution parameter for Louvain and greedy methods. Default is 1.25.
        k (int): Number of clusters for clustering methods. Only applicable for methods requiring a cluster count (e.g., 'asyn_fluidc', 'spectral', 'hclust', 'node2vec', 'deepwalk', 'nmf'). Default is 3.
        seed (int): Random seed for reproducibility. Default is 123.
        is_weighted (bool): Flag to indicate if the graph should consider edge weights. Use this flag if the graph is weighted. Default is False.
        is_gene_network (bool): Flag to assign biologically meaningful labels to communities. Use this flag for gene networks. Default is False.
        is_lossless (bool): Flag to perform lossless compression. Use this flag if lossless compression is required. Default is False.
        plot (bool): Flag to plot the original and compressed graphs. Default is False.
        is_interactive (bool): Flag to produce interactive plots in HTML format. Use this flag to enable interactive plots. Default is False.
        plot_disconnected (bool): Flag to plot all nodes in a disconnected graph, not just the largest connected component. Default is False.
        separate_communities (bool): Flag to separate communities in the graph plot. Default is False.
        title (str): Title for the graph plot. Default is ''.
        verbosity (int): Verbosity level for logging information (0: minimal, 1: moderate, 2: detailed). Default is 2.

    Returns:
        None

    Examples:
        >>> from graphpack.compression import *
        >>> import networkx as nx
        >>> G = nx.Graph()
        >>> G.add_edges_from([(1, 2), (2, 3), (3, 1), (4, 5), (5, 6), (6, 7), (5, 8), (4, 8)])
        >>> input_graph = 'simple_graph.txt'
        >>> save_graph(G, input_graph)
        >>> perform_compression(input_graph, output='results', method='greedy', plot=True, is_interactive=True, plot_disconnected=True, separate_communities=True, title='Greedy')
    """
    input_name = os.path.splitext(os.path.split(input)[-1])[0]
    output_folder = str(os.path.join(output, input_name))
    output_folder_method = str(os.path.join(output, input_name, method))

    # Print according to verbosity level
    if verbosity > 0:
        print(f"\n{'=' * 80}") if verbosity > 1 else None
        print(f"{BOLD}{TITLE:^80}{RESET}\n") if verbosity > 1 else None
        print(f"{BOLD}▶ Input graph file:{RESET}   {input}\t({input_name})")
        print(f"{BOLD}▶ Output folder:{RESET}      {output_folder}")
        print(f"{BOLD}▶ Compression method:{RESET} {method}\n")

        print("▶ ⚖️ The graph provided will be treated as a weighted graph") if is_weighted else None
        print("▶ 🧬 The graph provided will be treated as a gene network") if is_gene_network else None
        print("▶ 🔝 The compression method will be lossless") if is_lossless else None
        print(
            "▶ ➗  The plotting function is separating nodes belonging to different communities") if separate_communities and plot else None

        print(f"▶ 🔢 Running {method} with resolution {resolution}") if method in ['louvain', 'greedy'] else None
        print(f"▶ 📏 Running {method} with k = {k} minimum clique size") if method == 'cpm' else None

        if method in ['asyn_fluidc', 'spectral', 'hclust', 'node2vec', 'deepwalk']:
            print(f"▶ ℹ️ Running {method} with k = {k} clusters")

    # Read graph from input file
    print("\n📑 Reading input graph...") if verbosity > 0 else None
    graph = read_graph(input)

    # Disable interactive mode if the graph is too large
    if graph.number_of_nodes() > MAX_N_NODES_INTERACTIVE and plot:
        print(f"\n{ORANGE_BOLD}Warning: The graph is too large to be interactive. Interactive mode is disabled.{RESET}")
        print(
            f"{BOLD}You can enable interactive mode for bigger graphs by changing the MAX_N_NODES_INTERACTIVE constant in the script.{RESET}")
        print(f"{BOLD}The original graph will not be plotted.{RESET}\n")
        is_interactive_ = False
    else:
        is_interactive_ = is_interactive

    # Draw original graph
    if graph.number_of_nodes() < MAX_N_NODES_INTERACTIVE and plot:
        print("🖌️Drawing original graph...") if verbosity > 1 else None
        draw_graph(graph, title='Original graph',
                   file_path=output_folder,
                   is_interactive=is_interactive_,  # Draw the interactive plot, if requested and not overridden
                   plot_disconnected_components=plot_disconnected)

    # Compute k for community detection methods
    k_arg = k
    # For visualization purposes, no more than 30 supernodes are allowed
    if method in ['asyn_fluidc', 'spectral', 'hclust', 'node2vec', 'deepwalk']:
        k = min(graph.number_of_nodes(), k_arg, MAX_K) if plot else min(graph.number_of_nodes(), k_arg)
        if verbosity > 0 and k != k_arg:
            print(
                f"\n{ORANGE_BOLD}Warning: k cannot be higher than the number of nodes.{RESET}") if k > graph.number_of_nodes() else None
            print(
                f"\n{ORANGE_BOLD}Warning: for visualization purposes, the maximum number of cluster k is {MAX_K}.{RESET}") if k_arg > MAX_K and plot else None
            print(
                f"{BOLD}You can increase the maximum number of clusters by changing the MAX_K constant in the script.{RESET}") if k_arg > MAX_K and plot else None
            print(f"{BOLD}Running {method} with k = {k} clusters.{RESET}\n")

    else:
        k = k_arg

    # Set the number of components to the MAX_K constant, or provided k, if the graph is too large
    if method == 'nmf':
        if graph.number_of_nodes() > MAX_N_NODES_NMF:
            k = max(MAX_K, k)
            output_folder_method = str(os.path.join(output_folder, method) + f"_max_k_{k}")
            print(
                f"\n{ORANGE_BOLD}Warning: The graph is too large for automatic inference of the number of components for NMF.{RESET}")
            print(f"{BOLD}Running {method} with n_components = {k}{RESET}\n")
        else:
            k = 'auto'
            print(f"\n{ORANGE_BOLD}Warning: Automatic inference of the number of components for NMF.{RESET}\n")

    # Add the parameter to the output folder name
    if method in ['asyn_fluidc', 'spectral', 'hclust', 'node2vec', 'deepwalk', 'cpm']:
        output_folder_method = str(os.path.join(output_folder, method) + f"_k_{k}")

    if method in ['louvain', 'greedy']:
        output_folder_method = str(os.path.join(output_folder, method) + f"_r_{resolution}")

    # Compress network
    try:
        # Start the timeout countdown
        compressed_graph, compression_mapping, decompression_mapping = compress_graph(graph,
                                                                                      method=method,
                                                                                      seed=seed,
                                                                                      resolution=resolution,
                                                                                      k=k,
                                                                                      is_weighted=is_weighted,
                                                                                      is_gene_network=is_gene_network)

        if is_gene_network:
            labels = {community_id: str(label[0]) for community_id, label in compressed_graph.nodes(data='label')}
        else:
            labels = None

        if plot:
            # Draw compressed graph
            partition_colors = None
            node_colors = None
            color_map = None

            if is_gene_network:
                print("\n🧬 Drawing compressed graph with GSEA labels...") if verbosity > 0 else None
                partition_colors, node_colors, color_map = assign_community_colors(graph,
                                                                                   compressed_graph,
                                                                                   decompression_mapping,
                                                                                   labels)

                draw_graph(compressed_graph,
                           labels=labels,
                           title=f"Compressed graph (GSEA labels) {input_name} {title}",
                           file_path=output_folder_method,
                           is_interactive=is_interactive,  # Draw the interactive plot, if requested
                           plot_disconnected_components=plot_disconnected,
                           node_color=partition_colors)

                if graph.number_of_nodes() <= MAX_N_NODES:
                    print("🧬 Drawing original graph with partition colors (GSEA labels)...") if verbosity > 1 else None
                    draw_graph(graph,
                               title=(f"Original graph with partition colors (GSEA labels) {input_name} {title}" +
                                      (" - communities layout" if separate_communities else "")),
                               file_path=output_folder_method,
                               is_interactive=is_interactive_,
                               # Draw the interactive plot, if requested and not overridden
                               node_color=node_colors,
                               color_map=color_map,
                               plot_disconnected_components=plot_disconnected,
                               separate_communities=separate_communities)
                else:
                    print(f"\n{ORANGE_BOLD}Warning: The graph is too large to be plotted.{RESET}")
                    print(
                        f"{BOLD}You can enable plotting for bigger graphs by changing the MAX_N_NODES constant in the script.{RESET}\n")

            print("\n🖌️Drawing compressed graph...") if verbosity > 0 else None
            if node_colors is None and color_map is None:
                partition_colors, node_colors, color_map = assign_community_colors(graph,
                                                                                   compressed_graph,
                                                                                   decompression_mapping)

            draw_graph(compressed_graph,
                       title=f"Compressed graph {input_name} {title}",
                       file_path=output_folder_method,
                       is_interactive=is_interactive,  # Draw the interactive plot, if requested
                       plot_disconnected_components=plot_disconnected,
                       node_color=partition_colors,
                       color_map=color_map)

            if graph.number_of_nodes() <= MAX_N_NODES_INTERACTIVE and not is_gene_network:
                print("🖌️Drawing original graph with partition colors...") if verbosity > 1 else None
                draw_graph(graph,
                           title=(f"Original graph with partition colors {input_name} {title}" +
                                  (" - communities layout" if separate_communities else "")),
                           file_path=output_folder_method,
                           is_interactive=is_interactive,  # Draw the interactive plot, if requested
                           node_color=node_colors,
                           color_map=color_map,
                           plot_disconnected_components=plot_disconnected,
                           separate_communities=separate_communities)

        # Save network files
        print("\n💾 Saving network files...") if verbosity > 0 else None
        # Saving edge weights only if the original graph is weighted and the compression is not lossless
        save_network_files(graph,
                           compressed_graph,
                           compression_mapping,
                           decompression_mapping,
                           output_folder_method,
                           labels_mapping=labels,
                           save_data=is_weighted and not is_lossless,
                           file_format=output_format)

        # Compute and save edges to remove
        if is_lossless:
            print("\n🗑️ Computing and saving edges to remove...") if verbosity > 0 else None
            compute_and_save_edges_to_remove(graph, compressed_graph, output_folder_method)

            # Reconstruct the original network
            print("🛠️ Reconstructing the original network...") if verbosity > 0 else None
            reconstructed_network = reconstruct_original_network(
                output_folder_method,
                compressed_network_file=f'compressed_network.{output_format}'
            )

            # Draw reconstructed graph
            num_nodes = reconstructed_network.number_of_nodes()
            if num_nodes <= MAX_N_NODES_INTERACTIVE and plot:
                print("\n🖌️ Drawing reconstructed graph...") if verbosity > 1 else None
                draw_graph(reconstructed_network,
                           title="Reconstructed graph",
                           file_path=output_folder_method,
                           is_interactive=is_interactive_,  # Draw the interactive plot, if requested and not overridden
                           plot_disconnected_components=plot_disconnected)

        # Print compression efficacy and remove unnecessary files
        print("\n⚙️ Calculating compression efficacy...") if verbosity > 1 else None
        calculate_compression_efficacy(output_folder_method,
                                       is_lossless=is_lossless,
                                       compressed_network_file=f'compressed_network.{output_format}')

        remove_unnecessary_files(output_folder_method)

        print(f"\n{GREEN_BOLD}✅  Graph compression completed successfully!{RESET}") if verbosity > 0 else None


    except ValueError as e:
        print(f"Error: {e}")
        print(f"{RED_BOLD}❌  Graph compression failed!{RESET}")

    except TimeoutError as e:
        print(f"TimeoutError: {e}")
        print(f"{RED_BOLD}❌  Graph compression timed out!{RESET}")

    except Exception as e:
        print(f"Exception: {e}")
        print(f"{RED_BOLD}❌  Graph compression failed!{RESET}")

    print(f"\n{BOLD}▶ Finished processing {input}\n{RESET}") if verbosity > 0 else None


def parse_args():
    """
    Parse command-line arguments for graph compression.

    Command-line arguments:

    Args:
        --input (str): Path to the input graph file. This is a required argument.\n
        --output (str): Path to the output folder where results will be saved. Default is 'data/output'.\n
        --output-format (str): Output file format. Options: '.edgelist', '.txt', '.csv', '.tsv', '.json', '.gpickle', '.gml', '.graphml', '.net', '.pajek', '.gexf', '.yaml', '.yml'. Default is 'txt'.\n
        --method (str): Community detection method to use for graph compression. Options: 'louvain', 'greedy' (default), 'label_propagation', 'asyn_fluidc', 'spectral', 'hclust', 'node2vec', 'deepwalk', 'cpm', 'nmf'.\n
        --resolution (float): Resolution parameter for Louvain and greedy methods. Default is 1.25.\n
        --k (int): Number of clusters for clustering methods. Only applicable for methods requiring a cluster count (e.g., 'asyn_fluidc', 'spectral', 'hclust', 'node2vec', 'deepwalk', 'nmf'). Default is 3.\n
        --seed (int): Random seed for reproducibility. Default is 123.\n
        --is-weighted (bool): Flag to indicate if the graph should consider edge weights. Use this flag if the graph is weighted.\n
        --is-gene-network (bool): Flag to assign biologically meaningful labels to communities. Use this flag for gene networks.\n
        --is-lossless (bool): Flag to perform lossless compression. Use this flag if lossless compression is required.\n
        --plot (bool): Flag to plot the original and compressed graphs. Default is False.\n
        --is-interactive (bool): Flag to produce interactive plots in HTML format. Use this flag to enable interactive plots.\n
        --plot-disconnected (bool): Flag to plot all nodes in a disconnected graph, not just the largest connected component.\n
        --title (str): Title for the graph plot. Default is ''.\n
        --verbosity (int): Verbosity level for logging information (0: minimal, 1: moderate, 2: detailed). Default is 2.\n

    Returns:
        args (argparse.Namespace): Parsed command-line arguments.
    """
    from graphpack import __version__

    # Create the custom parser
    parser = CustomArgumentParser(
        description="Compress graphs using different community detection algorithms.",
        epilog="For more information, please refer to the documentation.",
        add_help=False
    )

    # Create groups for different sets of options
    help_group = parser.add_argument_group(f'{BOLD}Options{RESET}')
    input_output_group = parser.add_argument_group(f'{BOLD}Input/Output Options{RESET}')
    method_group = parser.add_argument_group(f'{BOLD}Method Options{RESET}')
    graph_options_group = parser.add_argument_group(f'{BOLD}Graph and Compression Options{RESET}')
    plotting_group = parser.add_argument_group(f'{BOLD}Plotting Options{RESET}')
    misc_group = parser.add_argument_group(f'{BOLD}Miscellaneous Options{RESET}')

    # Add arguments to the parser
    help_group.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit.')
    help_group.add_argument('--version', action='version', version=f'%(prog)s {__version__}',
                            help='Show the version number and exit.')

    # Input/Output options
    input_output_group.add_argument('-i', '--input', type=str, required=True,
                                    help='Input graph file path')
    input_output_group.add_argument('-o', '--output', type=str, default='data/output',
                                    help='Output folder path. Default is "data/output"')
    input_output_group.add_argument('-f', '--output-format', type=str, default='txt',
                                    help=f'Output file format. Options: {EXTENSIONS} Default is "txt"')

    # Method options
    method_group.add_argument('-m', '--method', type=str, default='greedy',
                              help='Community detection method. Default is "greedy"')
    method_group.add_argument('-r', '--resolution', type=float, default=1.25,
                              help='Resolution parameter for Louvain and greedy methods. Higher values lead to more communities. Default is 1.25')
    method_group.add_argument('-k', '--k', type=int, default=3,
                              help="Number of clusters for clustering methods. Relevant for methods requiring a predefined number of clusters. For 'cpm' is the size of the smallest clique. Default is 3.")

    # Graph options
    graph_options_group.add_argument('-w', '--is-weighted', action='store_true',
                                     help="Indicates that the graph is weighted. If not specified, the graph is considered unweighted.")
    graph_options_group.add_argument('-g', '--is-gene-network', action='store_true',
                                     help="Indicates that the input graph is a gene network. If not specified, the graph is considered a general network.")
    graph_options_group.add_argument('-l', '--is-lossless', action='store_true',
                                     help="Use lossless compression for the graph data. If not specified, lossy compression is used.")

    # Plotting options
    plotting_group.add_argument('-p', '--plot', action='store_true',
                                help='Plot the original and compressed graphs. Default is False.')
    plotting_group.add_argument('-x', '--is-interactive', action='store_true',
                                help="Produce and save the plots also as interactive html files. The static plots will still be saved as images.")
    plotting_group.add_argument('--plot-disconnected', action='store_true',
                                help="Plot all the nodes in a disconnected graph, not only the largest connected component.")
    plotting_group.add_argument('--separate-communities', action='store_true',
                                help="Enforces separation of the identified communities in the plots.")
    plotting_group.add_argument('--title', type=str, default='Original graph',
                                help='Title of the graph plots. Default is "Original graph"')

    # Miscellaneous options
    misc_group.add_argument('-s', '--seed', type=int, default=123,
                            help='Random seed for reproducibility')
    misc_group.add_argument('-v', '--verbosity', type=int, default=2,
                            help='Verbosity level (0: minimal, 1: moderate, 2: detailed)')

    return parser.parse_args()


def main():
    args = parse_args()
    perform_compression(**vars(args))

if __name__ == "__main__":
    main()
