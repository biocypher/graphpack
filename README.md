## Description
**GraphPack** is a Python tool engineered to facilitate the compression and visualization of large-scale networks, such
as protein-protein interaction (PPI) networks or metabolic pathways. It offers a user-friendly interface that enables
the application of diverse graph compression algorithms and the visualization of the original and compressed networks.

GraphPack provides flexibility in how you interact with it, supporting both command-line interface (CLI) usage with
arguments and integration into Python applications via an API.

The tool supports both weighted and unweighted graphs, allowing users to analyze a wide range of network types.
It is specifically designed to handle large-scale, biological networks  such as protein-protein interaction (PPI)
networks, gene regulatory networks,and metabolic pathways, but can be applied to any network data.

## Detailed Description
GraphPack includes a variety of graph compression algorithms to choose from:

- **Louvain Clustering**
- **Greedy Algorithm**
- **Label Propagation**
- **Asynchronous Fluid Communities**
- **Spectral Clustering**
- **Hierarchical Clustering**
- **Node2Vec**
- **DeepWalk**
- **Clique Percolation Method (CPM)**
- **Non-negative Matrix Factorization (NMF)**

GraphPack generates mapping files that maintain the relationship between the original and compressed nodes. This
ensures that the compressed network can be decompressed to the original network without any loss of information, in case
of lossless compression, or in any case that information about the relationship between the new nodes and the old nodes
is available.

GraphPack provides robust visualization options for both the original and compressed networks, facilitating 
easy comparison and in-depth analysis.

## Installation
Install GraphPack from PyPI via:

```bash
pip install graphpack
```

## License
This project is licensed under the MIT License.
