Overview
========

GraphPack provides flexibility in how you interact with it, supporting both command-line interface (CLI) usage with
arguments and integration into Python applications via an API.

To have a better understanding of the compression pipeline, please refer to the :ref:`Compression Pipeline`. Here
the focus is on the main features of the tool, and how to use them.

Input and Output
----------------

GraphPack requires as input a file containing a network, which can be specified as an edge list in various formats.
The tool then compresses the network and saves the compressed network, in the chosen format, as well as additional
files (see :ref:`Mapping and Annotation`) and plots (see :ref:`Visualization`), in the specified output directory.

.. admonition:: Input File, Output Directory and Output Format
        :class: note

        The input file can be specified using the command line argument ``--input <str>`` or, if you are using the
        Python API, as the first parameter (``input='path/to/graph.edgelist'``) of the
        main :func:`~compression.perform_compression` function.

        The output directory can be specified using the command line argument ``--output <str>`` or by setting the
        parameter ``optput='path/to/output_folder'``, if using the API. The output file format can be specified as well,
        via the command line argument ``--output-format <str>`` (or ``output_format='txt'``).

The input file should contain the network data. The nodes can be specified as strings or integers, and the edges can
be weighted or unweighted. For example, the following files represent the same toy network, in different formats (some
of which are weighted).

.. admonition:: Input/Output File Formats
    :class: hint

    The input network can be provided as an edge list, in various formats (``json``, ``edgelist``, ``txt``, ``csv``,
    ``tsv``, ``gpickle``, ``gml``, ``graphml``, ``net``, ``pajek``, ``gexf``, ``gml``, ``yaml``, ``yml``).

    .. code-block:: bash
        :caption: Example of different input file formats

            # Contents of 'graph.edgelist':
            A B
            B C

            # Contents of 'graph.txt':
            A B {'weight': 0.50}
            B C {'weight': 0.99}

            # Contents of 'graph.csv':
            source,target
            A,B
            B,C

            # Contents of 'graph.tsv':
            source\ttarget\tweight
            A\tB\t0.50
            B\tC\t0.99

            # Contents of 'graph.json':
            {
                "nodes": ["A", "B", "C"],
                "edges": [
                    {"source": "A", "target": "B", "weight": 0.50},
                    {"source": "B", "target": "C", "weight": 0.99}
                ]
            }


The tool supports both weighted and unweighted graphs, allowing users to analyze a wide range of network types.

It is specifically designed to handle large-scale, biological networks  such as protein-protein interaction (PPI)
networks, gene regulatory networks,and metabolic pathways, but can be applied to any network data. In case of a gene
network, GraphPack also generates a labels mapping file that maintains the relationship between the compressed node
IDs (numeric) and their associated labels (GSEA top-k most enriched terms), which is also saved in the output directory.

The compression can be lossy or lossless. In case of lossless compression, the tool ensures that the compressed
network can be decompressed to the original network without any loss of information. To achieve this, in accordance
to the :term:`Minimum Description Length` (MDL) principle, the tool generates a supplementary file that contains the
edges to be removed from the reconstructed network from the compressed one, in order to obtain the original network.
This file s saved in the output directory.

.. admonition:: Graph and Compression option flags
        :class: note

        To specify that the graph is weighted, use the ``--is-weighted`` flag (or ``is_weighted=True`` in the API).
        To specify that the graph is a gene network, use the ``--is-gene-network`` flag (or ``is_gene_network=True`` in
        the API). To use lossless compression, use the ``--is-lossless`` flag (or ``is_lossless=True`` in the API).

Compression Algorithms
----------------------

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

It is possible to specify the number of clusters to form, the minimum size of the cliques, or the number of dimensions
for the embedding methods (parameter ``k``). For Louvain and Greedy algorithms, the resolution parameter can be also
changed. Higher values lead to more communities. The tool supports the specification of the seed for reproducibility.

.. admonition:: Method Options
        :class: note

        The method can be specified using the command line argument ``--method <str>`` (or ``method='greedy'`` in
        the API). To change the default parameters of the compression methods, use ``--resolution <float>``
        (``resolution=1.25`` in the API) and ``--k <int>`` (``k=3`` in the API). To set a different seed for
        reproducibility, ``--seed <int>`` (``seed=123``).

For further details on each compression method, refer to the :ref:`Compression Algorithm Glossary`.

Mapping and Annotation
----------------------

GraphPack generates mapping files that maintain the relationship between the original and compressed nodes. This
ensures that the compressed network can be decompressed to the original network without any loss of information, in case
of lossless compression, or in any case that information about the relationship between the new nodes and the old nodes
is available. For example, it is possible to know that the compressed node ``2`` represents the original nodes
``"MAPK1"``, ``"MAP2K1"``, ``"RPS6KA1"``, etc. The mapping files are saved in the output directory.

.. raw:: html

   <style>
   .json-table {
       width: 100%;
       border-collapse: collapse;
       margin-bottom: 20px;
   }
   .json-table th, .json-table td {
       padding: 10px;
       text-align: left;
       vertical-align: top; /* Ensure vertical alignment at the top */
   }
   pre {
       background-color: #f0f0f0; /* Gray background for JSON blocks */
       padding: 10px;
       border-radius: 5px;
       overflow-x: auto; /* Enable horizontal scroll if needed */
       margin: 0; /* Remove default margin */
   }
   .key {
       color: green; /* Green color for keys */
       font-weight: bold; /* Bold font for keys */
   }
   .value {
       color: #B22222; /* Lighter bordeaux red color for values */
   }
   </style>

   <table class="json-table">
   <tr>
   <th>compression_mapping.json</th>
   <th>decompression_mapping.json</th>
   </tr>
   <tr>
   <td><pre>
   {
       <span class="key">"2"</span>: [
           <span class="value">"MAPK1"</span>,
           <span class="value">"MAP2K1"</span>,
           <span class="value">"RPS6KA1"</span>,
           ...
       ],
       <span class="key">"1"</span>: [
           <span class="value">"TP53"</span>,
           <span class="value">"SFN"</span>,
           <span class="value">"HIF1A"</span>,
           ...
       ],
       ...
   }
   </pre></td>
   <td><pre>
   {
       <span class="key">"MAPK1"</span>: <span class="value">2</span>,
       <span class="key">"MAP2K1"</span>: <span class="value">2</span>,
       <span class="key">"RPS6KA1"</span>: <span class="value">2</span>,
       ...,
       <span class="key">"TP53"</span>: <span class="value">1</span>,
       <span class="key">"SFN"</span>: <span class="value">1</span>,
       <span class="key">"HIF1A"</span>: <span class="value">1</span>,
       ...,
       <span class="key">"STAT3"</span>: <span class="value">4</span>,
       <span class="key">"JAK1"</span>: <span class="value">4</span>,
       <span class="key">"JAK2"</span>: <span class="value">4</span>,
       ...
   }

   </pre></td>
   </tr>
   </table>

If the graph is a gene network, the tool also generates a labels mapping file, that maintains the relationship between
the compressed nodes and associated labels (GSEA top-k most enriched terms), which is saved in the output directory.
This file is useful for interpreting the compressed network, as it provides information about the biological functions
associated with the supernodes (cluster of genes/nodes of the original network), for example that cluster ``2`` is
likely to represent the ``"MAPK signaling pathway"`` (i.e. the genes in compressed node ``2`` have that term enriched).

.. raw:: html

   <p><strong>labels_mapping.json</strong></p>

   <pre>
   {
       <span class="key">"2"</span>: <span class="value">"MAPK signaling pathway"</span>,
       <span class="key">"1"</span>: <span class="value">"p53 signaling pathway"</span>,
       <span class="key">"4"</span>: <span class="value">"JAK-STAT signaling pathway"</span>,
       <span class="key">"8"</span>: <span class="value">"Glioma"</span>,
       <span class="key">"6"</span>: <span class="value">"Wnt signaling pathway"</span>,
       <span class="key">"7"</span>: <span class="value">"C-type lectin receptor signaling pathway"</span>,
       <span class="key">"3"</span>: <span class="value">"PI3K-Akt signaling pathway"</span>,
       <span class="key">"9"</span>: <span class="value">"Hippo signaling pathway"</span>,
       <span class="key">"10"</span>: <span class="value">"Apoptosis"</span>,
       <span class="key">"0"</span>: <span class="value">"TGF-beta signaling pathway"</span>,
       <span class="key">"5"</span>: <span class="value">"Notch signaling pathway"</span>
   }
   </pre>

Compression Pipeline
====================

Here you can see the steps involved in the (standard) compression pipeline, as implemented in the
:func:`~compression.perform_compression` function:

.. figure:: _images/PerformCompression.png
    :alt: Compression Pipeline
    :align: center

