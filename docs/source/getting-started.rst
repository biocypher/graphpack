Introduction
============

.. figure:: _images/graph.png
    :width: 800px
    :alt: Adapted from Ahnert, S. Generalised power graph compression reveals dominant relationship patterns in complex networks. Sci Rep 4, 4385 (2014). https://doi.org/10.1038/srep04385

    Figure adapted from :cite:t:`power_compression`.

**GraphPack** is a Python package that provides a comprehensive set of tools for compressing and analyzing large-scale
networks. The tool offers a variety of graph compression algorithms, allowing users to reduce the size of a network
while preserving its structural properties. GraphPack also includes visualization options for both the original and
compressed networks, enabling users to gain deeper insights into network structures.

Installation
============

GraphPack can be installed using pip, the Python package installer. The tool requires Python 3.6 or higher.

.. code-block:: bash

    pip install graphpack

To verify the installation, run the following commands:

.. code-block:: bash

    graphpack --version
    graphpack --help

The second command will display the help message, which provides an overview of the tool's functionality and available
options.

Quick Start
===========

To compress a network (assuming it is saved in edgelist format in file ``graph.txt``) using the Greedy algorithm,
run the following command in the terminal:

.. code-block:: bash

    graphpack --input path/to/graph.txt --output result --method greedy --plot

This command will compress the network stored in the "graph.edgelist" file using the Greedy algorithm and save the
resulting compressed network in the "result" directory, together with the mapping files between the original and
compressed nodes. Some plots will also be generated to visualize the original and compressed networks.

For more detailed instructions on using GraphPack, refer to the :ref:`overview`.

On Ramp
=======

.. admonition:: Don't know where to find network data?
    :class: note

    You can find various network datasets on the `Network Repository <https://networkrepository.com/>`_ or use a
    `NetworkX graph generator <https://networkx.org/documentation/stable/reference/generators.html>`_ to generate
    synthetic networks. For biological networks, you can explore databases such as `STRING <https://string-db.org/>`_
    or `BioGRID <https://thebiogrid.org/>`_.

For convenience, GraphPack provides a demo script that showcases the compression of simple graphs using various methods.
The demo includes the graph download from STRING API, data preparation, compression, and visualization. It is installed
along with the package and can be run using the following command:

.. code-block:: bash

    gp-demo --help

This command will display the help message for the demo script, providing information on how to run the demo and
specify the input graph, output directory, compression methods, and visualization options.

Example Usage
-------------

To run a basic demo, you can use a command like this:

.. code-block:: bash

    gp-demo --demo 0 --use-gene-list --input ./ --output demo_results --method greedy

This command will run the demo with the first example graph, using the Greedy algorithm for compression and saving the
results in the "demo_results" directory. In particular, the ``--use-gene-list`` flag indicates that the demo will use
a predefined gene list for the graph, e.g. the one below:

.. code-block:: python

    gene_list = ['ATP1A1', 'BRCA1', 'CDK2', 'DNMT1', 'EGFR', 'FGFR2', 'GATA1', 'HIF1A', 'IGF1R', 'JAK2']

From this list of genes, the demo will retrieve the corresponding network from the STRING API, compress it using the
Greedy algorithm, and save the results in the specified output directory. To do so, it simply saves the to-be-compressed
graph and calls the main compression function with the specified parameters.

Example Output
--------------

The output of the demo script includes the following files:

* ``demo_gene_list_0.txt``: the original graph in edgelist format, as retrieved from the STRING API

* ``compressed_network.txt``: the compressed graph in edgelist format

* ``compression_mapping.json``: the mapping between the original and compressed nodes

* ``decompression_mapping.json``: the mapping between the compressed and decompressed nodes

And the following plots:

* ``original_network.pdf``: visualization of the original network (also as interactive HTML file)

* ``compressed_graph.pdf``: visualization of the compressed network (also as HTML file)

* ``original_graph_with_partition_colors.pdf``: visualization of the original network with  partition colors (also as
  HTML file)

The names of the plots are automatically set to keep track uniquely of the method and original graph.
For more information on the demo script and its options, refer to the :ref:`Demo Pipeline` section.

.. admonition:: Want to learn more about graph compression?
    :class: hint

    Check out the `Graph Compression Glossary <graph-compression-glossary>`_ for definitions of key terms and concepts
    related to graph compression. The glossary provides a comprehensive overview of the terminology used in the field
    of graph compression, specifically in the context of GraphPack tool.
