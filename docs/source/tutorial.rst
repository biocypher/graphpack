Demo Pipeline
=============

.. admonition:: Prerequisites and Troubleshooting
    :class: hint

    Before proceeding with the tutorial, make sure to install GraphPack by following the installation instructions in
    the `Getting Started <getting-started>`_ section.

    If you encounter any issues or have questions while following the tutorial, refer to the `API Reference <modules>`_
    or the `Glossaries <graph-compression-glossary>`_ for detailed information on the functions and algorithms used in
    GraphPack.

Demo Overview
-------------

Here is a high-level overview of what the demo script does:

* **Input Selection**: The script can use a demo graph, a custom gene list, or data from a custom API URL.
  Depending on the input type, the script constructs the appropriate file path for saving the graph.

* **Graph Retrieval**: For a demo graph, the script uses predefined examples (STRING DB urls to retrieve graphs).
  For a custom gene list, it generates a graph based on the specified genes. For a custom API URL, it fetches data from
  the specified endpoint. The script checks if the graph is weighted and saves it to the specified file path.

* **Compression Pipeline**: The script runs the compression methods specified by the user. It compresses the graph using
  the selected methods and saves the compressed graph to the output directory. The script iterates through the specified
  compression methods. For each method, it constructs the command to run the compression, including various options like
  input/output paths, compression method, seed for randomness, and whether to plot the results.

  For details on the compression pipeline, refer to the `Compression Pipeline <compression-pipeline>`_ section.

Input Selection and Graph Retrieval
-----------------------------------

The first step in the demo pipeline involves selecting and retrieving the input graph. Depending on the provided inputs,
the script can handle different scenarios as described below.

* **Using Demo Graphs**: When a demo number is specified, the script retrieves a predefined example graph from the
  STRING database. It saves the graph to a file, which will be used as the input for the compression pipeline.

* **Using Custom Gene Lists**: If a custom gene list is provided, the script generates a graph based on the specified
  genes. The gene list can be specified directly or selected based on the demo number. Additional parameters can be
  passed to customize the API request.

* **Using Custom API URLs**: For more flexibility, users can specify a custom API URL to fetch the graph data. This
  allows the retrieval of graphs from STRING DB endpoint.

The fetched data is parsed and saved to a file for compression.

.. admonition:: STRING DB API Usage
    :class: note

    The STRING API provides a rich set of options for fetching graph data. Users can specify custom parameters to
    filter the data based on species, required score, and other criteria. The API request can be customized using the
    ``custom_params`` argument. Please consult `STRING API documentation <https://string-db.org/help/api/>`_ for more
    information.

    Example of request url for the STRING API:

    .. code-block:: python

        api_url = "https://string-db.org/api/json/interaction_partners?identifiers=P53_HUMAN%0dMDM2_HUMAN%0dATM_HUMAN&species=9606&required_score=900"

    Example of API request with parameters and gene list:

    .. code-block:: python

        gene_list = ['ATP1A1', 'BRCA1', 'CDK2', 'DNMT1', 'EGFR', 'FGFR2', 'GATA1', 'HIF1A', 'IGF1R', 'JAK2']
        api_url = "https://string-db.org/api/json/network"
        params = {
            "identifiers": "%0d".join(gene_list),
            "species": species,
            'caller_identity': 'my_app',
            'required_score': required_score,
            'network_flavor': network_flavour,
            'network_type': network_type
        }

    Example of custom parameters for the STRING API request (they will update the default parameters):

    .. code-block:: python

        custom_params = {"species": 9606, "required_score": 100}

Demo Options
------------

To demonstrate the compression pipeline, we provide a demo that showcases the compression of simple graphs using
various methods. The demo script allows users to specify the input graph, output directory, compression methods, and
other parameters. The demo script can be run using with some command line arguments, here detailed, or by using the
API directly, calling the main function :func:`~demo.run_demo` with the desired parameters.

.. code-block:: none

    -d <int>, --demo <int>
        Demo number to run. Accepts values from 0 to 6.

    -D, --use-gene-list
        Use a predefined gene list corresponding to the demo number (0-3).
        Overrides --gene-list if both are provided.

    -L <str> [<str> ...], --gene-list <str> [<str> ...]
        List of genes to use for the query. Ignored if --use-gene-list and --demo are specified.

    -c <str> [<str> ...], --custom-params <str> [<str> ...]
        Custom parameters to be sent with the API request. Only applicable when --use-gene-list is set.
        Should be in the format 'key1=value1 key2=value2'.

    -a <str>, --api-url <str>
        Custom API URL to fetch the data from. If not specified, the default STRING API URL is used.


Input/Output Options
~~~~~~~~~~~~~~~~~~~~

.. code-block:: none

    -i <str>, --input <str>
        Path to the input file containing the graph data. Default is '../data/input/simple_graph.txt'.
        If --demo is specified, this argument is ignored.

    -o <str>, --output <str>
        Path to the output folder where results will be saved. Default is '../data/output'.

    -f <str>, --output-format <str>
        Output file format. Options: 'edgelist', 'txt', 'csv', 'tsv', 'json', 'gpickle',
        'gml', 'graphml', 'net', 'pajek', 'gexf', 'yaml', 'yml'. Default is 'txt'.

Method Options
~~~~~~~~~~~~~~

With the demo script, it is possible to run multiple compression methods at once. The method(s) can be specified using
the command line argument ``--method``, followed by at least one of the available methods mentioned above, or ``all``
to apply all of them.

.. code-block:: none

    -m <str> [<str> ...], --method <str> [<str> ...]
        List of compression methods to run. If 'all', all the implemented compression methods will be used.
        Default is ['greedy']. Accepts multiple methods, among {METHODS}.

    -r <float>, --resolution <float>
        Resolution parameter for Louvain and greedy methods. Higher values lead to more communities.
        If not specified, the default value (1.25) is used.

    -k <int>, --k <int>
        Number of clusters for clustering methods. Relevant for methods requiring a predefined number of
        clusters. Default is 30. For 'cpm' is the size of the smallest clique. Default is 3.

.. admonition:: Disclaimer for Compression Methods
    :class: error

    The compression methods implemented are heuristic and may not always provide the optimal solution. Users are
    encouraged to experiment with different methods and parameters to find the best compression strategy for their
    specific use case.

Graph and Compression Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: none

    -w, --is-weighted
        Indicates that the graph is weighted. If not specified, the graph is considered unweighted.

    -g, --is-gene-network
        Indicates that the input graph is a gene network. If not specified, the graph is considered a
        general network.

    -l, --is-lossless
        Use lossless compression for the graph data. If not specified, lossy compression is used.

Plotting Options
~~~~~~~~~~~~~~~~

.. code-block:: none

    --no-plot
        Disable plotting the original and compressed graphs.

    --no-interactive-plot
        Disable saving the plots also as interactive html files. The static plots will still be saved
        as images.

    --plot-disconnected
        Plot all the nodes in a disconnected graph, not only the largest connected component.

    --separate-communities
        Enforces separation of the identified communities in the plots.

