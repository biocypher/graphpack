Graph Compression Glossary
==========================

.. glossary::

    Network
        A network is a collection of entities (such as proteins or genes) and relationships between the entities
        (such as protein-protein interactions, gene co-expression, transcriptional regulation, or similarities).

        Networks are used to model complex systems in a wide range of fields, including biology, social science,
        and computer science. Their structure can be visualized as a :term:`graph`.

    Graph
        A graph is a mathematical representation of a network, consisting of nodes (vertices) and edges (links) that
        connect the nodes. Graphs are used to model any complex system that can be represented as a network.

    Graph Compression
        Graph compression (:cite:t:`compression`) is intended as the process of reducing the size of a network while
        preserving its :term:`structural properties`. This is achieved by identifying and merging similar nodes
        (:term:`community`) and/or edges.

        One of the main goals of graph compression is to simplify the network structure, making it easier to analyze
        and visualize. This process is particularly useful for large-scale networks, where the sheer number of nodes
        and edges can make it difficult to interpret the underlying relationships.

    Lossless Compression
        Lossless compression is a type of data compression that allows the original data to be perfectly reconstructed
        from the compressed data. This means that no information is lost during the compression process.

        Lossless compression algorithms are commonly used in graph compression to reduce the size of a network without
        losing any of the underlying structure. It is commonly achieved by producing, in addition to the compressed
        network, a set of instructions that can be used to reconstruct the original network (for example, a mapping of
        the nodes in the compressed network to the nodes in the original network, and a list of edges to add to/remove
        from the reconstructed network; :cite:t:`mdl`).

        .. note::
            Lossless compression algorithms are typically more computationally intensive than lossy compression
            algorithms, as they must preserve all of the original data, especially in the case of large-scale networks.

    Minimum Description Length
        The minimum description length (MDL) principle is a method for selecting the best model from a set of
        competing models. The MDL principle states that the best model is the one that minimizes the total length
        of the model description and the data needed to reconstruct the original network.

        In the context of graph compression, the MDL principle can be used to select the best compression algorithm
        for a given network. By comparing the total length of the compressed network and the data needed to reconstruct
        the original network, it is possible to identify the most effective compression algorithm.

    Community
        A community is a group of nodes that are more densely connected to each other than to the rest of the network.
        Community is a key concept in graph theory and network analysis, as it serves as a proxy for potentially
        relevant functional modules of the original network, such as protein complexes or metabolic pathways.

    Community Detection
        Community detection is the process of identifying clusters of nodes in a network.

        There are many different community detection algorithms available. Some algorithms are based on the structure
        of the network (for example, the density of connections between nodes), while others use statistical methods
        to identify communities (such as modularity optimization).

    Compression Algorithm
        A compression algorithm is a set of rules and procedures used to reduce the size of a network. These algorithms
        typically identify and merge similar nodes and/or edges, resulting in a simplified network structure that tries
        to retain the essential relationships between nodes.

        There are many different compression algorithms available, each with its own strengths and weaknesses. Some
        algorithms are designed to work with specific types of networks, such as social networks or biological networks,
        while others are more general-purpose and can be applied to a wide range of network types. For further details,
        see :ref:`compression algorithm glossary`.

    Node Embedding
        Node embedding is a technique used to represent nodes in a network as vectors in a low-dimensional space. These
        vector representations capture the topological properties of the network, allowing for efficient analysis and
        visualization. Node embedding is often used in conjunction with compression algorithms to simplify the network
        structure while preserving the essential relationships between nodes.

    Random Walk
        A random walk is a mathematical process that describes a path through a network in which each step is chosen
        randomly. Random walks are often used to explore the structure of a network and identify important nodes and
        edges. Random walks can be used in compression algorithms to generate node embeddings.

    Structural Properties
        The structural properties of a network refer to the arrangement of nodes and edges within the network. These
        properties include the degree distribution, clustering coefficient, degree centrality, and modularity.

        The structural properties of a network are important for understanding its behavior and function. By analyzing
        these properties, it is possible to identify key nodes and edges, detect communities, and predict how the network
        will evolve over time.

    NetworkX
        NetworkX is a Python library for the creation, manipulation, and analysis of complex
        networks. It provides tools for generating graphs, computing network properties, and visualizing network
        structures.

        In NetworkX, graphs are represented as collections of nodes and edges, with each node and edge capable of
        storing arbitrary data as Python dictionaries. This flexibility makes NetworkX a powerful tool for working
        with a wide range of network types and structures.

        For more information, see the `NetworkX Documentation <https://networkx.org/documentation/stable/index.html>`_.

    Matplotlib
        Matplotlib is a Python library for creating static, animated, and interactive visualizations in Python. It
        provides tools for creating a wide range of plots, including line plots, scatter plots, bar plots, and
        histograms.

        Matplotlib is designed to work seamlessly with NumPy, Pandas, and other Python libraries, making it easy to
        create complex visualizations from data stored in these formats.

        For more information, see the `Matplotlib Documentation <https://matplotlib.org/stable/contents.html>`_.

    PyVis
        PyVis is a Python library for creating interactive network visualizations in Jupyter notebooks. It provides
        tools for creating dynamic visualizations of complex networks, including interactive graphs, node-link diagrams,
        and tree maps.

        PyVis is built on top of the :term:`VisJS` JavaScript library, which provides a wide range of interactive features
        for visualizing networks. PyVis makes it easy to create interactive visualizations of networks without having
        to write complex JavaScript code.

        For more information, see the `PyVis Documentation <https://pyvis.readthedocs.io/en/latest/>`_.

    VisJS
        VisJS is a JavaScript library for creating interactive network visualizations in web browsers. It provides
        tools for creating dynamic visualizations of complex networks, including interactive graphs, node-link diagrams,
        and tree maps.

        VisJS is designed to work seamlessly with Python libraries such as PyVis, making it easy to create interactive
        visualizations of networks in Jupyter notebooks. VisJS provides a wide range of interactive features, including
        zooming, panning, and highlighting nodes and edges.

        For more information, see the `VisJS Documentation <https://visjs.org/>`_.

    Plotly
        Plotly is a Python library for creating interactive visualizations in Python. It provides tools for creating
        a wide range of interactive plots, including line plots, scatter plots, bar plots, and 3D plots.

        Plotly is designed to work seamlessly with Pandas, NumPy, and other Python libraries, making it easy to create
        complex visualizations from data stored in these formats.

        For more information, see the `Plotly Documentation <https://plotly.com/python/>`_.

    Plotly Sankey Diagram
        The Plotly Sankey Graph Object is a specialized visualization tool for displaying flow diagrams. It is
        particularly useful for visualizing the flow of resources, energy, or information through a network.

        The Sankey Graph Object provides tools for creating interactive Sankey diagrams, including the ability to
        customize the appearance of nodes and edges, add annotations, and create dynamic visualizations.

        For more information, see the `Plotly Sankey Graph Object Documentation <https://plotly.com/python/sankey-diagram/>`_.
