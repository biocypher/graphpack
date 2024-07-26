Compression Algorithm Glossary
==============================

.. glossary::

   Louvain Clustering
        Louvain Community Detection algorithm (:cite:t:`louvain`) is used in compressing large-scale gene networks by
        identifying densely connected subgroups of nodes. This is a heuristic method based on modularity optimization.

        .. toggle::

            .. figure:: _images/1.png
                :width: 800px
                :alt: Louvain Clustering

            .. admonition:: Parameters and Hints

                 ``resolution``: Controls the size of the communities. Higher values lead to smaller communities.

                 ✅ Efficient for large networks.

                 ✅ Provides high modularity scores.

                 ⚠️ May result in overly large communities if the resolution parameter is not well-tuned.


   Greedy Algorithm
        This method uses Clauset-Newman-Moore greedy modularity maximization (:cite:t:`greedy`) to find the community
        partition with the largest modularity. It iteratively merges nodes to minimize the graph size while preserving
        its structure.

        .. toggle::

            .. figure:: _images/2.png
                :width: 800px
                :alt: Greedy Algorithm

            .. admonition:: Parameters and Hints

                 ``resolution``: Controls the size of the communities. Higher values lead to smaller communities.

                 ✅ Simple and fast implementation.

                 ✅ Can be applied to weighted and unweighted graphs.

                 ⚠️ May not provide the optimal solution compared to more sophisticated methods.


   Label Propagation
        Label Propagation (:cite:t:`label_propagation`) can compress networks by propagating labels across the network
        to identify and merge similar node communities. The algorithm is probabilistic and the found communities may
        vary in different executions.

        .. toggle::

            .. figure:: _images/3.png
                :width: 800px
                :alt: Label Propagation

            .. admonition:: Hints

                 ✅ Fast and scalable.

                 ⚠️ Can produce different results on different runs due to its stochastic nature.


   Asynchronous Fluid Communities
        The asynchronous fluid communities algorithm (:cite:t:`asyn_fluidc`) is based on the simple idea of fluids
        interacting in an environment, expanding and pushing each other. Its initialization is random, so found
        communities may vary on different executions.

        .. toggle::

            .. figure:: _images/4.png
                :width: 800px
                :alt: Asynchronous Fluid Communities

            .. admonition:: Parameters and Hints

                 ``k``: Number of clusters to form.

                 ✅ Suitable for dynamic and evolving networks.

                 ✅ Captures fluid and overlapping communities.

                 ⚠️ Requires the graph to be connected.

                 ⚠️ May produce different results on different runs due to random initialization.


   Spectral Clustering
        Spectral clustering (:cite:t:`spectral`), applied to community detection task by using the adjacency matrix as
        affinity, utilizes graph Laplacian eigenstructure to partition the graph into clusters. This method is effective
        for identifying complex, non-convex cluster shapes like nested circles on a 2D plane. By specifying
        ``affinity='precomputed'``, and providing the adjacency matrix as input, it is possible to accurately identify
        densely connected subgraphs (communities) within the network.

        .. toggle::

            .. figure:: _images/5.png
                :width: 800px
                :alt: Spectral Clustering

            .. admonition:: Parameters and Hints

                 ``k``: Number of clusters to form.

                 ✅ Effective for small to medium-sized networks.

                 ✅ Captures non-linear structures in the data.

                 ⚠️ Computationally expensive for large networks.


   Hierarchical Clustering
        Agglomerative clustering, applied to community detection tasks by utilizing the adjacency
        matrix as the metric, merges the most similar communities iteratively. This method is effective for identifying
        hierarchical structures within networks, where nodes are progressively merged based on their pairwise distances.

        .. toggle::

            .. figure:: _images/6.png
                :width: 800px
                :alt: Hierarchical Clustering

            .. admonition:: Parameters and Hints

                 ``k``: Number of clusters to form.

                 ✅ Provides a multi-level representation of the network.

                 ✅ Useful for visualizing the hierarchical structure.

                 ⚠️ May be computationally intensive for large datasets.


   Node2Vec
        Node2Vec (:cite:t:`node2vec`), applied to community detection tasks, embeds nodes through biased random walks to
        capture complex relationships within the network. These embeddings can then be clustered to identify communities.

        .. toggle::

            .. figure:: _images/7.png
                :width: 800px
                :alt: Node2Vec

            .. admonition:: Parameters and Hints

                 ``k``: Number of clusters to form.

                 ✅ Captures complex relationships within the network.

                 ✅️ Generates high-quality embeddings through random walks.

                 ⚠️ Requires careful tuning of hyperparameters for optimal performance.



   DeepWalk
        DeepWalk (:cite:t:`deepwalk`), applied to community detection tasks, learns node representations via truncated
        random walks. These representations are effective for large networks, preserving both local and global network
        structures.

        .. toggle::

            .. figure:: _images/8.png
                :width: 800px
                :alt: DeepWalk

            .. admonition:: Parameters and Hints

                 ``k``: Number of clusters to form.

                 ✅ Effective for large networks.

                 ✅ Preserves local and global network structures.

                 ⚠️ May require significant computational resources for large-scale applications.


   Clique Percolation Method (CPM)
        Clique Percolation Method (CPM) (:cite:t:`cpm`) is a community detection algorithm designed to identify
        overlapping communities by finding k-cliques that share (k-1) nodes. This method is particularly useful for
        detecting functionally significant modules within networks.

        .. toggle::

            .. figure:: _images/9.png
                :width: 800px
                :alt: Clique Percolation Method (CPM)

            .. admonition:: Parameters and Hints

                 ``k``: Size of the smallest clique.

                 ✅ Captures overlapping community structures.

                 ✅ Useful for detecting functionally significant modules.

                 ⚠️ May not scale well with very large networks.


   Non-negative Matrix Factorization (NMF)
        Non-negative Matrix Factorization (NMF) (:cite:t:`nmf`) is a dimensionality reduction technique that factorizes
        the adjacency matrix of a graph to detect communities. This method effectively preserves community structure in
        a reduced space, making it useful for identifying clusters within networks.

        .. toggle::

            .. figure:: _images/10.png
                :width: 800px
                :alt: Non-negative Matrix Factorization (NMF)

            .. admonition:: Parameters and Hints

                 ``k``: Number of components to factorize. If ``'auto'``, the algorithm will determine the optimal number.

                 ✅ Effective for dimensionality reduction.

                 ✅ Preserves community structure in the reduced space.

                 ⚠️ May require multiple runs to achieve stable results.
