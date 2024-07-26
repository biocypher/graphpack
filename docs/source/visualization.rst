Visualization
=============

GraphPack provides robust visualization options for both the original and compressed networks. This facilitates easy
comparison and in-depth analysis, helping users gain deeper insights into network structures.

The compression pipeline includes visualization of the original and compressed networks, as well as the original
graph colored by the communities found by the compression method. The visualization can be static or interactive,
depending on the user's preference and on the graph size.

Interactive visualizations are generated in HTML format, which can be opened in any web browser. The static
plots are saved as PDF images.

.. admonition:: Visualization Options
        :class: note

        * To generate plots, use the ``--plot``/``plot`` flag.

        * To enable the interactive visualization, use the ``--is-interactive``/``is_interactive`` flag.

        * The user can also specify a title for the plots using the ``--title <str>`` command line argument (or
          ``title='My Graph'`` in the API).

.. admonition:: Disconnect Graphs
        :class: warning

        In case of a disconnected graph, the default choice is to plot only the largest connected component, but it is
        possible to plot all nodes using the ``--plot-disconnected``/``plot_disconnected`` flag.

Static visualization Example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is an example of a static visualization of a compressed biological network, with GSEA labels as node names.

.. image:: _images/static_graph.png
    :width: 800px

It is possible to open the `compressed_graph with GSEA labels <compressed_graph_gsea_labels.pdf>`_ as a pdf file.

Please note that if the `original graph <original_graph.pdf>`_  was not considered as a gene network, the
`compressed graph <compressed_graph.pdf>`_ with the Greedy method, resolution 1.25 would have been the same, but
with only node IDs.

It is also possible to plot the
`original graph, with the nodes colored by their communities <original_graph_with_partition_colors.pdf>`_,
as found by the compression method.

Interactive visualization Example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The same network can be visualized interactively, allowing users to explore the network structure and relationships.

.. raw:: html
   :file: _files/interactive_graph.html

For an interactive visualization, with some additional features, open the html of this
`interactive compressed graph <compressed_graph_gsea_labels.html>`_. The interactive visualization of the
`original network, with nodes color-coded by community <original_graph_with_partition_colors.html>`_ is also available.

.. admonition:: Large Networks
        :class: error

        The interactive visualization is recommended for small to medium-sized graphs, as it may be slow for large
        networks. In such cases, the static visualization is recommended. It is anyway possible to save the interactive
        visualization of the compressed graph, because with a suitable compression strength any large network can be
        compressed to a way smaller one.

        For large networks the plotting is time-consuming and mostly ineffective, also in the static case. In such
        cases, it is recommended to disable the plotting by setting ``plot=False``.


**Interactive Visualization Features**

.. |zoom_in| image:: _images/zoom_in.png
   :height: 2ex

.. |zoom_out| image:: _images/zoom_out.png
   :height: 2ex

.. |reset| image:: _images/reset.png
   :height: 2ex

.. |left| image:: _images/left.png
   :height: 2ex

.. |right| image:: _images/right.png
   :height: 2ex

.. |up| image:: _images/up.png
   :height: 2ex

.. |down| image:: _images/down.png
   :height: 2ex

.. |add_node| image:: _images/add_node.png
   :height: 2ex

.. |remove| image:: _images/remove.png
   :height: 2ex

.. |add_edge| image:: _images/add_edge.png
   :height: 2ex

.. |edit| image:: _images/edit.png
   :height: 2ex

.. |enable_physics| image:: _images/enable_physics.png
   :height: 2ex

.. |disable_physics| image:: _images/disable_physics.png
   :height: 2ex

In the interactive visualization, users can:

* Zoom in and out of the network (mouse wheel, keyboard shortcuts, or the zoom buttons |zoom_out| |zoom_in| in the
  bottom right).

* Reset the zoom level (using the reset button |reset| in the bottom right).

* Pan the network (click in any background region and drag, keyboard arrows, or pan buttons |left| |right| |up| |down|
  in the bottom left).

* Drag nodes to reposition them (if physics is enabled, the simulated forces will try to find a new stable
  configuration).

* Hover over nodes to display their names/group labels.

* Click on nodes to highlight their neighbors (multiple selection is allowed). After selecting one or multiple nodes,
  a blue button will appear to allow the user to hide the unselected nodes/neighbours.

* Click on the background to deselect nodes/edges, or to restore the visualization with the full network after node
  hiding.

* Add or remove nodes/edges through the edit toolbar on the top (add node |add_node|, add edge |add_edge|, edit edge
  |edit|, delete selected nodes/edges |remove|).

* Select (multiple) node(s) by their IDs, or any network item (node/edge) that satisfy certain properties.

  .. admonition:: Select a Community of Nodes
       :class: tip

       To select all the nodes belonging to the same community (e.g. p53 signaling pathway), e.g. in the interactive
       `original network, with nodes color-coded by community <original_graph_with_partition_colors.html>`_, use the
       second filer and set respectively the drop-down menus to ``node``, ``title``, ``p53 signaling pathway``, then
       click on the blue `Filter` button.

* Enable or disable the physics via the switch button |disable_physics| / |enable_physics| in the top right.

  .. admonition:: Physics Options
       :class: warning

        Please note that physics is automatically b disabled, once stabilisation is concluded, for large networks to
        prevent the nodes from wiggling. For smaller graphs, the ``forceAtlas2Based`` physics solver is kept active for
        a smoother interaction. If the user enable physics for a big network, the solver is instead ``barnesHut``,
        customised with some constants that are more suitable for a fast stabilisation and responsive interaction.

Sankey Diagram Example
~~~~~~~~~~~~~~~~~~~~~~

An useful visualization for comparing different compression results is the Sankey diagram (see
:func:`~sankey.produce_sankey` function for an example of usage). This kind of diagram shows the flow of genes in
their respective supernodes across compressed networks, for different resolutions of the Louvain algorithm, or for
different ``k`` values in the Hierarchical Clustering algorithm.

This advanced plotting functionality is available through the command line interface using the following command:

.. code:: bash

    gp-sankey --help

This command will display the help message for the Sankey diagram functionality, providing information on how to run
the demo and specify the input graph, output directory, compression method and parameters, and visualization options.

**Example Usage**

.. admonition:: Required Directories' Structure
        :class: note

        The Sankey diagram functionality requires a specific directory structure to work properly. The output directory
        should contain the compressed networks generated by the compression methods, with the following structure:

        .. code:: bash

                data
                ├── output
                    ├── {graph_name}
                        ├── {method}_{param_name}_{param_value_1}
                            ├── compressed_graph.txt
                            ├── compression_mapping.json
                            ├── decompression_mapping.json
                        ├── {method}_{param_name}_{param_value_2}
                            ├── compressed_graph.txt
                            ├── compression_mapping.json
                            ├── decompression_mapping.json
                        ├── {method}_{param_name}_{param_value_3}
                            ├── compressed_graph.txt
                            ├── compression_mapping.json
                            ├── decompression_mapping.json

        Where ``{method}_{param_name}_{param_value_X}`` can be for example, ``louvain_r_1.25`` or ``louvain_r_3.0``.

        Please note that, after running the same compression method on the same input graph, with different parameters,
        all the necessary files are saved under ``data/output/{graph_name}/{method}_{param_name}_{param_value}``, so
        there is no need to move or copy the files around. The Sankey diagram functionality will automatically find the
        files and generate the interactive diagram.

        Of course, the user can specify a different input directory, using the ``--input-path``/``input_path``
        parameter, and change the output directory, using the ``--output-folder``/``output_folder`` parameter.

For example, one might run the Node2Vec algorithm with different values of the number of dimensions (``k``), and then
run the Sankey plotting functionality with the following command:

.. code:: bash

    gp-sankey --graph graph_name --method node2vec --parameter k --values 100 1000 3000

**Example Output**

This is an example of Sankey diagram, showing the flow of nodes between compressed networks, for different values of
the number of dimensions in the Node2Vec algorithm. It is also possible to

.. image:: _images/sankey_diagram.png
    :width: 800px

For the interactive visualization, open the `interactive diagram <sankey_diagram.html>`_.  The interactive Sankey
diagram allows users to explore the flow of nodes between compressed networks, by hovering over the nodes.

.. admonition:: How to deal with Small Clusters?
        :class: hint

        Sometimes, the Sankey diagram may show small clusters that are not very informative. In such cases, it is
        possible to filter out small clusters by setting a minimum size for the clusters. This can be done using the
        ``--min-size``/``min_size`` parameter, that by default is set to 100 original nodes.

.. admonition:: Interactive Sankey Diagram
        :class: warning

        If the original graph contains a large number of nodes, the interactive Sankey diagram may be slow to load, and
        the responsiveness of the hover effect may be affected.
