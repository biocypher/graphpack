import sys

import plotly.graph_objects as go

from graphpack.utils import *


MARGIN = 50  # Margin for the plot
MAX_MARGIN = 500  # Maximum margin for the plot (increase if last layer nodes' labels are cut off)


def lighten_color(color, amount=0.5):
    """
    Lighten the given color by the specified amount.

    Args:
        color (str): A color name or RGBA string.
        amount (float, optional): The amount to lighten the color (0.0 to 1.0).

    Returns:
        str: The lightened color as an RGBA string.
    """
    # Check if the color is already in RGBA format
    rgba_pattern = re.compile(r'rgba\((\d+),\s*(\d+),\s*(\d+),\s*(\d*\.?\d+)\)')
    match = rgba_pattern.match(color)

    if match:
        # Extract RGBA components from the string
        r, g, b, a = map(float, match.groups())
        r, g, b = int(r), int(g), int(b)
        # Apply the amount to the alpha component
        return f'rgba({r},{g},{b},{amount})'
    else:
        # Convert named color to RGBA
        c = mcolors.to_rgba(color)
        return f'rgba({c[0] * 255:.0f},{c[1] * 255:.0f},{c[2] * 255:.0f},{amount})'


def load_data(input_path, graph, method, parameter, parameters):
    """
    Load compression mappings and group labels for each parameter.

    Args:
        input_path (str): Path to the input files.
        graph (str): Input graph identifier.
        method (str): Clustering method used.
        parameter (str): Clustering parameter name.
        parameters (list): List of parameters to be analyzed.

    Returns:
        tuple: Two dictionaries, one for compression mappings and one for group labels, and the (eventually updated) list of parameter's values.
    """
    # Initialize dictionaries to store compression mappings and groups for each parameter
    compression_mappings = {}
    groups = {}

    # Create a copy of parameters to iterate over
    parameters_copy = parameters[:]

    # Load compression mappings and group labels for each parameter
    for param in parameters_copy:
        try:
            with open(f'{input_path}/{graph}/{method}_{parameter}_{param}/compression_mapping.json', 'r') as f:
                compression_mappings[param] = json.load(f)

            labels_file = f'{input_path}/{graph}/{method}_{parameter}_{param}/labels_mapping.json'
            if os.path.exists(labels_file):
                with open(labels_file, 'r') as f:
                    groups[param] = json.load(f)
            else:
                # If labels_mapping.json doesn't exist, use names from compression_mapping.json as group labels
                compression_data = compression_mappings[param]
                groups[param] = {key: key for key in compression_data}

        except FileNotFoundError:
            print(f"Files for parameter {param} not found. Skipping this value.")
            parameters.remove(param)

    return compression_mappings, groups, parameters


def map_transitions(compression_mappings, groups, parameters, min_size):
    """
    Map transitions between consecutive parameters.

    Args:
        compression_mappings (dict): Compression mappings for each parameter.
        groups (dict): Group labels for each parameter.
        parameters (list): List of parameters to be analyzed.
        min_size (int): Minimum cluster size to be considered significant.

    Returns:
        list: List of transitions between clusters.
    """
    # Initialize a list to store the gene transitions between resolutions
    transitions = []
    small_cluster_label = "small clusters"

    # Iterate through each pair of consecutive resolutions to map transitions
    for i in range(len(parameters) - 1):
        param_from = parameters[i]
        param_to = parameters[i + 1]

        if param_from in compression_mappings and param_to in compression_mappings:
            mapping_from = compression_mappings[param_from]
            mapping_to = compression_mappings[param_to]
            group_from = groups[param_from]
            group_to = groups[param_to]

            # Create a reverse mapping from gene to cluster for the target resolution
            reverse_mapping_to = {gene: cluster for cluster, genes in mapping_to.items() for gene in genes}

            # Map transitions from the current resolution to the next
            for cluster_from, genes_from in mapping_from.items():
                source_label = f'{param_from} - {group_from.get(str(cluster_from), small_cluster_label)}'
                if len(genes_from) > min_size:  # Main cluster handling
                    for gene in genes_from:
                        if gene in reverse_mapping_to:
                            cluster_to = reverse_mapping_to[gene]
                            target_label = f'{param_to} - {group_to.get(str(cluster_to), small_cluster_label)}'
                            if len(mapping_to[cluster_to]) > min_size:
                                transitions.append([source_label, target_label, gene])
                            else:
                                transitions.append([source_label, f'{param_to} - {small_cluster_label}', gene])
                else:
                    # Handle small clusters at the source resolution
                    for gene in genes_from:
                        if gene in reverse_mapping_to:
                            cluster_to = reverse_mapping_to[gene]
                            target_label = f'{param_to} - {group_to.get(str(cluster_to), small_cluster_label)}'
                            if len(mapping_to[cluster_to]) > min_size:
                                transitions.append([f'{param_from} - {small_cluster_label}', target_label, gene])
                            else:
                                transitions.append(
                                    [f'{param_from} - {small_cluster_label}', f'{param_to} - {small_cluster_label}',
                                     gene])
        else:
            print(f"Skipping transition from resolution {param_from} to {param_to}.")

    return transitions


def create_sankey_plot(transitions, min_size, method, input_graph, output_folder):
    """
    Create and save the Sankey plot.

    Args:
        transitions (list): List of transitions between clusters.
        min_size (int): Minimum cluster size to be considered significant.
        parameter (str): Clustering parameter name.
        method (str): Clustering method used.
        input_graph (str): Knowledge graph identifier.
        output_folder (str): Path to the output folder.

    Returns:
        None
    """
    # Create a DataFrame from the transitions list
    transitions_df = pd.DataFrame(transitions, columns=['Source', 'Target', 'Gene'])

    # Extract unique stages from Source and Target columns
    stages = sorted(list(set(transitions_df['Source'].str.split(' - ').str[0]).union(
        set(transitions_df['Target'].str.split(' - ').str[0]))))

    # Initialize a dictionary to store nodes grouped by stage
    stage_nodes = {stage: [] for stage in stages}

    # Populate stage_nodes dictionary with nodes, sorted alphabetically within each stage
    for stage in stage_nodes:
        stage_nodes[stage].extend(transitions_df[transitions_df['Source'].str.contains(stage)]['Source'].unique())
        stage_nodes[stage].extend(transitions_df[transitions_df['Target'].str.contains(stage)]['Target'].unique())
        stage_nodes[stage] = sorted(list(set(stage_nodes[stage])))

    # Flatten stage_nodes dictionary into a sorted list of nodes
    nodes = [node for stage in stages for node in stage_nodes[stage]]

    # Create a dictionary to map node labels to their indices
    node_indices = {node: idx for idx, node in enumerate(nodes)}

    # Initialize lists to store links with updated indices
    links = []

    # Iterate through each row in transitions_df to create links with updated indices
    for index, row in transitions_df.iterrows():
        source_index = node_indices[row['Source']]
        target_index = node_indices[row['Target']]
        links.append({
            'source': source_index,
            'target': target_index,
            'value': 1
        })

    # Prepare custom data for each node with limited gene list
    node_customdata = []
    gene_display_limit = 10  # Limit the number of genes displayed in hover info

    for node in nodes:
        node_stage = node.split(' - ')[0]
        group = node.split(' - ')[1]
        genes = transitions_df[(transitions_df['Source'] == node) | (transitions_df['Target'] == node)]['Gene'].unique()
        gene_list = ', '.join(genes[:gene_display_limit])  # Limit the number of genes displayed
        if len(genes) > gene_display_limit:
            gene_list += f', ... (+{len(genes) - gene_display_limit} more)'
        custom_data = f"Parameter: {node_stage}<br>Community: {group}<br>Genes: {gene_list}"
        node_customdata.append(custom_data)

    # Extract unique community names
    communities = sorted(list(set(node.split(' - ')[1] for node in nodes)))

    # Generate a larger color palette using matplotlib
    colors = [mcolors.CSS4_COLORS[name] for name in COLORS]
    color_map = {community: colors[i % len(colors)] for i, community in enumerate(communities)}

    # Assign colors to nodes based on their community
    node_colors = []
    for node in nodes:
        community = node.split(' - ')[1]
        if node.endswith('small clusters'):
            node_colors.append('rgba(0, 0, 0, 0)')  # Transparent color for small clusters
        else:
            node_colors.append(color_map[community])

    # Create a list to store line colors based on whether the target node is a small cluster
    line_colors = []
    for link in links:
        target_node = nodes[link['target']]
        source_color = node_colors[link['source']]
        target_color = node_colors[link['target']]
        if target_node.endswith('small clusters'):
            line_colors.append(mcolors.CSS4_COLORS['whitesmoke'])  # Lighter gray for links going into small clusters
        else:
            # Lighter version of the source node color, if big communities, else lightgray for other links
            line_colors.append(lighten_color(source_color, 0.2)) if min_size >= 100 else line_colors.append(
                mcolors.CSS4_COLORS['lightgray'])

    # Create Sankey plot with custom hover information, node colors, and link colors
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,  # Padding between nodes
            thickness=20,  # Thickness of the links
            line=dict(color=mcolors.CSS4_COLORS['lightgray'], width=0.5),
            label=nodes,  # Node labels
            color=node_colors,  # Assign node colors
            customdata=node_customdata,  # Use the prepared custom data
            hovertemplate='%{customdata}<extra></extra>',  # Use the custom data for hover text
            hoverinfo='all',  # Enable hover information for nodes
        ),
        link=dict(
            source=[link['source'] for link in links],  # Indices of source nodes
            target=[link['target'] for link in links],  # Indices of target nodes
            value=[link['value'] for link in links],  # Link values
            hoverinfo='none',  # Disable hover interaction for links
            line=dict(width=0.0005),
            color=line_colors,  # Specify link colors based on conditions
        )
    )])

    # Update the layout of the plot
    fig.update_layout(
        title_text=f"Genes' community membership for {input_graph} - {method}",
        margin=dict(t=MARGIN, l=MARGIN, r=MAX_MARGIN, b=MARGIN),  # Adjusted margins
        font_size=10,  # Adjusted font size for better visibility
        width=2000,  # Increased width for better visibility
        height=2000,  # Increased height for better visibility
        plot_bgcolor='white',
        paper_bgcolor='white',
        font_color='black'
    )

    # JavaScript for aligning node labels to the right
    js = '''
    const TEXTPAD = 3; // constant used by Plotly.js
    
    function sankeyNodeLabelsAlign(position, forcePos) {
      const textAnchor = {left: 'end', right: 'start', center: 'middle'}[position];
      const nodes = gd.getElementsByClassName('sankey-node');
    
      for (const node of nodes) {
        const d = node.__data__;
        const label = node.getElementsByClassName('node-label').item(0);
    
        // Ensure to reset any previous modifications
        label.setAttribute('x', 0);
    
        if (!d.horizontal)
          continue;
    
        // This is how Plotly's default text positioning is computed (coordinates
        // are relative to that of the corresponding node).
        const padX = d.nodeLineWidth / 2 + TEXTPAD;
        const posX = padX + d.visibleWidth;
        let x;
    
        switch (position) {
          case 'left':
            if (d.left || (d.node.originalLayer === 0 && !forcePos))
              continue;
            x = -posX - padX;
            break;
    
          case 'right':
            if (!d.left || !forcePos)
              continue;
            x = posX + padX;
            break;
    
          case 'center':
            if (!forcePos && (d.left || d.node.originalLayer === 0))
              continue;
            x = (d.nodeLineWidth + d.visibleWidth) / 2 + (d.left ? padX : -posX);
            break;
        }
    
        // Ensure last layer nodes' labels are inside the plot area
        if (d.node.originalLayer === d.layerLength - 1) {
          x = Math.min(x, gd.layout.width - label.getBBox().width - padX);
        }
    
        label.setAttribute('x', x);
        label.setAttribute('text-anchor', textAnchor);
      }
    }
    const gd = document.getElementById('{plot_id}');
    const position = 'right'; // Set position to 'right', 'left', or 'center'
    const forcePos = true;
    
    gd.on('plotly_afterplot', sankeyNodeLabelsAlign.bind(gd, position, forcePos));
    gd.emit('plotly_afterplot');
    '''

    # Create output folder if it does not exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    fig.write_html(f"{output_folder}/{input_graph}_sankey_{method}_s_{min_size}.html", post_script=js)


def produce_sankey(graph, input_path='data/output', output_folder='sankey', min_size=100,
                   method='louvain', parameter='r', values=[1.25, 3.0, 5.0, 10.0, 20.0, 30.0]):
    """
    Main function for the GraphPack tool Sankey plot script.

    This script generates a Sankey plot to visualize gene community membership transitions across
    different clustering resolutions for a given network.

    Args:
        graph (str): Input graph identifier. Required parameter.
        input_path (str, optional): Path to the input files. Default is "data/output".
        output_folder (str, optional): Path to the output folder. Default is "sankey".
        min_size (int, optional): Minimum cluster size to be considered significant. Default is 100.
        method (str, optional): Clustering method used. Default is "louvain".
        parameter (str, optional): Clustering parameter name, as it appears in the subfolders' names. Default is "r".
        values (list of float, optional): List of parameters to be analyzed. Default is [1.25, 3.0, 5.0, 10.0, 20.0, 30.0].

    Examples:
        >>> from graphpack.demo.sankey import *
        >>> produce_sankey(input_path="./", output_folder="results", min_size=50, graph='simple_graph', method='hclust', parameter='k', values=[10, 50, 100, 250])
    """
    # Print information about the arguments
    print(f"\n{'=' * 80}")
    print(f"{'GraphPack Tool Sankey plot script':^80}")
    print(f"{'=' * 80}")

    print(f"\n▶ Input graph:       {graph}")
    print(f"▶ Method:           {method}")
    print(f"▶ Parameters:       {parameter} in [ {', '.join(map(str, values))} ]")
    print(f"\n▶ Min cluster size: {min_size}")
    print(f"\n▶ Input folder:     {input_path}")
    print(f"▶ Output folder:    {output_folder}")

    print("\n📑 Loading data...")
    compression_mappings, groups, values = load_data(input_path,
                                                     graph,
                                                     method,
                                                     parameter,
                                                     values)

    print("⚙️ Mapping transitions...")
    transitions = map_transitions(compression_mappings, groups, values, min_size)

    print("📊 Creating Sankey plot...")
    create_sankey_plot(transitions, min_size, method, graph, output_folder)

    print("\n✅  Done!")
    print(f"Sankey plot for gene community membership transitions has been saved in '{output_folder}'.")


def parse_args():
    """
    Parse command-line arguments for Sankey plot script.

    Command-line arguments:

    Args:
        --graph (str): Input graph identifier. Required argument.
        --input-path (str): Path to the input files. Default is "data/output".
        --output-folder (str): Path to the output folder. Default is "sankey".
        --min-size (int): Minimum cluster size to be considered significant. Default is 100.
        --method (str): Clustering method used. Default is "louvain".
        --parameter (str): Clustering parameter name. Default is "r".
        --parameters (list of float): List of parameters to be analyzed. Default is [1.25, 3.0, 5.0, 10.0, 20.0, 30.0].

    Returns:
        args (argparse.Namespace): Parsed command-line arguments.
    """
    from graphpack import __version__


    # Create the custom parser
    parser = CustomArgumentParser(
        description="Generate a Sankey plot for gene community membership transitions.",
        epilog="For more information, please refer to the documentation.",
        add_help=False
    )

    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                        help='Show this help message and exit.')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}',
                        help='Show the version of the program.')

    parser.add_argument('-g', '--graph', type=str, required=True,
                        help='Input graph identifier.')

    parser.add_argument('-i', '--input-path', type=str, default='data/output',
                        help='Path to ihe input files. Default is "data/output".')
    parser.add_argument('-o', '--output-folder', type=str, default='sankey',
                        help='Path to the output folder. Default is "sankey".')

    parser.add_argument('-s', '--min-size', type=int, default=100,
                        help='Minimum cluster size to be considered significant. Default is 100.')

    parser.add_argument('-m', '--method', type=str, default='louvain',
                        help='Clustering method used. Default is "louvain".')
    parser.add_argument('-p', '--parameter', type=str, default='r',
                        help='Clustering parameter name, as it appears in the output subfolder. Default is "r".')
    parser.add_argument('-V', '--values', type=float, nargs='+', default=[1.25, 3.0, 5.0, 10.0, 20.0, 30.0],
                        help='List of parameters to be analyzed. Default is [1.25, 3.0, 5.0, 10.0, 20.0, 30.0].')

    args = parser.parse_args()

    if len(sys.argv) - 1 == 0:
        print(LONG_DESCR)
        print(f"{ORANGE_BOLD}Warning: no parameters provided. To display the help, run the script with --help{RESET}")
        sys.exit(0)

    # Validate input  path
    if not os.path.exists(args.input_path):
        print(f"{RED_BOLD}Error: Input file '{args.input_path}' does not exist.")
        sys.exit(1)


    # Validate output folder path
    output_path = os.path.abspath(args.output_path)

    if not os.path.isdir(output_path):
        print(f"{ORANGE_BOLD}Warning: Output folder '{output_path}' does not exist.")
        print(f"{BOLD}Creating it now.{RESET}")

        try:
            os.makedirs(output_path, exist_ok=True)
        except Exception as e:
            print(f"{RED_BOLD}Error: Failed to create output folder '{output_path}'. {str(e)}{RESET}")
            sys.exit(1)

    # Validate the method argument
    if args.method not in METHODS:
        print(f"{RED_BOLD}Error: Unsupported compression method '{args.method}'.{RESET}")
        print(f"{BOLD}Supported methods are: {', '.join(METHODS)}{RESET}")
        sys.exit(1)

    # Cast to integer the parameter values if the method is not Louvain (k must be integer)
    if args.method != 'louvain':
        args.values = [int(param) for param in args.values]

    return args


def main():
    args = parse_args()
    produce_sankey(**vars(args))

if __name__ == "__main__":
    main()
