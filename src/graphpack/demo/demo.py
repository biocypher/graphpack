import os
import sys

import subprocess
import requests
import time

from graphpack.utils import *


def fetch_api_data(request_url, params=None):
    """
    Fetches data from an API.

    Args:
        request_url (str): The URL of the API endpoint.
        params (dict, optional): Parameters to be sent with the request. Defaults to None.

    Returns:
        dict: JSON response from the API.

    Raises:
        Exception: If there's an issue with the request or parsing the response.
    """
    try:
        if params is not None:
            response = requests.post(request_url, data=params)
        else:
            response = requests.post(request_url)

        response.raise_for_status()  # Raise an exception for 4xx and 5xx status codes
        json_data = response.json()
        return json_data

    except Exception as e:
        raise Exception(f"Error fetching API data: {e}")


def parse_graph_data(json_data):
    """
    Parses data from the STRING API response into a NetworkX graph.

    Args:
        json_data (list): List of dictionaries containing data from the API response.

    Returns:
        nx.Graph: A NetworkX graph representing the parsed data.

    Raises:
        ValueError: If the JSON data is invalid or missing required fields.
    """
    try:
        edges = [(data['preferredName_A'], data['preferredName_B'], data['score']) for data in json_data]
        graph = nx.Graph()
        graph.add_weighted_edges_from(edges)
        return graph

    except KeyError as e:
        raise ValueError(f"Invalid JSON data: missing key {e}")

    except Exception as e:
        raise ValueError(f"Error parsing JSON data: {e}")


def get_example_graph(demo_number=0, use_gene_list=False, gene_list=None, custom_params=None, api_url=None):
    """
    Retrieves and parses an example graph from the STRING database based on the specified demo number,
    gene list usage, and additional API parameters.

    Args:
        demo_number (int): The demo number specifying which example query to use.
        use_gene_list (bool): Whether to use a predefined gene list for the query.
        gene_list (list, optional): List of genes to use for the query. Defaults to None.
        custom_params (dict, optional): Custom API parameters to override or supplement defaults. Defaults to None.
        api_url (str, optional): Custom API URL to fetch the data from. Defaults to None.

    Returns:
        nx.Graph: A NetworkX graph representing the retrieved data.
    """
    # Get API parameters and request URL
    params, request_url = get_api_params(demo_number, use_gene_list, gene_list, custom_params, api_url)

    # Fetch STRING data
    json_data = fetch_api_data(request_url, params)

    # Parse data into NetworkX graph and return it
    return parse_graph_data(json_data)


def get_api_params(demo_number=0, use_gene_list=False, gene_list=None, custom_params=None, api_url=None):
    """
    Generates API request parameters and URL based on the specified demo number, gene list usage,
    and custom parameters.

    Args:
        demo_number (int): The demo number specifying which example query to use.
        use_gene_list (bool): Whether to use a predefined gene list for the query.
        gene_list (list, optional): List of genes to use for the query. Defaults to None.
        custom_params (dict, optional): Custom API parameters to override or supplement defaults. Defaults to None.
        api_url (str, optional): Custom API URL to fetch the data from. Defaults to None.

    Returns:
        tuple: A tuple containing the API parameters and request URL.
    """
    if custom_params is None:
        custom_params = {}

    params = None

    if use_gene_list:
        if gene_list is None:
            # Select gene list based on the demo number
            if demo_number == 1:
                gene_list = ['ATP1A1', 'BRCA1', 'CDK2', 'DNMT1', 'EGFR', 'FGFR2', 'GATA1', 'HIF1A', 'IGF1R', 'JAK2']
            elif demo_number == 2:
                gene_list = ['TP53', 'EGFR', 'CDK2', 'BRCA1', 'PTEN', 'PIK3CA', 'RB1', 'ATM', 'MDM2', 'CTNNB1',
                             'AKT1', 'BRAF', 'KRAS', 'SMAD4', 'NOTCH1', 'MYC', 'CCND1', 'NRAS', 'SRC', 'ERBB2',
                             'APC', 'FBXW7', 'GATA3', 'MTOR', 'JAK2', 'PIK3R1', 'FGFR1', 'VHL', 'AR', 'ABL1']
            elif demo_number == 3:
                gene_list = ['ATP1A1', 'BRCA1', 'CDK2', 'DNMT1', 'EGFR', 'FGFR2', 'GATA1', 'HIF1A', 'IGF1R', 'JAK2',
                             'KCNQ1', 'LMO2', 'MAPK1', 'NFKB1', 'ODC1', 'PDGFRA', 'PTEN', 'RARA', 'SMAD4', 'TGFB1',
                             'UBE2C', 'VHL', 'WNT1', 'XBP1', 'YAP1', 'ZEB1', 'ACTB', 'BRAF', 'CBL', 'DDR2',
                             'EZH2', 'FGF23', 'GNAQ', 'HDAC2', 'IDH1', 'JUN', 'KDM5A', 'LATS1', 'MMP9', 'NFATC1',
                             'PAX5', 'PRDM1', 'RAD51', 'SRC', 'TERT', 'TGFBR2', 'USP9X', 'VEGFA', 'YY1', 'ZNF217']
            else:
                gene_list = ['TP53', 'BRCA1', 'EGFR', 'MYC', 'AKT1', 'PIK3CA', 'PTEN', 'RB1', 'KRAS', 'MAPK1',
                             'CDKN1A', 'MDM2', 'CTNNB1', 'APC', 'GAPDH', 'CDK2', 'FOXP3', 'STAT3', 'SMAD4', 'VEGFA',
                             'HIF1A', 'CCND1', 'ERBB2', 'BCL2', 'JAK2']

        # Set API parameters
        species = 9606
        required_score = 700
        string_api_url = "https://string-db.org/api"
        output_format = "json"
        method = "network"
        network_flavour = 'confidence'
        network_type = 'functional'
        request_url = "/".join([string_api_url, output_format, method])

        # Update params with custom_params
        params = {
            "identifiers": "%0d".join(gene_list),
            "species": species,
            'caller_identity': 'my_app',
            'required_score': required_score,
            'network_flavor': network_flavour,
            'network_type': network_type
        }
        params.update(custom_params)

    elif api_url is not None:
        request_url = api_url

    else:
        # Example queries based on demo_number
        if demo_number == 1:
            request_url = "https://string-db.org/api/json/interaction_partners?identifiers=TP53%0dCDK2&limit=20"
        elif demo_number == 2:
            request_url = "https://string-db.org/api/json/interaction_partners?identifiers=P53_HUMAN%0dMDM2_HUMAN%0dATM_HUMAN&species=9606&required_score=900"
        elif demo_number == 3:
            request_url = "https://string-db.org/api/json/interaction_partners?identifiers=TP53%0dBRCA1%0dEGFR%0dMYC%0dAKT1%0dPIK3CA%0dPTEN%0dRB1%0dKRAS%0dMAPK1&species=9606&required_score=950"
        elif demo_number == 4:
            request_url = "https://string-db.org/api/json/interaction_partners?identifiers=AKT1%0dMAPK1%0dCTNNB1%0dTGFB1%0dNFKB1%0dBAX%0dTP53%0dNOTCH1%0dSTAT3%0dKRAS%0dYAP1&species=9606&required_score=990"
        elif demo_number == 5:
            request_url = "https://string-db.org/api/json/interaction_partners?identifiers=AKT1%0dMAPK1%0dCTNNB1%0dTGFB1%0dNFKB1%0dBAX%0dTP53%0dNOTCH1%0dSTAT3%0dKRAS%0dYAP1&species=9606&limit=20"
        else:
            request_url = "https://string-db.org/api/json/network?identifiers=all"

    return params, request_url


def run_demo(demo=None, use_gene_list=False, gene_list=None, custom_params=None, api_url=None,
             input="data/input/simple_graph.txt",
             output="data/output", output_format="txt", resolution=None, k=None, is_weighted=False,
             is_gene_network=False, is_lossless=False, no_plot=False, no_interactive_plot=False,
             plot_disconnected=False, separate_communities=False, method=["greedy"]):
    """
    Main function for the GraphPack tool demo script.

    This script demonstrates the usage of the GraphPack tool by fetching data from the STRING database,
    parsing it into a NetworkX graph, and performing various compression methods on the graph.

    Args:
        demo (int, optional): Demo number to run (0-6). Specifies which example query to use for fetching data from the STRING database.
        use_gene_list (bool, optional): Flag to use a predefined gene list for the query.
        gene_list (list of str, optional): List of genes to use for the query.
        custom_params (dict, optional): Custom parameters to be sent with the API request. Should be in the format '{"key1": "value1", "key2": "value2"}'.
        api_url (str, optional): Custom API URL to fetch the data from.
        input (str, optional): Path to the input file. Default is "data/input/simple_graph.txt".
        output (str, optional): Path to the output folder. Default is "data/output".
        output_format (str): File format to save the network files. Options: '.edgelist', '.txt' (default), '.csv', '.tsv', '.json', '.gpickle', '.gml', '.graphml', '.net', '.pajek', '.gexf', '.yaml', '.yml'.
        resolution (float, optional): Resolution parameter for Louvain and greedy methods.
        k (int, optional): Number of clusters for clustering methods. Only applicable for methods requiring a cluster count (e.g., 'asyn_fluidc', 'spectral', 'hclust', 'node2vec', 'deepwalk', 'nmf').
        is_weighted (bool, optional): Flag to indicate the use of a weighted graph.
        is_gene_network (bool, optional): Flag to indicate the use of a gene network.
        is_lossless (bool, optional): Flag to indicate the use of lossless compression.
        no_plot (bool, optional): Flag to disable plotting the original and compressed graphs. Default is False.
        no_interactive_plot (bool, optional): Flag to disable plotting the graph in interactive mode.
        plot_disconnected (bool, optional): Flag to plot all nodes in a disconnected graph, not just the largest connected component.
        separate_communities (bool, optional): Flag to enforce separation of identified communities in the plots.
        method (list of str, optional): List of compression methods to run. Options: 'louvain', 'greedy' (default), 'label_propagation', 'asyn_fluidc', 'spectral', 'hclust', 'node2vec', 'deepwalk', 'cpm', 'nmf'.

    Returns:
        None

    Examples:
        >>> from graphpack.demo.demo import *

        >>> # Run a default demo, greedy method
        >>> run_demo(demo=5, input='./', output='results')
        >>> # Run a demo with predefined gene list, weighted, non default method and parameters
        >>> run_demo(demo=3, use_gene_list=True, method=['hclust'], k=11, is_gene_network=True, is_weighted=True, input='./', output='results')
        >>> # Run demo pipeline with custom gene list, multiple methods
        >>> run_demo(use_gene_list=True, gene_list=['TP53', 'BRCA1', 'EGFR', 'MYC', 'AKT1', 'PIK3CA', 'PTEN', 'RB1', 'KRAS', 'MAPK1'], is_gene_network=True, method=['louvain', 'greedy', 'hclust'], resolution=1.0, k=4, input='./', output='results')
        >>> # Run demo pipeline with custom gene network, suppress plots
        >>> run_demo(input="path/to/edgelist.txt", is_gene_network=True, is_weighted=True, no_plot=True, output='results')
        >>> # Run demo pipeline with custom API URL and non defgault output format
        >>> run_demo(api_url="https://string-db.org/api/json/interaction_partners?identifiers=TP53%0dBRCA1%0dEGFR%0dMYC&required_score=990", output_format="json", input='./', output='results')
        >>> # Run demo pipeline with predefined gene list, custom parameters
        >>> run_demo(demo=1, use_gene_list=True, custom_params={"required_score": 100}, input='./', output='results')
        >>> # Run demo pipeline with custom gene list, custom parameters
        >>> run_demo(demo=1, use_gene_list=True, custom_params={"required_score": 100}, input='./', output='results')
    """
    # Get the input graph from STRING DB via API
    if demo is not None or use_gene_list or gene_list is not None or api_url is not None:

        graph = None

        if custom_params is not None:
            print(f"{BOLD}Using custom parameters: {custom_params}{RESET}")

        # Get and save the example graph
        if demo is not None:
            if use_gene_list:
                input = os.path.dirname(input) + f"/demo_gene_list_{demo}"
                input = input + "_custom_params.txt" if custom_params is not None else input + ".txt"
            else:
                input = os.path.dirname(input) + f"/demo_{demo}.txt"

            graph = get_example_graph(demo, use_gene_list, custom_params=custom_params)

        if gene_list is not None:
            print(f"{BOLD}Using custom gene list: {gene_list}{RESET}")
            input = os.path.dirname(input) + f"/custom_gene_list"
            input = input + "_custom_params.txt" if custom_params is not None else input + ".txt"

            graph = get_example_graph(use_gene_list=use_gene_list, gene_list=gene_list, custom_params=custom_params)

        if api_url is not None:
            print(f"{BOLD}Using custom API URL: {api_url}{RESET}")
            input = os.path.dirname(input) + f"/custom_api_data.txt"
            graph = get_example_graph(api_url=api_url)

        if graph.number_of_nodes() == 0:
            print(f"{RED_BOLD}Error: No data found for the provided custom parameters/API url/gene list.{RESET}")
            sys.exit(1)

        # Save the graph
        def is_weighted(graph):
            for u, v, data in graph.edges(data=True):
                if 'weight' in data:
                    return True
            return False

        save_graph(graph, input, save_data=any('weight' in data for _, _, data in graph.edges(data=True)))
        print(f"{BOLD}Graph saved in {input}{RESET}")

    # Run the main compression script
    compression_methods = method

    for compression_method in compression_methods:

        call_args = [
            'graphpack', #'python3', '../compression.py',
            '--input', input,
            '--output', output,
            '--output-format', output_format,
            '--method', compression_method,
            '--seed', str(SEED),
            '--title', compression_method
        ]

        if is_weighted:
            call_args.append('--is-weighted')

        if is_gene_network:
            call_args.append('--is-gene-network')

        if is_lossless:
            call_args.append('--is-lossless')

        if not no_plot:
            call_args.append('--plot')

        if not no_interactive_plot:
            call_args.append('--is-interactive')

        if compression_method in ['louvain', 'greedy'] and resolution is not None:
            call_args.extend(['--resolution', str(resolution)])

        if compression_method in ['asyn_fluidc', 'spectral', 'hclust', 'node2vec', 'deepwalk', 'nmf']:
            call_args.extend(['--k', str(k) if k is not None else '30'])

        if compression_method == 'cpm':
            call_args.extend(['--k', str(k) if k is not None else '3'])

        if plot_disconnected:
            call_args.append('--plot-disconnected')

        if separate_communities:
            call_args.append('--separate-communities')

        # Run the compression method and measure runtime
        start_time = time.time()
        subprocess.run(call_args)
        end_time = time.time()

        runtime = end_time - start_time
        print(f"\n🕑 Runtime for {compression_method}: {runtime:.2f} seconds")


def parse_args():
    """
    Parse command-line arguments for graph compression demo script.

    Command-line arguments:

    Args:
        --demo (int): Demo number to run (0-6). Specifies which example query to use for fetching data from the STRING database.\n
        --use-gene-list (bool): Flag to use a predefined gene list for the query.\n
        --gene-list (list of str): List of genes to use for the query. Ignored if --use-gene-list is specified.\n
        --custom-params (list of str): Custom parameters to be sent with the API request. Only applicable when --use-gene-list is set. Should be in the format 'key1=value1 key2=value2'.\n
        --api-url (str): Custom API URL to fetch the data from. If not specified, the default STRING API URL is used.\n
        --input (str): Path to the input file. Default is "data/input/simple_graph.txt".\n
        --output (str): Path to the output folder. Default is "data/output".\n
        --output-format (str): Output file format. Options: '.edgelist', '.txt', '.csv', '.tsv', '.json', '.gpickle', '.gml', '.graphml', '.net', '.pajek', '.gexf', '.yaml', '.yml'. Default is 'txt'.\n
        --resolution (float): Resolution parameter for Louvain and greedy methods.\n
        --k (int): Number of clusters for clustering methods. Only applicable for methods requiring a cluster count (e.g., 'asyn_fluidc', 'spectral', 'hclust', 'node2vec', 'deepwalk', 'nmf').\n
        --is-weighted (bool): Flag to indicate the use of a weighted graph.\n
        --is-gene-network (bool): Flag to indicate the use of a gene network.\n
        --is-lossless (bool): Flag to indicate the use of lossless compression.\n
        --no-plot (bool): Flag to disable plotting the original and compressed graphs. Default is False.\n
        --no-interactive-plot (bool): Flag to disable plotting the graph in interactive mode.\n
        --plot-disconnected (bool): Flag to plot all nodes in a disconnected graph, not just the largest connected component.\n
        --method (list of str): List of compression methods to run. Options: 'louvain', 'greedy' (default), 'label_propagation', 'asyn_fluidc', 'spectral', 'hclust', 'node2vec', 'deepwalk', 'cpm', 'nmf'. If 'all', all the implemented compression methods will be used.\n

    Returns:
        args (argparse.Namespace): Parsed command-line arguments.
    """
    from graphpack import __version__

    # Print the title and visualization of GraphPack tool
    print(f"{BOLD}{TITLE}{RESET}{VIZ}\n")

    # Create the custom parser
    parser = CustomArgumentParser(
        description=f"{DESCR}",
        epilog="For more information, please refer to the documentation.",
        add_help=False
    )

    # Create groups for different sets of options
    help_group = parser.add_argument_group(f'{BOLD}Options{RESET}')
    demo_group = parser.add_argument_group(f'{BOLD}Demo Options{RESET}')
    input_output_group = parser.add_argument_group(f'{BOLD}Input/Output Options{RESET}')
    method_group = parser.add_argument_group(f'{BOLD}Method Options{RESET}')
    graph_options_group = parser.add_argument_group(f'{BOLD}Graph and Compression Options{RESET}')
    plotting_group = parser.add_argument_group(f'{BOLD}Plotting Options{RESET}')

    # Add arguments to the parser
    help_group.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit.')
    help_group.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}',
                            help='Show the version of the program.')

    # Demo options
    demo_group.add_argument('-d', '--demo', type=int,
                            help="Demo number to run. Accepts values from 0 to 6.")
    demo_group.add_argument('-D', '--use-gene-list', action="store_true",
                            help="Use a predefined gene list corresponding to the demo number (0-3). Overrides --gene-list if both are provided.")
    demo_group.add_argument('-L', '--gene-list', type=str, nargs="+",
                            help="List of genes to use for the query. Ignored if --use-gene-list and --demo are specified.")
    demo_group.add_argument('-c', '--custom-params', type=str, nargs="+",
                            help="Custom parameters to be sent with the API request. Should be in the format 'key1=value1 key2=value2'.")
    demo_group.add_argument('-a', '--api-url', type=str,
                            help="Custom API URL to fetch the data from. If not specified, the default STRING API URL is used. It overrides --demo and --use-gene-list.")

    # Input/Output options
    input_output_group.add_argument('-i', '--input', type=str, default="data/input/simple_graph.txt",
                                    help="Path to the input file containing the graph data. Default is 'data/input/simple_graph.txt'. If --demo is specified, this argument is ignored.")
    input_output_group.add_argument('-o', '--output', type=str, default="data/output",
                                    help="Path to the output folder where results will be saved. Default is 'data/output'.")
    input_output_group.add_argument('-f', '--output-format', type=str, default='txt',
                                    help=f'Output file format. Options: {EXTENSIONS} Default is "txt"')

    # Method options
    method_group.add_argument('-m', '--method', type=str, nargs="+", default=["greedy"],
                              help=f"List of compression methods to run. If 'all', all the implemented compression methods will be used. Default is ['greedy']. Accepts multiple methods, among {METHODS}.")
    method_group.add_argument('-r', '--resolution', type=float,
                              help="Resolution parameter for Louvain and greedy methods. Higher values lead to more communities. If not specified, the default value (1.25) is used.")
    method_group.add_argument('-k', '--k', type=int,
                              help="Number of clusters for clustering methods. Relevant for methods requiring a predefined number of clusters. Default is 30. For 'cpm' is the size of the smallest clique. Default is 3.")

    # Graph options
    graph_options_group.add_argument('-w', '--is-weighted', action="store_true",
                                     help="Indicates that the graph is weighted. If not specified, the graph is considered unweighted.")
    graph_options_group.add_argument('-g', '--is-gene-network', action="store_true",
                                     help="Indicates that the input graph is a gene network. If not specified, the graph is considered a general network.")
    graph_options_group.add_argument('-l', '--is-lossless', action="store_true",
                                     help="Use lossless compression for the graph data. If not specified, lossy compression is used.")

    # Plotting options
    plotting_group.add_argument('--no-plot', action="store_true",
                                help="Disable plotting the original and compressed graphs.")
    plotting_group.add_argument('--no-interactive-plot', action="store_true",
                                help="Disable saving the plots also as interactive html files. The static plots will still be saved as images.")
    plotting_group.add_argument('--plot-disconnected', action='store_true',
                                help="Plot all the nodes in a disconnected graph, not only the largest connected component.")
    plotting_group.add_argument('--separate-communities', action='store_true',
                                help="Enforces separation of the identified communities in the plots.")

    # Parse arguments
    args = parser.parse_args()

    if len(sys.argv) - 1 == 0:
        print(LONG_DESCR)
        print(f"{ORANGE_BOLD}Warning: no parameters provided. To display the help, run the script with --help{RESET}")
        print(f"{BOLD}\nRunning the default demo?{RESET}")

        ans = input(f"Enter 'y' to continue or 'n' to exit: ")
        if ans.lower() != 'y':
            sys.exit(0)

        # Running the default demo
        print(f"\n{GREEN_BOLD}Running the default demo (0) with the greedy method.{RESET}\n")
        args.demo = 0

    # Validate arguments
    if args.demo is None and args.use_gene_list is False and args.gene_list is None and args.api_url is None and args.input is None:
        print(f"{RED_BOLD}Error: At least one of --input, --demo, --gene-list, or --api-url must be provided.{RESET}")
        print(f"{BOLD}--demo can be used alone, followed by the demo number, or with --use-gene-list.{RESET}")
        print(f"{BOLD}--gene-list must be used with --use-gene-list.{RESET}")
        print(f"{BOLD}--api-url can be used alone, with a custom STRING DB API URL.{RESET}")
        sys.exit(1)

    # Vallidate api_url argument
    if args.api_url is not None:
        if args.demo is not None or args.use_gene_list or args.gene_list is not None:
            print(
                f"{ORANGE_BOLD}Warning: --api-url is only applicable when --demo is not set and --use-gene-list and --gene-list are not provided.{RESET}")
            print(f"{BOLD}Ignoring the custom API URL anyway.{RESET}")
            args.api_url = None

    # Validate demo argument
    if args.demo is not None:
        if args.use_gene_list:
            if not (0 <= args.demo <= 3):
                print(
                    f"{ORANGE_BOLD}Warning: when --use-gene-list is set, --demo must be an integer between 0 and 3.{RESET}")
                print(f"{BOLD}Running the demo with default gene list (0).{RESET}")
                args.demo = 0

        elif args.custom_params is not None:
            print(f"{ORANGE_BOLD}Warning: --custom-params is only applicable when --use-gene-list is set.{RESET}")
            print(f"{BOLD}Ignoring the custom parameters.{RESET}")
            args.custom_params = None

        else:
            if not (0 <= args.demo <= 6):
                print(f"{ORANGE_BOLD}Warning: --demo must be an integer between 0 and 6.{RESET}")
                print(f"{BOLD}Running the default demo (0).{RESET}")
                args.demo = 0

    # Validate gene list argument
    if args.gene_list is not None:
        if not args.use_gene_list:
            print(f"{ORANGE_BOLD}Warning: --gene-list is only applicable when --use-gene-list is set.{RESET}")
            print(f"{BOLD}Attempting to use the provided gene list anyway.{RESET}")
            args.use_gene_list = True

    # Validate custom parameters argument
    if args.custom_params is not None:
        try:
            args.custom_params = {param.split('=')[0]: param.split('=')[1] for param in args.custom_params}
        except ValueError:
            print(f"{RED_BOLD}Error: --custom-params must be in the format 'key1=value1 key2=value2'.{RESET}")
            sys.exit(1)

    # Validate input file path
    if not os.path.isfile(args.input) and args.demo is None:
        print(f"{RED_BOLD}Error: Input file '{args.input}' does not exist.")
        sys.exit(1)

    # Validate input folder path
    if not os.path.isdir(os.path.dirname(args.input)):
        print(f"{ORANGE_BOLD}Warning: Input folder '{os.path.dirname(args.input)}' does not exist.{RESET}")
        print(f"{BOLD}Creating it now.{RESET}")
        os.makedirs(os.path.dirname(args.input))

    # Validate output folder path
    output_path = os.path.abspath(args.output)

    if not os.path.isdir(output_path):
        print(f"{ORANGE_BOLD}Warning: Output folder '{output_path}' does not exist.")
        print(f"{BOLD}Creating it now.{RESET}")

        try:
            os.makedirs(output_path, exist_ok=True)
        except Exception as e:
            print(f"{RED_BOLD}Error: Failed to create output folder '{output_path}'. {str(e)}{RESET}")
            sys.exit(1)

    # Validate output format
    if args.output_format not in EXTENSIONS:
        print(f"{RED_BOLD}Error: Unsupported output format '{args.output_format}'.{RESET}")
        print(f"{BOLD}Supported formats are: {', '.join(EXTENSIONS)}{RESET}")
        print(f"{ORANGE_BOLD}Warning: Running the compression with the default output format ('txt').{RESET}")

    # Validate resolution parameter
    if args.resolution is not None:
        if args.resolution <= 0:
            print(f"{RED_BOLD}Error: --resolution must be a positive float.{RESET}")
            sys.exit(1)

    # Validate number of clusters (k)
    if args.k is not None:
        if args.k <= 0:
            print(f"{RED_BOLD}Error: --k must be a positive integer.{RESET}")
            sys.exit(1)

    # Validate compression methods
    for method in args.method:
        if method == 'all':
            args.method = METHODS
            break
        if method not in METHODS:
            print(f"{RED_BOLD}Error: Unsupported compression method '{method}'.{RESET}")
            print(f"{BOLD}Supported methods are: {', '.join(METHODS)}{RESET}")
            sys.exit(1)

    return args


def main():
    args = parse_args()
    run_demo(**vars(args))

if __name__ == "__main__":
    main()
