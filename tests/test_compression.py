import io
import shutil
import sys
import tempfile
import unittest
from contextlib import contextmanager
from unittest.mock import patch

from graphpack.compression import *


@contextmanager
def suppress_output():
    """
    A context manager that redirects stdout to a mock_stdout.
    """
    # Helper function to suppress print statements during testing
    with patch('sys.stdout', new=io.StringIO()) as mock_stdout:
        yield mock_stdout


class TestCompressionFunctions(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory for testing
        self.temp_dir = tempfile.mkdtemp()

        # Create a sample graph for testing
        self.graph = nx.Graph()
        self.graph.add_weighted_edges_from([(1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 6, 1)])
        self.graph.add_edge(1, 3, weight=2)
        self.graph.add_edge(2, 4, weight=3)

        # Mapping dictionary from numbers to letters
        mapping = {1: 'A', 2: 'B', 3: 'C', 4: 'D', 5: 'E', 6: 'F'}

        # Relabel the nodes using the mapping
        self.graph = nx.relabel_nodes(self.graph, mapping)

    def tearDown(self):
        # Remove the temporary directory
        shutil.rmtree(self.temp_dir)

    def test_detect_communities(self):

        with suppress_output():
            # Test community detection using all methods
            for method in ['louvain', 'greedy', 'label_propagation', 'asyn_fluidc', 'spectral', 'hclust', 'node2vec',
                           'deepwalk', 'cpm', 'nmf']:
                partition = detect_communities(self.graph, method=method)
                self.assertIsInstance(partition, dict)
                self.assertEqual(len(partition), len(self.graph))

                # Convert partition dictionary to list of sets
                communities = {}
                for node, community in partition.items():
                    if community not in communities:
                        communities[community] = set()
                    communities[community].add(node)
                partition_list = list(communities.values())

                self.assertTrue(nx.algorithms.community.is_partition(self.graph, partition_list))

    def test_unsupported_detection_method(self):
        with self.assertRaises(ValueError) as context:
            detect_communities(self.graph, 'unsupported_method')

        self.assertEqual(str(context.exception), "Unsupported detection method: unsupported_method")

    def test_disconnected_graph_detection(self):
        # Test community detection on a disconnected graph
        graph = nx.Graph()
        graph.add_edges_from([(1, 2), (3, 4)])

        with self.assertRaises(Exception) as context:
            _ = detect_communities(graph, method='asyn_fluidc')

        self.assertIn("Graph must be connected for community detection.", str(context.exception))

    def test_compress_network(self):
        # Test compressing the network
        communities = louvain_communities(self.graph)
        partition = {node: cid for cid, community in enumerate(communities) for node in community}
        community_graph, compression_mapping, decompression_mapping = compress_graph_partition_based(self.graph,
                                                                                                     partition,
                                                                                                     is_gene_network=False)
        self.assertIsInstance(community_graph, nx.Graph)
        self.assertIsInstance(compression_mapping, dict)
        self.assertIsInstance(decompression_mapping, dict)
        self.assertEqual(len(compression_mapping), len(set(partition.values())))

    def test_save_network_files(self):
        # Test saving network files and mappings
        community_graph, compression_mapping, decompression_mapping = compress_graph(self.graph, method='louvain',
                                                                                     is_gene_network=False)
        output_folder = self.temp_dir
        save_network_files(self.graph, community_graph, compression_mapping, decompression_mapping, output_folder)
        original_network_file = os.path.join(output_folder, 'original_network.txt')
        compressed_network_file = os.path.join(output_folder, 'compressed_network.txt')
        compression_mapping_file = os.path.join(output_folder, 'compression_mapping.msgpack')
        decompression_mapping_file = os.path.join(output_folder, 'decompression_mapping.msgpack')
        self.assertTrue(os.path.exists(original_network_file))
        self.assertTrue(os.path.exists(compressed_network_file))
        self.assertTrue(os.path.exists(compression_mapping_file))
        self.assertTrue(os.path.exists(decompression_mapping_file))

    def test_calculate_compression_efficacy(self):
        # Test calculating compression efficacy
        community_graph, compression_mapping, _ = compress_graph(self.graph, method='louvain', is_gene_network=False)
        output_folder = self.temp_dir
        save_network_files(self.graph, community_graph, compression_mapping, {}, output_folder)
        compute_and_save_edges_to_remove(self.graph, community_graph, output_folder)

        captured_output = io.StringIO()
        sys.stdout = captured_output
        calculate_compression_efficacy(output_folder)
        sys.stdout = sys.__stdout__
        for item in ['Original size:', 'Compressed size:', 'Percentage compression:',
                     "Lossy compression size:", "Percentage lossy compression:",
                     "Lossless compression size:", "Percentage lossless compression:"]:
            self.assertIn(item, captured_output.getvalue())

    def test_error_calculate_compression_efficacy(self):
        # Test calculating compression efficacy
        community_graph, _, _ = compress_graph(self.graph, method='louvain', is_gene_network=False)
        output_folder = self.temp_dir
        save_network_files(self.graph, community_graph, {}, {}, output_folder)

        captured_output = io.StringIO()
        sys.stdout = captured_output

        with self.assertRaises(Exception) as context:
            calculate_compression_efficacy(output_folder)

        self.assertIn("Error calculating compression efficacy", str(context.exception))

        sys.stdout = sys.__stdout__

    def test_compute_and_save_edges_to_remove(self):
        # Create a sample graph for testing
        graph = nx.Graph()
        graph.add_edges_from([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)])
        graph.add_edge(1, 3, weight=2)
        graph.add_edge(2, 4, weight=3)

        # Create a sample community graph for testing
        community_graph = nx.Graph()
        community_graph.add_nodes_from([(1, {'nodes': [1, 2]}), (2, {'nodes': [3, 4, 5, 6]})])

        # Create a sample partition for testing
        partition = {1: 1, 2: 1, 3: 2, 4: 2, 5: 2, 6: 2}

        # Create a temporary directory for testing
        output_folder = tempfile.mkdtemp()

        try:
            # Call the function being tested
            compute_and_save_edges_to_remove(graph, community_graph, output_folder)

            # Verify that the edges to remove file was successfully created
            edges_to_remove_file = os.path.join(output_folder, 'edges_to_remove.txt')
            self.assertTrue(os.path.exists(edges_to_remove_file))

            # Verify the content of the edges to remove file
            with open(edges_to_remove_file, 'r') as f:
                lines = f.readlines()
                # Check if the number of lines matches the expected number based on the test scenario
                self.assertEqual(len(lines), 3)

        finally:
            # Clean up: Remove the temporary directory
            shutil.rmtree(output_folder)

    def test_reconstruct_original_network(self):
        # Create a temporary directory for testing
        output_folder = tempfile.mkdtemp()

        try:
            # Create sample files for testing
            compressed_network_file = os.path.join(output_folder, 'compressed_network.txt')
            compression_mapping_file = os.path.join(output_folder, 'compression_mapping.msgpack')
            edges_to_remove_file = os.path.join(output_folder, 'edges_to_remove.txt')
            output_file = os.path.join(output_folder, 'reconstructed_network.txt')

            # Create a sample compressed network file
            with open(compressed_network_file, 'w') as f:
                f.write("0\t1\n0\t2\n0\t3\n1\t2\n1\t3\n2\t3\n")

            # Create a sample compression mapping file
            compression_mapping = {0: ["A", "B"], 1: ["C", "D", "E"], 2: ["F", "G", "H", "I"], 3: ["J"]}
            with open(compression_mapping_file, 'wb') as f:
                f.write(msgpack.packb(compression_mapping, use_bin_type=True))

            # Create a sample edges to remove file
            with open(edges_to_remove_file, 'w') as f:
                f.write("A\tI\nB\tI\nC\tF\nC\tI\nD\tF\nD\tI\nE\tF\nE\tI\nI\tJ\n")

            # Call the function being tested
            _ = reconstruct_original_network(output_folder, compressed_network_file,
                                             compression_mapping_file, edges_to_remove_file,
                                             output_file)

            # Verify that the output file was successfully created
            self.assertTrue(os.path.exists(output_file))

            # Verify the content of the output file
            with open(output_file, 'r') as f:
                lines = f.readlines()
                # Check if the number of lines matches the expected number based on the test scenario
                self.assertEqual(len(lines), 36)

                lines = set([line.strip() for line in lines])
                self.assertEqual(lines, {'A\tB', 'A\tC', 'A\tD', 'A\tE', 'A\tF', 'A\tG', 'A\tH', 'A\tJ', 'B\tC',
                                         'B\tD', 'B\tE', 'B\tF', 'B\tG', 'B\tH', 'B\tJ', 'C\tD', 'C\tE', 'C\tG',
                                         'C\tH', 'C\tJ', 'D\tE', 'D\tG', 'D\tH', 'D\tJ', 'E\tG', 'E\tH', 'E\tJ',
                                         'F\tG', 'F\tH', 'F\tI', 'F\tJ', 'G\tH', 'G\tI', 'G\tJ', 'H\tI', 'H\tJ'})

        finally:
            # Clean up: Remove the temporary directory
            shutil.rmtree(output_folder)

    def test_reconstruct_original_network_missing_file(self):
        # Create a temporary directory for testing
        output_folder = tempfile.mkdtemp()

        try:
            # Create sample files for testing
            compressed_network_file = os.path.join(output_folder, 'compressed_network.txt')
            edges_to_remove_file = os.path.join(output_folder, 'edges_to_remove.txt')
            output_file = os.path.join(output_folder, 'reconstructed_network.txt')

            # Create a sample compressed network file
            with open(compressed_network_file, 'w') as f:
                f.write("0\t1\n0\t2\n0\t3\n1\t2\n1\t3\n2\t3\n")

            # Create a sample edges to remove file
            with open(edges_to_remove_file, 'w') as f:
                f.write("A\tI\nB\tI\nC\tF\nC\tI\nD\tF\nD\tI\nE\tF\nE\tI\n")

            # Attempt to call the function with a missing compression mapping file
            with self.assertRaises(ValueError) as context:
                reconstruct_original_network(output_folder, compressed_network_file,
                                             'non_existent_file.msgpack', edges_to_remove_file,
                                             output_file)

            self.assertIn("File 'non_existent_file.msgpack' not found", str(context.exception))

        finally:
            # Clean up: Remove the temporary directory
            shutil.rmtree(output_folder)

    def test_reconstruct_original_network_warnings(self):
        # Create a temporary directory for testing
        output_folder = tempfile.mkdtemp()

        try:
            # Create sample files for testing
            compressed_network_file = os.path.join(output_folder, 'compressed_network.txt')
            compression_mapping_file = os.path.join(output_folder, 'compression_mapping.msgpack')
            edges_to_remove_file = os.path.join(output_folder, 'edges_to_remove.txt')
            output_file = os.path.join(output_folder, 'reconstructed_network.txt')

            # Create a sample compressed network file
            with open(compressed_network_file, 'w') as f:
                f.write("0\t1\n0\t2\n0\t3\n1\t2\n1\t3\n2\t3\n")

            # Create a sample compression mapping file with a non-existent community
            compression_mapping = {0: ["A", "B"], 1: ["C", "D", "E"], 2: ["F", "G", "H", "I"], 3: ["J"], 4: ["K"]}
            with open(compression_mapping_file, 'wb') as f:
                f.write(msgpack.packb(compression_mapping, use_bin_type=True))

            # Create a sample edges to remove file
            with open(edges_to_remove_file, 'w') as f:
                f.write("A\tI\nB\tI\nC\tF\nC\tI\nD\tF\nD\tI\nE\tF\nE\tI\nI\tJ\n")

            # Capture the printed warning
            with suppress_output() as mock_stdout:
                reconstruct_original_network(output_folder, compressed_network_file,
                                             compression_mapping_file, edges_to_remove_file,
                                             output_file)
                # Get the printed output
                printed_output = mock_stdout.getvalue().strip()

            # Check for warning message
            expected_warning = f"{ORANGE_BOLD}Warning: Community 4 not found in the compressed network.{RESET}"
            self.assertIn(expected_warning, printed_output)

            # Create a sample compression mapping file with a non-existent community
            compression_mapping = {0: ["A", "B"], 1: ["C", "D", "E"], 2: ["F", "G", "H", "I"], 3: ["J"]}
            with open(compression_mapping_file, 'wb') as f:
                f.write(msgpack.packb(compression_mapping, use_bin_type=True))

            # Create a sample edges to remove file with a non-existent edge
            with open(edges_to_remove_file, 'w') as f:
                f.write("A\tI\nB\tI\nC\tF\nC\tI\nD\tF\nD\tI\nE\tF\nE\tI\nI\tJ\nC\tX\n")

            # Capture the printed warning
            with suppress_output() as mock_stdout:
                reconstruct_original_network(output_folder, compressed_network_file,
                                             compression_mapping_file, edges_to_remove_file,
                                             output_file)
                # Get the printed output
                printed_output = mock_stdout.getvalue().strip()

            # Check for warning message
            expected_warning = f"{ORANGE_BOLD}Warning: Edge C-X does not exist in the graph.{RESET}"
            self.assertIn(expected_warning, printed_output)
        finally:
            # Clean up: Remove the temporary directory
            shutil.rmtree(output_folder)


class TestMainFunction(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory for testing
        self.temp_dir = tempfile.mkdtemp()

        # Create a sample graph for testing
        self.graph = nx.Graph()
        self.graph.add_weighted_edges_from([(1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 6, 1)])

    def tearDown(self):
        # Remove the temporary directory
        shutil.rmtree(self.temp_dir)

    @patch('sys.stdout', new_callable=io.StringIO)
    def test_main(self, mock_stdout):
        # Create a simple graph for testing
        input_graph_file = os.path.join(self.temp_dir, 'simple_graph.txt')
        save_graph(self.graph, input_graph_file, save_data=True)

        # Create a mock command-line argument list
        argv = ['graphpack',
                '--input', input_graph_file,
                '--output', self.temp_dir,
                '--method', 'greedy',
                '--is-lossless',
                '--is-gene-network']

        # Patch sys.argv to use our mock argument list
        with patch('sys.argv', argv):
            args = parse_args()
            perform_compression(**vars(args))

        # Capture the printed output
        printed_output = mock_stdout.getvalue()

        # Check if the output is as expected
        self.assertIn(input_graph_file, printed_output)
        self.assertIn('/tmp', printed_output)
        self.assertIn('greedy', printed_output)

    @patch('sys.stdout', new_callable=io.StringIO)
    def test_main_with_invalid_input_file(self, mock_stdout):
        # Mock read_graph to raise FileNotFoundError
        with patch('graphpack.utils.read_graph') as mock_read_graph:
            mock_read_graph.side_effect = FileNotFoundError("File not found")

            # Call main function with invalid input file
            with self.assertRaises(FileNotFoundError):
                with patch('sys.argv', ['graphpack', '--input', 'test.txt', '--output', self.temp_dir]):
                    with suppress_output():
                        args = parse_args()
                        perform_compression(**vars(args))


if __name__ == '__main__':
    unittest.main()
