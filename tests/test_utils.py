import json
import os
import shutil
import tempfile
import unittest

import networkx as nx

from graphpack.utils import read_graph, save_graph, draw_graph


class TestUtilsFunctions(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory for testing
        self.temp_dir = tempfile.mkdtemp()

        # Create a sample graph for testing
        self.graph = nx.Graph()
        self.graph.add_edge("A", "B", weight=1)
        self.graph.add_edge("B", "C", weight=2)
        self.graph.add_edge("A", "C", weight=2)
        self.graph.add_edge("D", "A", weight=4)

    def tearDown(self):
        # Remove the temporary directory
        shutil.rmtree(self.temp_dir)

    def test_draw_graph(self):
        # Test drawing the graph
        draw_graph(self.graph, title="test", is_interactive=False)

    def test_read_graph_from_edge_list(self):
        # Create a sample edge list file
        edge_list_path = os.path.join(self.temp_dir, 'test.edgelist')
        nx.write_edgelist(self.graph, edge_list_path)

        # Read the graph from the edge list file
        loaded_graph = read_graph(edge_list_path)

        # Assert that the loaded graph is the same as the original graph
        self.assertTrue(nx.is_isomorphic(self.graph, loaded_graph))

    def test_read_graph_from_json(self):
        # Create a sample JSON file
        json_path = os.path.join(self.temp_dir, 'test.json')
        data = {
            'nodes': ["A", "B", "C", "D"],
            'edges': [{'source': "A", 'target': "B"},
                      {'source': "B", 'target': "C", 'weight': 2},
                      {'source': "A", 'target': "C", 'weight': 3},
                      {'source': "D", 'target': "A", 'weight': 4},
                      ]
        }
        with open(json_path, 'w') as f:
            json.dump(data, f)

        # Read the graph from the JSON file
        loaded_graph = read_graph(json_path)

        # Assert that the loaded graph is the same as the original graph
        self.assertTrue(nx.is_isomorphic(self.graph, loaded_graph))

    def test_save_graph_as_edge_list(self):
        # Save the sample graph as an edge list file
        edge_list_path = os.path.join(self.temp_dir, 'test.edgelist')
        save_graph(self.graph, edge_list_path, save_data=True)

        # Read the saved graph from the edge list file
        loaded_graph = read_graph(edge_list_path)

        # Assert that the loaded graph is the same as the original graph
        self.assertTrue(nx.is_isomorphic(self.graph, loaded_graph))

        # Same, with CSV and TSV
        edge_list_path = os.path.join(self.temp_dir, 'test.csv')
        save_graph(self.graph, edge_list_path, save_data=True)
        loaded_graph = read_graph(edge_list_path)
        self.assertTrue(nx.is_isomorphic(self.graph, loaded_graph))

        edge_list_path = os.path.join(self.temp_dir, 'test.tsv')
        save_graph(self.graph, edge_list_path, save_data=True)
        loaded_graph = read_graph(edge_list_path)
        self.assertTrue(nx.is_isomorphic(self.graph, loaded_graph))

    def test_save_graph_as_csv(self):
        # Save the sample graph as an edge list file
        edge_list_path = os.path.join(self.temp_dir, 'test.csv')
        save_graph(self.graph, edge_list_path, save_data=True)

        # Read the saved graph from the edge list file
        loaded_graph = read_graph(edge_list_path)

        # Assert that the loaded graph has the same nodes and edges
        self.assertTrue(nx.is_isomorphic(self.graph, loaded_graph))

        # If the graph has edge attributes, check if they are preserved
        if self.graph.number_of_edges() > 0:
            for u, v, data in self.graph.edges(data=True):
                self.assertIn((u, v), loaded_graph.edges())
                self.assertEqual(loaded_graph[u][v]['weight'], data.get('weight', 1))

        # If the graph has node attributes, check if they are preserved
        if self.graph.number_of_nodes() > 0:
            for node, data in self.graph.nodes(data=True):
                self.assertIn(node, loaded_graph.nodes())
                # Compare node attributes if present
                for key, value in data.items():
                    self.assertEqual(loaded_graph.nodes[node].get(key), value)

    def test_save_graph_as_json(self):
        # Save the sample graph as a JSON file
        json_path = os.path.join(self.temp_dir, 'test.json')
        save_graph(self.graph, json_path, save_data=True)

        # Read the saved graph from the JSON file
        with open(json_path, 'r') as f:
            loaded_data = json.load(f)
        loaded_graph = nx.Graph()
        loaded_graph.add_nodes_from(loaded_data['nodes'])
        for edge in loaded_data['edges']:
            loaded_graph.add_edge(edge['source'], edge['target'], weight=edge.get('weight', 1))

        # Assert that the loaded graph is the same as the original graph
        self.assertTrue(nx.is_isomorphic(self.graph, loaded_graph))

    def test_invalid_file_extension(self):
        # Attempt to read from a file with an unsupported extension
        unsupported_file_path = os.path.join(self.temp_dir, 'test.invalid')
        with self.assertRaises(ValueError):
            read_graph(unsupported_file_path)

    def test_invalid_json_format(self):
        # Create a sample JSON file with invalid format
        invalid_json_path = os.path.join(self.temp_dir, 'invalid.json')
        with open(invalid_json_path, 'w') as f:
            f.write("This is not valid JSON")

        # Attempt to read from the invalid JSON file
        with self.assertRaises(json.JSONDecodeError):
            read_graph(invalid_json_path)

    def test_missing_edges_key_in_json(self):
        # Create a sample JSON file missing the 'edges' key
        missing_edges_json_path = os.path.join(self.temp_dir, 'missing_edges.json')
        data = {
            'nodes': ["A", "B", "C", "D"]
        }
        with open(missing_edges_json_path, 'w') as f:
            json.dump(data, f)

        # Attempt to read from the JSON file with missing 'edges' key
        with self.assertRaises(KeyError):
            read_graph(missing_edges_json_path)


if __name__ == '__main__':
    unittest.main()
