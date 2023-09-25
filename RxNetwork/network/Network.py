import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

class Network():
    def __init__(self, reaction):
        # network_definition is the dictionary that encodes the connectivity and labels of the network
        # The node index is the first value for each state, and the second value is a list of connections
        # connections are drawn from the current node to the next one.
        self.network_definition = {"NRR": {"initial":(0, [(0,1)]),
                                      "N2": (1, [(1,2)]),
                                      "N2H": (2, [(2,3), (2,9)]),
                                      "NNH2": (3, [(3,4)]),
                                      "N": (4, [(4,5)]),
                                      "NH": (5, [(5,6)]),
                                      "NH2": (6, [(6,7)]),
                                      "NH3": (7, [(7,8)]),
                                      "final": (8, [(8,0)]),
                                      "NHNH": (9, [(9,10)]),
                                      "NHNH2": (10, [(10,11)]),
                                      "NH2NH2": (11, [(11,6)]),}
                                      }
        self.reaction = reaction

        def build_master_graph(network_definition, reaction):
            # Builds the network x digraph (directed graph) object given the network definition
            graph = nx.DiGraph()
            for state_string, state_data in network_definition[reaction].items():
                graph.add_node(state_data[0], state=state_string)
                for edge in state_data[1]:
                    graph.add_edge(edge[0], edge[1])
            return graph
        
        self.master_graph = build_master_graph(self.network_definition, reaction)
        # print(nx.adjacency_matrix(self.master_graph).todense())
        node_list = list(range(0, len(self.master_graph.nodes())))
        print(self.master_graph.nodes())
        # self.draw_graph(self.master_graph)

    def add_data_to_nodes(self, FED_energies:dict):
        '''
        Adds a dictionary to each state that looks like this:
        {index: {"label":intermediate, "energy":energy}}
        for example, the starting node would look like this:
        {0: {"label":"initial", "energy":1.9}}  (This assumes NRR)
        '''
        for state, energy in FED_energies.items():
            index = self.state_to_index(state)
            self.master_graph.nodes[index].update({"label":state, "energy":energy})

    def states_to_indices(self, states:list) -> list:
        # Converts a list of states to a list of indices
        indices = []
        for state in states:
            indices.append(self.state_to_index(state))
        return indices
    
    def state_to_index(self, state:str) -> int:
        return self.network_definition[self.reaction][state][0]
    
    def get_loop_indices(self, states) -> list:
        """
        This functiopn takes a list of states and returns a list of their indices within their respective loops
        This is used to plot multiple pathways on the same FED plot
        """
        for cycle in nx.simple_cycles(self.master_graph):
            print(cycle)
            remap = {} # dictionary to remap indices of the cycle
            for old_index in cycle:
                remap[old_index] = None
            
            nx.relabel_nodes(cycle ,remap, copy=False)
            self.draw_graph(cycle)
        return None


    def connected_subgraph(self, intermediates:list) -> nx.DiGraph:
        # Returns a subgraph of the master graph that contains only the intermediates
        # that have been converged
        # Employs cooper's graph algorithm to draw connections to dangling nodes
        indices = self.states_to_indices(intermediates)
        indices.extend([0, self.network_definition[self.reaction]["final"][0]]) #always need initial and final states
        subgraph = self.master_graph.subgraph(indices)
        subgraph = self.reconnect(subgraph)
        # self.draw_graph(subgraph)
        # return subgraph

    def reconnect(self, sub_graph):
        sub_graph = sub_graph.copy()
        node_list = list(range(0, len(self.master_graph.nodes()))) # need to do this or networkx will return the wrong adjacency matrix
        A_master = nx.adjacency_matrix(self.master_graph, nodelist=node_list).todense()
        # odd_graph = nx.DiGraph(A_master)
        # self.draw_graph(odd_graph)
        print(A_master)
        sub_node_list = list(sub_graph.nodes()).sort()
        A_sub = nx.adjacency_matrix(sub_graph, nodelist=sub_node_list).todense()
        odd_graph = nx.DiGraph(A_sub)
        # self.draw_graph(odd_graph)
        print(A_sub)
        subgraph_node_indices = list(sub_graph.nodes())
        subgraph_node_indices.sort()
        print(subgraph_node_indices)
        # new_array is an array full of zeros of the same size as the master graph adjacency matrix
        new_array = np.zeros(A_master.shape)
        # here I'm taking 2D a slice of new_array that corresponds to the indices of the subgraph
        # I then replace that slice with the subgraph adjacency matrix
        # This effectively creates an adjacency matrix that has the same size as the master graph, but with
        # zeros in the places where there are no nodes in the sub graph.
        new_array[np.ix_(subgraph_node_indices, subgraph_node_indices)] = A_sub
        print(new_array)
        # here I'm creating a matrix that shows the differences between the adjacency matrix of the master graph and the new array
        # This will show where the subgraph is missing connections, enabling us to redraw them.
        mismatch = np.equal(new_array, A_master).astype(int)
        print(mismatch)
        # I now get the indices of the places where the mismatch matrix is zero, signifying a missing connection
        mismatch_indices = np.where(mismatch == 0)
        print("mismatch indices", mismatch_indices)
        pairs_to_add = []
        for x,y in zip(mismatch_indices[0], mismatch_indices[1]): 
            target_indices = np.where(mismatch[y,:] == 0)
            node_pairs = [(x, yy) for yy in target_indices[0]]
            pairs_to_add.append(node_pairs)
        print(pairs_to_add)
        for pairs in pairs_to_add:
            for pair in pairs:
                sub_graph.add_edge(pair[0], pair[1])

        return sub_graph


    def draw_graph(self, graph):
        # node_labels = {0:'1', 1:'2', 2:'3', 3:'4', 4:'5'}
        # nx.set_node_attributes(G, node_labels, 'label')
        pos = nx.spring_layout(graph, seed=2)
        edge_labels = nx.get_edge_attributes(graph, "weight")
        nx.draw_networkx_edge_labels(graph, pos, edge_labels)
        nx.draw_networkx_edges(graph, pos, width=1)
        nx.draw_networkx_nodes(graph, pos, node_size=200)
        nx.draw_networkx_labels(graph, pos)
        plt.show()

        
