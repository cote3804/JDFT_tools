import networkx as nx
import numpy as np
import matplotlib.pyplot as plt


def draw_graph(graph):
    # node_labels = {0:'1', 1:'2', 2:'3', 3:'4', 4:'5'}
    # nx.set_node_attributes(G, node_labels, 'label')
    pos = nx.spring_layout(graph, seed=2)
    edge_labels = nx.get_edge_attributes(graph, "weight")
    nx.draw_networkx_edge_labels(graph, pos, edge_labels)
    nx.draw_networkx_edges(graph, pos, width=1)
    nx.draw_networkx_nodes(graph, pos, node_size=200)
    nx.draw_networkx_labels(graph, pos)
    plt.show()

''' ##### Algorithm for calculating the cycle span of a graph given a vector of energies #####
# G = nx.DiGraph()

E = np.array([0, 0.5, 0.6, 0.3, 1]).T


# DiA = np.array([[0,1,0,0,0],
#              [0,0,1,1,0],
#              [0,0,0,0,1],
#              [0,0,0,0,1],
#              [1,0,0,0,0]])

# G = nx.from_numpy_array(DiA, create_using=nx.DiGraph())

G = nx.DiGraph()
G.add_nodes_from([1,2,3,4,5])
G.add_edges_from([(1,2), (2,3), (2,4), (3,5), (4,5), (5,1)])
# print(nx.adjacency_matrix(G).todense())
I = nx.incidence_matrix(G, oriented=True).todense()
E_e = np.matmul(I.T, E) #edge energy weights
E_e[-1] = 0

def add_edge_weights(G:nx.DiGraph, edge_weights:np.array):
    for i, (u,v) in enumerate(G.edges()):
        G[u][v]['weight'] = edge_weights[i]
    return G



def cycle_spans(graph):
    cycles = nx.simple_cycles(graph)
    span_list = []
    for cycle in cycles:
        subgraph = graph.subgraph(cycle)
        FW = nx.floyd_warshall_numpy(subgraph)
        print(FW, np.argmax(FW))

G = add_edge_weights(G, E_e)

cycles = nx.simple_cycles(G)
cycles_list = sorted(cycles)
subgraph = G.subgraph(cycles_list[1])
FW = nx.floyd_warshall_numpy(subgraph)

cycle_spans(G)

# print(FW)
draw_graph(G)
'''

#### Algorithm for extracting subgraph of complete graph ####

def reconnect(G, subgraph):
    # master graph is G
    # subgraph needs to reconnected following master graph 
    A_s = nx.adjacency_matrix(subgraph).todense()
    A_m = nx.adjacency_matrix(G).todense()
    new_array = np.zeros(A_m.shape)
    new_array[np.ix_([0, 2, 3, 4], [0,2,3,4])] = A_s
    mismatch = np.equal(new_array, A_m).astype(int) # slice to exclude 2nd row and 2nd column from new_array
    mismatch_indices = np.where(mismatch == 0)
    print(mismatch_indices)
    pairs_to_add = []
    for x,y in zip(mismatch_indices[0], mismatch_indices[1]):
        target_indices = np.where(mismatch[y,:] == 0)
        node_pairs = [(x, yy) for yy in target_indices[0]]
        pairs_to_add.append(node_pairs)

    print(pairs_to_add[0][0])
    subgraph.add_edge(pairs_to_add[0][0][0], pairs_to_add[0][0][1])
    subgraph.add_edge(pairs_to_add[0][1][0], pairs_to_add[0][1][1])
    # for pair in pairs_to_add:
    #     print(pair)
    #     subgraph.add_edge(pair)


    return subgraph
    # for node in subgraph.nodes():
    #     sub_connections = A_s[node,:]
    #     master_connections = A_m[node,:]
    #     missing_connections = np.where(sub_connections != master_connections)[0]
    #     print(missing_connections)

G = nx.DiGraph()
G.add_nodes_from([0,1,2,3,4])
G.add_edges_from([(0,1), (1,2), (1,3), (2,4), (3,4), (4,0)])
# print(nx.adjacency_matrix(G).todense())
subgraph = G.subgraph([0,2,3,4])
unfrozen_subgraph = nx.DiGraph(subgraph)
subgraph = reconnect(G, unfrozen_subgraph)
draw_graph(G)
draw_graph(subgraph)