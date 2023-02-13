from dimod import Binary, ConstrainedQuadraticModel, quicksum, BinaryQuadraticModel
from dwave.system import LeapHybridCQMSampler, DWaveSampler, EmbeddingComposite
import networkx as nx

import matplotlib

try:
    import matplotlib.pyplot as plt
except ImportError:
    matplotlib.use("agg")
    import matplotlib.pyplot as plt

def build_graph(num_nodes):
    """Build graph."""

    print("\nBuilding graph...")

    G = nx.powerlaw_cluster_graph(num_nodes, 3, 0.4)
    
    G = nx.Graph([(0, 3), (0, 4), (0, 6), (0, 7), (1, 3), (1, 4), (1, 5), (2, 3), (2, 5), (3, 4), (3, 6), (3, 7), (4, 5), (4, 6), (6, 7)])
    pos = nx.spring_layout(G)
    nx.draw(G, pos=pos, node_size=50, edgecolors='k', cmap='hsv')
    plt.savefig("original_graph.png")

    print(G.edges)

    return G, pos

# function to build the Constrained Quadratic Model, it takes in input the graph G and the number of partitions, it returns the cqm for the problem
def build_cqm(G,k):
    
    # create array of the partitions, each partition is identified by the index in the array
    partitions = range(k)
    
    # Initialize the CQM object
    cqm = ConstrainedQuadraticModel()

    # Add binary variable, the binary variable x_(i,k) is set ot 1 if the node i is assigned to partition k
    v =[[Binary(f'v_{i},{p}') for p in partitions] for i in G.nodes]

    # Constraint 1: each node is assigned only to one partition
    #for i in G.nodes:
    #    cqm.add_constraint(quicksum(v[i][k] for k in partitions) == 1, label=f'one-hot-node-{i}')

    # is equivalent to this? this is the original from dwave
    for i in G.nodes:
        cqm.add_discrete([f'v_{i},{k}' for k in partitions], label=f'one-hot-node-{i}')

    for p in partitions:
        #print("\nAdding partition size constraint for partition", k)
        cqm.add_constraint(quicksum(v[n][p] for n in G.nodes) == G.number_of_nodes()/k, label='partition-size-{}'.format(p))

    # Objective function: minimize edges between partitions
    min_edges = []
    for i,j in G.edges:
        for p in partitions:
            min_edges.append(v[i][p]+v[j][p]-2*v[i][p]*v[j][p])
    cqm.set_objective(sum(min_edges))
    return cqm

def run_resolution(cqm, sampler_cqm,):
    """Send the CQM to the sampler and return the best sample found.
    Args:
        cqm (ConstrainedQuadraticModel): The CQM for our problem
        sampler: The CQM sampler to be used. Must have sample_cqm function.
    
    Returns:
        dict: The first feasible solution found
    """

    # Initialize the solver
    print("\nSending to the cqm solver...")
    
    # Solve the CQM problem using the solver
    sampleset_cqm = sampler_cqm.sample_cqm(cqm, 4*sampler_cqm.min_time_limit(cqm), label='Example - Graph Partitioning cqm')

    # Return the first feasible solution

    return sampleset_cqm





# ------- Main program -------
if __name__ == "__main__":

    num_nodes = 8

    G, pos = build_graph(num_nodes)
    print("\nCreating cqm model...")
    cqm = build_cqm(G, 2)

    sampler_cqm = LeapHybridCQMSampler()

    
    sample_cqm = run_resolution(cqm, sampler_cqm)

    print(sample_cqm.first)


