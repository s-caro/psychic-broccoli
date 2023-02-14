from dimod import Binary, ConstrainedQuadraticModel, quicksum, BinaryQuadraticModel
from dwave.system import LeapHybridCQMSampler, DWaveSampler, EmbeddingComposite
import networkx as nx
from dwave.cloud import Client

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
    
    pos = nx.spring_layout(G)
    nx.draw(G, pos=pos, node_size=300, edgecolors='k', cmap='hsv', with_labels=True)
    print(G.edges())
    plt.savefig("original_graph_partitioning.png")

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

    print("qpu_access_time: ", sampleset_cqm.info['qpu_access_time'])
    print("run_time: ", sampleset_cqm.info['run_time'])
    # Return the first feasible solution
    soln = {key[0]: key[1] for key, val in sampleset_cqm.first.sample.items() if val == 1.0}

    return soln

def plot_soln(sample, pos):
    """Plot results and save file.
    
    Args:
        sample (dict):
            Sample containing a solution. Each key is a node and each value 
            is an int representing the node's color.
        pos (dict):
            Plotting information for graph so that same graph shape is used.
    """

    print("\nProcessing sample...")

    node_colors_num = [sample[i] for i in G.nodes()]
    numToColor = {0: 'yellow', 1: 'blue', 2: 'red', 3:'green'}
    node_colors = []
    for i in node_colors_num:
        node_colors.append(numToColor[i])
    print(node_colors_num)
    print(node_colors)
    nx.draw(G, pos=pos, node_color=node_colors, node_size=300, edgecolors='k', cmap='hsv',with_labels=True)
    fname = 'graph_result_partitioning.png'
    plt.savefig(fname)

    print("\nSaving results in {}...".format(fname))



# ------- Main program -------
if __name__ == "__main__":

    num_nodes = 20

    G, pos= build_graph(num_nodes)
    print("\nCreating cqm model...")
    cqm = build_cqm(G, 4)

    sampler_cqm = LeapHybridCQMSampler()

    
    sample_cqm = run_resolution(cqm, sampler_cqm)
    

    plot_soln(sample_cqm, pos)
