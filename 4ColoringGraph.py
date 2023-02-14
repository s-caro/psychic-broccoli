import matplotlib
import networkx as nx
from dimod import ConstrainedQuadraticModel, Binary, quicksum
from dwave.system import LeapHybridCQMSampler

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
    #[(0, 3), (0, 4), (0, 6), (0, 7), (0, 8), (0, 10), (0, 12), (0, 13), (1, 3), (1, 13), (2, 3), (2, 4), (2, 5), (2, 7), (2, 8), (2, 9), (2, 11), (2, 14), (3, 4), (3, 5), (3, 6), (3, 10), (3, 11), (3, 13), (4, 5), (4, 6), (4, 7), (4, 8), (4, 9), (6, 10), (7, 12), (7, 14), (8, 9), (8, 14), (9, 11), (11, 12)]
    nx.draw(G, pos=pos, node_size=300, edgecolors='k', cmap='hsv', with_labels=True)
    print(G.edges())
    plt.savefig("original_graph_coloring.png")

    return G, pos

def build_cqm(G, num_colors):
    """Build CQM model."""

    print("\nBuilding constrained quadratic model...")

    # Initialize the CQM object
    cqm = ConstrainedQuadraticModel()

    # Build CQM variables
    colors = {n: {c: Binary((n, c)) for c in range(num_colors)} for n in G.nodes}

    # Add constraint to make variables discrete
    for n in G.nodes():
        cqm.add_discrete([(n, c) for c in range(num_colors)])
  
    # Build the constraints: edges have different color end points
    for u, v in G.edges:
        for c in range(num_colors):
            cqm.add_constraint(colors[u][c]*colors[v][c] == 0)

    return cqm

def run_hybrid_solver(cqm):
    """Solve CQM using hybrid solver."""

    print("\nRunning hybrid sampler...")

    # Initialize the CQM solver
    sampler = LeapHybridCQMSampler()

    # Solve the problem using the CQM solver
    sampleset = sampler.sample_cqm(cqm, label='Example - Graph Coloring')
    feasible_sampleset = sampleset.filter(lambda row: row.is_feasible)

    print("qpu_access_time: ", feasible_sampleset.info['qpu_access_time'])
    print("run_time: ", feasible_sampleset.info['run_time'])

    print(feasible_sampleset.first)



    try:
        sample = feasible_sampleset.first.sample
    except:
        print("\nNo feasible solutions found.")
        exit()

    soln = {key[0]: key[1] for key, val in sample.items() if val == 1.0}

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
    fname = 'graph_result_coloring.png'
    plt.savefig(fname)

    print("\nSaving results in {}...".format(fname))

# ------- Main program -------
if __name__ == "__main__":

    num_nodes = 15

    G, pos = build_graph(num_nodes)
    num_colors = 4
    
    cqm = build_cqm(G, num_colors)

    sample = run_hybrid_solver(cqm)
    

    plot_soln(sample, pos)