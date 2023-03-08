import matplotlib
import networkx as nx
from dimod import ConstrainedQuadraticModel, Binary, quicksum, cqm_to_bqm, BinaryQuadraticModel
from dwave.system import LeapHybridCQMSampler, DWaveSampler, EmbeddingComposite
import dimod
import dwave.inspector

import dwave.inspector

try:
    import matplotlib.pyplot as plt
except ImportError:
    matplotlib.use("agg")
    import matplotlib.pyplot as plt

def build_graph(num_nodes):
    """Build graph."""

    print("\nBuilding graph...")

    G = nx.powerlaw_cluster_graph(num_nodes, 2, 0.4)
    if(nx.is_planar(G)):
        pos = nx.planar_layout(G)
    else:
        pos = nx.spring_layout(G)
    #[(0, 3), (0, 4), (0, 6), (0, 7), (0, 8), (0, 10), (0, 12), (0, 13), (1, 3), (1, 13), (2, 3), (2, 4), (2, 5), (2, 7), (2, 8), (2, 9), (2, 11), (2, 14), (3, 4), (3, 5), (3, 6), (3, 10), (3, 11), (3, 13), (4, 5), (4, 6), (4, 7), (4, 8), (4, 9), (6, 10), (7, 12), (7, 14), (8, 9), (8, 14), (9, 11), (11, 12)]

    nx.draw(G, pos=pos, node_size=300, edgecolors='k', cmap='hsv', with_labels=True)
    #print(G.edges())
    plt.savefig("original_graph_coloring.png")

    return G, pos


    

def build_cqm(G, num_colors):
    """Build CQM model."""

    print("\nBuilding constrained quadratic model...")

    # Initialize the CQM object
    cqm = ConstrainedQuadraticModel()

    # Build CQM variables
    colors = {n: {c: Binary((n, c)) for c in range(num_colors)} for n in G.nodes}

    # Add constraint each node can have only one color
    for n in G.nodes():
        cqm.add_discrete([(n, c) for c in range(num_colors)], label = f'1_color_node_{n}')
  
    # Build the constraints: edges have different color end points
    for u, v in G.edges:
        for c in range(num_colors):
            cqm.add_constraint(colors[u][c]+colors[v][c] == 1, label = f'no_adjacent_node_{u},{v}_color_{c}')

    return cqm

def run_hybrid_solver(cqm):
    """Solve CQM using hybrid solver."""

    print("\nRunning hybrid sampler...")

    # Initialize the CQM solver
    sampler = LeapHybridCQMSampler()

    # Solve the problem using the CQM solver
    sampleset = sampler.sample_cqm(cqm, label='Example - Graph Coloring')
    feasible_sampleset = sampleset.filter(lambda row: row.is_feasible)

    print("qpu_access_time: ", feasible_sampleset.info['qpu_access_time'], "us")
    print("run_time: ", feasible_sampleset.info['run_time'], "us")

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
    #print(node_colors_num)
    #print(node_colors)
    nx.draw(G, pos=pos, node_color=node_colors, node_size=300, edgecolors='k', cmap='hsv',with_labels=True)
    fname = 'graph_result_coloring.png'
    plt.savefig(fname)

    print("\nSaving results in {}...".format(fname))

def solve_and_print_bqm(G, num_colors ,pos):
    bqm = BinaryQuadraticModel('BINARY')

    x =[[f'x_{n}{c}' for c in range(num_colors)] for n in G.nodes]

    for n in G.nodes:
        for c in range(num_colors):
            bqm.add_variable(x[n][c])

    c1 = [(x[n][c],1) for c in range(num_colors) for n in G.nodes]
    bqm.add_linear_equality_constraint(c1, 10, -1)

    
    for c in range(num_colors):
        c3 = [(x[u][c] + x[v][c], 1) for u,v in G.edges]
        bqm.add_linear_inequality_constraint(c3, constant=-1, lagrange_multiplier=10, label=f'c3_color_{c}')
        
    sampler = EmbeddingComposite(DWaveSampler())
    sampleset = sampler.sample(bqm, num_reads=1000)


    print(sampleset.first)


    try:
        sample = sampleset.first.sample
    except:
        print("\nNo feasible solutions found.")
        exit()

    soln = {key[0]: key[1] for key, val in sample.items() if val == 1}

    print("\nProcessing sample...")

    node_colors_num = [soln[i] for i in G.nodes()]
    numToColor = {0: 'yellow', 1: 'blue', 2: 'red', 3:'green'}
    node_colors = []
    for i in node_colors_num:
        node_colors.append(numToColor[i])
    #print(node_colors_num)
    #print(node_colors)
    nx.draw(G, pos=pos, node_color=node_colors, node_size=300, edgecolors='k', cmap='hsv',with_labels=True)
    fname = 'graph_result_coloring_bqm.png'
    plt.savefig(fname)

    print("\nSaving results in {}...".format(fname))


# ------- Main program -------
if __name__ == "__main__":

    #num_nodes = 5

    #G, pos = build_graph(num_nodes)
    G = nx.from_edgelist([(0,1),(1,2),(2,3),(3,0)])
    num_colors = 2
    #solve_and_print_bqm(G, num_colors, pos)

   

    
    cqm = build_cqm(G, num_colors)

    print(cqm)
    bqm, invert = dimod.cqm_to_bqm(cqm)
    #sample = run_hybrid_solver(cqm)
    #print(bqm)
    sampleset2 = dimod.ExactSolver().sample(bqm)
    #print(sampleset.first.sample)
    print(invert(sampleset2.first.sample))

    Q = bqm.to_qubo()

    sampler = EmbeddingComposite(DWaveSampler(solver=dict(qpu=True)))

    sampleset = sampler.sample(bqm, num_reads=1000, label='Coloring example - Inspector')
    print(f'Energy: {sampleset.first.energy}')
    dwave.inspector.show(sampleset)
    #plot_soln(sample, pos)