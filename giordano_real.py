import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import dimod
from dimod import ConstrainedQuadraticModel, Binary, quicksum, Real, SampleSet, Integer
from dwave.system import LeapHybridCQMSampler

dimod.REAL_INTERACTIONS =True

class Variables:
    """Class that collects all CQM model variables for the Tutte's barycenter method

    Args:
        nodes: the nodes of the graph considered
        lower_bound: the minumum value for the coordinate of a node
        upper_bound: the maximum value for the coordinate of a node
    """

    def __init__(self, nodes: nx.Graph.nodes, lowerBound: float, upperBound: float):
        self.x = {i: Real(f'x_{i}',lower_bound = lowerBound, upper_bound = upperBound) for i in nodes}
        self.y = {i: Real(f'y_{i}',lower_bound = lowerBound, upper_bound = upperBound) for i in nodes}


def build_graph(num_nodes):
    """Build graph."""

    print("\nBuilding graph...")

    G = nx.powerlaw_cluster_graph(num_nodes, 2, 0.4)

    return G

def _add_variable(cqm: ConstrainedQuadraticModel, lowerBound: float, upperBound: float, g: nx.Graph()):
    for i in G.nodes():
        cqm.add_variable('REAL',f'x_{i}',lower_bound = lowerBound, upper_bound = upperBound)
        cqm.add_variable('REAL',f'y_{i}',lower_bound = lowerBound, upper_bound = upperBound)


def _add_constraint(cqm: ConstrainedQuadraticModel, vars: Variables, g: nx.Graph(), no_nodes: list):
    """constraint on the x and y position of the nodes, each node must be on the barycenter of its neighbour

    Args:
        cqm: the model
        vars: the x and y coordinates of the points
        g: the graph that we want to draw
    """
    nodes = list(set(g.nodes)-set(no_nodes))
    for v in nodes:
        cqm.add_constraint((g.degree(v)*vars.x[v]-quicksum(vars.x[u] for u in g.neighbors(v)))==0,label=f'x_constraint_node_{v}')
        cqm.add_constraint((g.degree(v)*vars.y[v]-quicksum(vars.y[u] for u in g.neighbors(v)))==0,label=f'y_constraint_node_{v}')

def _define_objective(cqm: ConstrainedQuadraticModel, vars: Variables, g: nx.Graph(), no_nodes: list):
    """objective function of the problem, minimize the distance of each point from the barycenter position

    Args:
        cqm: the model
        vars: the x and y coordinates of the points
        g: the graph that we want to draw
    """
    nodes = list(set(g.nodes)-set(no_nodes))
    x_obj_term = quicksum((g.degree(v)*vars.x[v]-quicksum(vars.x[u] for u in g.neighbors(v))) for v in nodes)
    y_obj_term = quicksum((g.degree(v)*vars.y[v]-quicksum(vars.y[u] for u in g.neighbors(v))) for v in nodes)
    x_obj_coefficient = 1
    y_obj_coefficient = 1
    quadratic_obj_fun = (x_obj_coefficient*x_obj_term + y_obj_coefficient*y_obj_term)**2
    cqm.set_objective(quadratic_obj_fun)

def build_cqm(vars: Variables, g: nx.Graph(), fixed_points: list, upperBound: int) -> ConstrainedQuadraticModel:
    """objective function of the problem, minimize the distance of each point from the barycenter position

    Args:
        vars: the x and y coordinates of the points
        g: the graph that we want to draw
        fixed_points: a list of tuple, each contains the node, the x-coordinate, the y-coordinate where the node needs to be fixed
    
    Returns:
        A ''dimod.CQM'' object that defines the Titto's barycenter method
    """
    cqm = ConstrainedQuadraticModel()
    _add_variable(cqm, 0, upperBound, g)
    #print(vars)
    no_nodes = []
    for el in fixed_points:
        #print(el)
        no_nodes.append(el[0])
        cqm.add_constraint(vars.x[el[0]] == el[1],label=f'x_constraint_node_{el[0]}')
        cqm.add_constraint(vars.y[el[0]] == el[2], label=f'y_constraint_node_{el[0]}')
    #print(cqm.variables)
    #print(cqm.constraints)
    #_define_objective(cqm, vars, g, no_nodes)
    #print(cqm.variables)
    _add_constraint(cqm, vars, g, no_nodes)
    
    #cqm.substitute_self_loops()
    
    print(cqm.objective)
    #print(cqm.constraints)
    
    return cqm

def call_solver(cqm: ConstrainedQuadraticModel) -> SampleSet:
    sampler = LeapHybridCQMSampler()
    #4*sampler.min_time_limit(cqm),
    res = sampler.sample_cqm(cqm, label="Titto's barycenter method")
    return res.first

def create_pos_for_drawing(node_to_pos : dict) -> dict:
    res = {}
    print(node_to_pos)
    for coor in node_to_pos.keys():
        num = coor.split('_')[1]
        num = int(num)
        if num in res.keys():
            res[num].append(node_to_pos[coor])
        else:
            res[num] = [node_to_pos[coor]]
    
    print(res)
    return res

if __name__ == '__main__':

    num_nodes = 4
    G = build_graph(num_nodes)
    #G = nx.from_edgelist([(0,1),(1,2),(2,0),(0,3),(1,3),(2,3)])
    #G = nx.from_edgelist([(0,1),(0,4),(0,3),(1,5),(1,2),(2,3),(2,6),(3,7),(4,5),(4,7),(5,6),(6,7)])
    #G = nx.from_edgelist([(0,1),(0,2),(0,3),(0,4),(0,6),(1,3),(1,4),(1,2),(1,5),(2,3),(2,5),(2,6),(3,4),(3,5),(3,6)])
    #G =nx.from_edgelist([(0,1),(0,7),(0,4),(1,9),(1,2),(2,3),(2,19),(3,4),(3,17),(4,5),(5,6),(5,16),(6,7),(6,13),(7,8),(8,9),(8,12),(9,10),(10,11),(10,19),(11,12),(11,15),(12,13),(13,14),(14,15),(14,16),(15,18),(16,17),(17,18),(18,19)])
    #G = nx.from_edgelist([(0,1),(0,2),(0,3),(0,4),(0,5),(0,6),(1,2),(1,3),(2,3),(2,4),(2,5),(2,6),(3,4),(4,5),(5,6)])
    G=nx.from_edgelist([(0,1),(0,10),(0,9),(1,2),(1,11),(2,3),(2,11),(3,4),(3,12),(4,5),(4,12),(5,6),(5,13),(6,7),(6,13),(7,8),(7,14),(8,9),(8,14),(9,10),(10,15),(11,17),(12,18),(13,19),(14,16),(15,16),(15,17),(16,19),(17,18),(18,19)])
    upperBound = 100
    lowerBound = 0
    vars = Variables(G.nodes(), lowerBound, upperBound)
    
    #fixed_points = [(0,0,0),(1,2,6),(2,4,0)]
    #fixed_points = [(0,0,0),(1,0,6),(2,6,6),(3,6,0)]
    #fixed_points = [(0,0,0),(1,6,10),(2,12,0)]
    #fixed_points = [(0,0,45),(1,50,100),(2,100,45),(3,80,0),(4,20,0)]
    #fixed_points = [(0,0,0),(1,50,100),(2,100,0)]
    fixed_points = [(0,0,50),(1,10,80),(2,40,100),(3,60,100),(4,90,80),(5,100,50),(6,90,20),(7,60,0),(8,40,0),(9,10,20)]
    cqm = build_cqm(vars, G, fixed_points, upperBound)
    
    best_feasible = call_solver(cqm)
    print(best_feasible.info)
    #print(best_feasible)
    #best_feasible = {'x_0': 0.0, 'x_1': 6.0, 'x_2': 12.0, 'x_3': 6.0, 'x_4': 4.0, 'x_5': 8.0, 'x_6': 6.0, 'y_0': 0.0, 'y_1': 10.0, 'y_2': 0.0, 'y_3': 3.0, 'y_4': 4.0, 'y_5': 4.0, 'y_6': 1.0}
    pos = create_pos_for_drawing(best_feasible[0])
    
    nx.draw(G, pos=pos, node_size=300, edgecolors='k', cmap='hsv', with_labels=True)
    #print(G.edges())s
    plt.savefig("titto_draw_real_degree_3_decagono_no_obj_larger.png")
