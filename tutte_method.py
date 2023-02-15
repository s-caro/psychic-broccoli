import matplotlib
import networkx as nx
from dimod import ConstrainedQuadraticModel, Binary, quicksum, Real, SampleSet, Integer
from dwave.system import LeapHybridCQMSampler

class Variables:
    """Class that collects all CQM model variables for the Tutte's barycenter method

    Args:
        nodes: the nodes of the graph considered
        lower_bound: the minumum value for the coordinate of a node
        upper_bound: the maximum value for the coordinate of a node
    """

    def __init__(self, nodes: nx.Graph.nodes, lowerBound: float, upperBound: float):
        self.x = {i: Integer(f'x_{i}',lower_bound = lowerBound, upper_bound = upperBound) for i in nodes}
        self.y = {i: Integer(f'y_{i}',lower_bound = lowerBound, upper_bound = upperBound) for i in nodes}


def build_graph(num_nodes):
    """Build graph."""

    print("\nBuilding graph...")

    G = nx.powerlaw_cluster_graph(num_nodes, 2, 0.4)

    return G

def _add_variable(cqm: ConstrainedQuadraticModel, lowerBound: float, upperBound: float, g: nx.Graph()):
    for i in G.nodes():
        cqm.add_variable('INTEGER',f'x_{i}',lower_bound = lowerBound, upper_bound = upperBound)
        cqm.add_variable('INTEGER',f'y_{i}',lower_bound = lowerBound, upper_bound = upperBound)


def _add_constraint(cqm: ConstrainedQuadraticModel, vars: Variables, g: nx.Graph()):
    """constraint on the x and y position of the nodes, each node must be on the barycenter of its neighbour

    Args:
        cqm: the model
        vars: the x and y coordinates of the points
        g: the graph that we want to draw
    """
    for v in g.nodes():
        cqm.add_constraint((g.degree(v)*vars.x[v]-quicksum(vars.x[u] for u,v in g.edges()))==0,label=f'x_constraint_node_{v}')
        cqm.add_constraint((g.degree(v)*vars.y[v]-quicksum(vars.y[u] for u,v in g.edges()))==0,label=f'y_constraint_node_{v}')

def _define_objective(cqm: ConstrainedQuadraticModel, vars: Variables, g: nx.Graph()):
    """objective function of the problem, minimize the distance of each point from the barycenter position

    Args:
        cqm: the model
        vars: the x and y coordinates of the points
        g: the graph that we want to draw
    """
    
    x_obj_term = quicksum((g.degree(v)*vars.x[v]-quicksum(vars.x[u] for u,v in g.edges()))**2 for v in g.nodes())
    y_obj_term = quicksum((g.degree(v)*vars.y[v]-quicksum(vars.y[u] for u,v in g.edges()))**2 for v in g.nodes())
    print(x_obj_term)
    x_obj_coefficient = 1
    y_obj_coefficient = 1
    cqm.set_objective(x_obj_coefficient*x_obj_term + y_obj_coefficient*y_obj_term)

def build_cqm(vars: Variables, g: nx.Graph(), fixed_points: list) -> ConstrainedQuadraticModel:
    """objective function of the problem, minimize the distance of each point from the barycenter position

    Args:
        vars: the x and y coordinates of the points
        g: the graph that we want to draw
        fixed_points: a list of tuple, each contains the node, the x-coordinate, the y-coordinate where the node needs to be fixed
    
    Returns:
        A ''dimod.CQM'' object that defines the Tutte's barycenter method
    """
    cqm = ConstrainedQuadraticModel()
    #_add_variable(cqm, 0, 1, g)

    _add_constraint(cqm, vars, g)

    
    _define_objective(cqm, vars, g)
    cqm.substitute_self_loops()
    for el in fixed_points:
        cqm.fix_variable(f'x_{el[0]}',el[1])
        cqm.fix_variable(f'y_{el[0]}', el[2])
    return cqm

def call_solver(cqm: ConstrainedQuadraticModel) -> SampleSet:
    sampler = LeapHybridCQMSampler()
    res = sampler.sample_cqm(cqm, label="Tutte's barycenter method")
    return res.first

if __name__ == '__main__':

    num_nodes = 4
    G = build_graph(num_nodes)
    G = nx.from_edgelist([(0,1),(0,4),(0,3),(1,5),(1,2),(2,6),(2,3),(3,7),(4,5),(4,7),(5,6),(6,7)])
    #print(type(G.nodes()))
    upperBound = 3
    lowerBound = 0
    vars = Variables(G.nodes(), lowerBound, upperBound)
    
    fixed_points = [(0,0,3),(1,3,3),(2,3,0),(3,0,0)]
    cqm = build_cqm(vars, G, fixed_points)
    
    best_feasible = call_solver(cqm)
    #print(best_feasible.info)
    print(best_feasible)