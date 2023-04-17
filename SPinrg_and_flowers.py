import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import math 
import dimod
from dimod import ConstrainedQuadraticModel, Binary, quicksum, Real, SampleSet, Integer, ExactCQMSolver
from dwave.system import LeapHybridCQMSampler
import pprint

dimod.REAL_INTERACTIONS =True

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


def _add_variable(cqm: ConstrainedQuadraticModel, lowerBound: float, upperBound: float, g: nx.Graph()):
    for i in G.nodes():
        cqm.add_variable('INTEGER',f'x_{i}',lower_bound = lowerBound, upper_bound = upperBound)
        cqm.add_variable('INTEGER',f'y_{i}',lower_bound = lowerBound, upper_bound = upperBound)

def _add_constraint_alt(cqm: ConstrainedQuadraticModel, vars: Variables, g: nx.Graph(), l:float):
    l
    l_1 = 1/l
    print(cqm.variables[0])
    for u in g.nodes:
        nodes = list(set(g.nodes)-set([u]))
        F_rep_x = 0
        F_rep_y = 0
        F_att_x = 0
        F_att_y = 0
        for v in nodes:
            print(f'node {v}')
            print(vars.y[u])
            print(type(vars.x[u]-vars.x[v] >= 0))

            F_rep_x += (((vars.x[u]-vars.x[v])+(vars.y[u]-vars.y[v]))*(vars.x[v]-vars.x[u])*(l**2))
            F_rep_x += (((vars.x[v]-vars.x[u])+(vars.y[u]-vars.y[v]))*(vars.x[v]-vars.x[u])*(l**2))
            F_rep_x += (((vars.x[u]-vars.x[v])+(vars.y[v]-vars.y[u]))*(vars.x[v]-vars.x[u])*(l**2))
            F_rep_x += (((vars.x[v]-vars.x[u])+(vars.y[v]-vars.y[u]))*(vars.x[v]-vars.x[u])*(l**2))

            F_rep_y += (((vars.x[u]-vars.x[v])+(vars.y[u]-vars.y[v]))*(vars.y[v]-vars.y[u])*(l**2))
            F_rep_y += (((vars.x[v]-vars.x[u])+(vars.y[u]-vars.y[v]))*(vars.y[v]-vars.y[u])*(l**2))
            F_rep_y += (((vars.x[u]-vars.x[v])+(vars.y[v]-vars.y[u]))*(vars.y[v]-vars.y[u])*(l**2))
            F_rep_y += (((vars.x[v]-vars.x[u])+(vars.y[v]-vars.y[u]))*(vars.y[v]-vars.y[u])*(l**2))

        for v in g.neighbors(u):

            F_att_x += (((vars.x[u]-vars.x[v])+(vars.y[u]-vars.y[v]))*(vars.x[v]-vars.x[u])*l_1 )
            F_att_x += (((vars.x[u]-vars.x[v])+(vars.y[v]-vars.y[u]))*(vars.x[v]-vars.x[u])*l_1 )
            F_att_x += (((vars.x[v]-vars.x[u])+(vars.y[u]-vars.y[v]))*(vars.x[v]-vars.x[u])*l_1 )
            F_att_x += (((vars.x[v]-vars.x[u])+(vars.y[v]-vars.y[u]))*(vars.x[v]-vars.x[u])*l_1 )

            F_att_y += (((vars.x[u]-vars.x[v])+(vars.y[u]-vars.y[v]))*(vars.y[v]-vars.y[u])*l_1 )
            F_att_y += (((vars.x[u]-vars.x[v])+(vars.y[v]-vars.y[u]))*(vars.y[v]-vars.y[u])*l_1 )
            F_att_y += (((vars.x[v]-vars.x[u])+(vars.y[u]-vars.y[v]))*(vars.y[v]-vars.y[u])*l_1 )
            F_att_y += (((vars.x[v]-vars.x[u])+(vars.y[v]-vars.y[u]))*(vars.y[v]-vars.y[u])*l_1 )


        cqm.add_constraint((-F_rep_x + F_att_x)==0, label=f'x_constraint_node_{u}')
        cqm.add_constraint((-F_rep_y + F_att_y)==0, label=f'y_constraint_node_{u}')
            

def _add_constraint(cqm: ConstrainedQuadraticModel, vars: Variables, g: nx.Graph(), l: float):
    """constraint on the x and y position of the nodes, the sum of all forces acting on a point must be zero
    Args:
        cqm: the model
        vars: the x and y coordinates of the points
        g: the graph that we want to draw
    """
    l
    l_1 = 1/l
    for u in g.nodes:
        nodes = list(set(g.nodes)-set([u]))
        # non posso usare quadrati o valori assoluti per cui per ogni differenza devo cercare i vari casi
        F_rep_x = quicksum(((vars.x[u]-vars.x[v])+(vars.y[u]-vars.y[v]))*(vars.x[v]-vars.x[u])*(l**2) for v in nodes if vars.x[u]-vars.x[v] >= 0 and vars.y[u]-vars.y[v] >= 0)
        F_rep_x = quicksum(((vars.x[v]-vars.x[u])+(vars.y[u]-vars.y[v]))*(vars.x[v]-vars.x[u])*(l**2) for v in nodes if vars.x[u]-vars.x[v] <= -1 and vars.y[u]-vars.y[v] >= 0)
        F_rep_x = quicksum(((vars.x[u]-vars.x[v])+(vars.y[v]-vars.y[u]))*(vars.x[v]-vars.x[u])*(l**2) for v in nodes if vars.x[u]-vars.x[v] >= 0 and vars.y[u]-vars.y[v] <= -1)
        F_rep_x = quicksum(((vars.x[v]-vars.x[u])+(vars.y[v]-vars.y[u]))*(vars.x[v]-vars.x[u])*(l**2) for v in nodes if vars.x[u]-vars.x[v] <= -1 and vars.y[u]-vars.y[v] <= -1)


        F_att_x = quicksum(((vars.x[u]-vars.x[v])+(vars.y[u]-vars.y[v]))*(vars.x[v]-vars.x[u])*l_1 for v in g.neighbors(u) if vars.x[u]-vars.x[v] >= 0 and vars.y[u]-vars.y[v] >= 0)
        F_att_x = quicksum(((vars.x[u]-vars.x[v])+(vars.y[v]-vars.y[u]))*(vars.x[v]-vars.x[u])*l_1 for v in g.neighbors(u) if vars.x[u]-vars.x[v] >= 0 and vars.y[u]-vars.y[v] <= -1)
        F_att_x = quicksum(((vars.x[v]-vars.x[u])+(vars.y[u]-vars.y[v]))*(vars.x[v]-vars.x[u])*l_1 for v in g.neighbors(u) if vars.x[u]-vars.x[v] <= -1 and vars.y[u]-vars.y[v] >= 0)
        F_att_x = quicksum(((vars.x[v]-vars.x[u])+(vars.y[v]-vars.y[u]))*(vars.x[v]-vars.x[u])*l_1 for v in g.neighbors(u) if vars.x[u]-vars.x[v] <= -1 and vars.y[u]-vars.y[v] <= -1)

        F_rep_y = quicksum(((vars.x[u]-vars.x[v])+(vars.y[u]-vars.y[v]))*(vars.y[v]-vars.y[u])*(l**2) for v in nodes if vars.x[u]-vars.x[v] >= 0 and vars.y[u]-vars.y[v] >= 0)
        F_rep_y = quicksum(((vars.x[v]-vars.x[u])+(vars.y[u]-vars.y[v]))*(vars.y[v]-vars.y[u])*(l**2) for v in nodes if vars.x[u]-vars.x[v] <= -1 and vars.y[u]-vars.y[v] >= 0)
        F_rep_y = quicksum(((vars.x[u]-vars.x[v])+(vars.y[v]-vars.y[u]))*(vars.y[v]-vars.y[u])*(l**2) for v in nodes if vars.x[u]-vars.x[v] >= 0 and vars.y[u]-vars.y[v] <= -1)
        F_rep_y = quicksum(((vars.x[v]-vars.x[u])+(vars.y[v]-vars.y[u]))*(vars.y[v]-vars.y[u])*(l**2) for v in nodes if vars.x[u]-vars.x[v] <= -1 and vars.y[u]-vars.y[v] <= -1)


        F_att_y = quicksum(((vars.x[u]-vars.x[v])+(vars.y[u]-vars.y[v]))*(vars.y[v]-vars.y[u])*l_1 for v in g.neighbors(u) if vars.x[u]-vars.x[v] >= 0 and vars.y[u]-vars.y[v] >= 0)
        F_att_y = quicksum(((vars.x[u]-vars.x[v])+(vars.y[v]-vars.y[u]))*(vars.y[v]-vars.y[u])*l_1 for v in g.neighbors(u) if vars.x[u]-vars.x[v] >= 0 and vars.y[u]-vars.y[v] <= -1)
        F_att_y = quicksum(((vars.x[v]-vars.x[u])+(vars.y[u]-vars.y[v]))*(vars.y[v]-vars.y[u])*l_1 for v in g.neighbors(u) if vars.x[u]-vars.x[v] <= -1 and vars.y[u]-vars.y[v] >= 0)
        F_att_y = quicksum(((vars.x[v]-vars.x[u])+(vars.y[v]-vars.y[u]))*(vars.y[v]-vars.y[u])*l_1 for v in g.neighbors(u) if vars.x[u]-vars.x[v] <= -1 and vars.y[u]-vars.y[v] <= -1)


        cqm.add_constraint((-F_rep_x + F_att_x)==0, label=f'x_constraint_node_{u}')
        cqm.add_constraint((-F_rep_y + F_att_y)==0, label=f'y_constraint_node_{u}')
        #cqm.add_constraint((quicksum((-l**2)*(vars.y[v]-vars.y[u])*((vars.x[v]-vars.x[u])**2+(vars.y[v]-vars.y[u])**(2)) for v in nodes) + quicksum((math.sqrt((vars.x[u]-vars.x[v])**2+(vars.y[u]-vars.y[v])**2))*(vars.y[v]-vars.y[u])*(l**(-1)) for v in g.neighbors(u)))==0, label=f'y_constraint_node_{u}')

        #cqm.add_constraint((g.degree(v)*vars.x[v]-quicksum(vars.x[u] for u in g.neighbors(v)))==0,label=f'x_constraint_node_{u}')
        #cqm.add_constraint((g.degree(v)*vars.y[v]-quicksum(vars.y[u] for u in g.neighbors(v)))==0,label=f'y_constraint_node_{u}')

def build_cqm(vars: Variables, g: nx.Graph(), upperBound: int, C: int) -> ConstrainedQuadraticModel:
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
    l = C*math.sqrt(upperBound**2/len(list(g.nodes)))
    _add_constraint_alt(cqm, vars, g, l)
    #_add_objective(cqm, vars, g, l)
    

    
    return cqm

def call_solver(cqm: ConstrainedQuadraticModel) -> SampleSet:
    sampler = LeapHybridCQMSampler()
    #4*sampler.min_time_limit(cqm),
    res = sampler.sample_cqm(cqm,  label="Titto's flower method")
    return res.first

def create_pos_for_drawing(node_to_pos : dict) -> dict:
    res = {}
    #print(node_to_pos)
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
    #1
    G = nx.from_edgelist([(0,1),(1,2),(2,0),(0,3),(1,3),(2,3)])

    pos = nx.spring_layout(G)

    nx.draw(G, pos=pos, node_size=300, edgecolors='k', cmap='hsv', with_labels=True)
    #print(G.edges())s
    plt.savefig("nx_draw_spring_layout.png")

    #2
    #G = nx.from_edgelist([(0,1),(0,4),(0,3),(1,5),(1,2),(2,3),(2,6),(3,7),(4,5),(4,7),(5,6),(6,7)])
    #3
    #G = nx.from_edgelist([(0,1),(0,2),(0,3),(0,4),(0,6),(1,3),(1,4),(1,2),(1,5),(2,3),(2,5),(2,6),(3,4),(3,5),(3,6)])
    #4
    #G =nx.from_edgelist([(0,1),(0,7),(0,4),(1,9),(1,2),(2,3),(2,19),(3,4),(3,17),(4,5),(5,6),(5,16),(6,7),(6,13),(7,8),(8,9),(8,12),(9,10),(10,11),(10,19),(11,12),(11,15),(12,13),(13,14),(14,15),(14,16),(15,18),(16,17),(17,18),(18,19)])
    #5
    #G = nx.from_edgelist([(0,1),(0,2),(0,3),(0,4),(0,5),(0,6),(1,2),(1,3),(2,3),(2,4),(2,5),(2,6),(3,4),(4,5),(5,6)])
    #6
    #G=nx.from_edgelist([(0,1),(0,10),(0,9),(1,2),(1,11),(2,3),(2,11),(3,4),(3,12),(4,5),(4,12),(5,6),(5,13),(6,7),(6,13),(7,8),(7,14),(8,9),(8,14),(9,10),(10,15),(11,17),(12,18),(13,19),(14,16),(15,16),(15,17),(16,19),(17,18),(18,19)])
    upperBound = 30
    lowerBound = 0
    C = 1
    vars = Variables(G.nodes(), lowerBound, upperBound)
    cqm = build_cqm(vars, G, upperBound, C)
    best_feasible = call_solver(cqm)
    print(best_feasible)
    #best_feasible = {'x_0': 0, 'x_1': 6.0, 'x_2': 12.0, 'x_3': 6.0, 'x_4': 4.0, 'x_5': 8.0, 'x_6': 6.0, 'y_0': 0, 'y_1': 10, 'y_2': 0, 'y_3': 3.0, 'y_4': 4.0, 'y_5': 4.0, 'y_6': 1.0}
    pos = create_pos_for_drawing(best_feasible[0])
    
    nx.draw(G, pos=pos, node_size=300, edgecolors='k', cmap='hsv', with_labels=True)
    #print(G.edges())s
    #plt.savefig("titto_draw_second_spring_1C.png")
    #plt.clf()
    pos = nx.spring_layout(G)

    nx.draw(G, pos=pos, node_size=300, edgecolors='k', cmap='hsv', with_labels=True)
    #print(G.edges())s
    #plt.savefig("nx_draw_second_spring.png")
