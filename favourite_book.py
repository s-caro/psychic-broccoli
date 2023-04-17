import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
from dimod import ConstrainedQuadraticModel, Binary, quicksum, Real, SampleSet, Integer
from dimod import ConstrainedQuadraticModel, Binary, quicksum, cqm_to_bqm, BinaryQuadraticModel
from dwave.system import LeapHybridCQMSampler, DWaveSampler, EmbeddingComposite
import dimod
import numpy as np
import pprint
import pandas as pd
from datetime import datetime


class Variables:
    def __init__(self, L1: dict):
        self.L1 = {(i,j):('BINARY', f'x_{i,j}') for (i,j) in L1.keys()}


class BQM_Variables:
    def __init__(self, L1: dict):
        self.L1 = {(i,j):Binary(f'x_{i,j}') for (i,j) in L1.keys()}


def _add_variable(cqm: ConstrainedQuadraticModel, g: nx.Graph()):
    for index, (i,k) in enumerate(list(G.edges)[:-1]):
        for (j,n) in list(G.edges)[index+1:]:
            if i != j and k != n and i != n and k != j:
                cqm.add_variable('BINARY', f'x_{i,j}')
                cqm.add_variable('BINARY', f'x_{j,i}')
                cqm.add_variable('BINARY', f'x_{k,n}')
                cqm.add_variable('BINARY', f'x_{n,k}')
                cqm.add_variable('BINARY', f'x_{i,n}')
                cqm.add_variable('BINARY', f'x_{n,i}')
                cqm.add_variable('BINARY', f'x_{k,j}')
                cqm.add_variable('BINARY', f'x_{j,k}')
                cqm.add_variable('BINARY', f'x_{i,k}')
                cqm.add_variable('BINARY', f'x_{k,i}')
                cqm.add_variable('BINARY', f'x_{n,j}')
                cqm.add_variable('BINARY', f'x_{j,n}')


def _add_constraint_discrete(cqm: ConstrainedQuadraticModel, vars: Variables, G: nx.Graph):
    inserted = set()
    for index, (i,k) in enumerate(list(G.edges)[:-1]):
        for (j,n) in list(G.edges)[index+1:]:
            if i != j and k != n and i != n and k != j:
                if (i,j) not in inserted:
                    cqm.add_constraint(vars.L1[i,j] + vars.L1[j,i] == 1, label=f'L1: discrete on nodes {i,j}')
                    inserted.add((i,j))
                    inserted.add((j,i))
                if (k,n) not in inserted: 
                    cqm.add_constraint(vars.L1[k,n] + vars.L1[n,k] == 1, label=f'L1: discrete on nodes {k,n}')
                    inserted.add((k,n))
                    inserted.add((n,k))
                if (i,n) not in inserted:
                    cqm.add_constraint(vars.L1[i,n] + vars.L1[n,i] == 1, label=f'L1: discrete on nodes {i,n}')
                    inserted.add((i,n))
                    inserted.add((n,i))
                if (k,j) not in inserted:
                    cqm.add_constraint(vars.L1[k,j] + vars.L1[j,k] == 1, label=f'L1: discrete on nodes {k,j}')
                    inserted.add((k,j))
                    inserted.add((j,k))
                if (i,k) not in inserted:
                    cqm.add_constraint(vars.L1[k,i] + vars.L1[i,k] == 1, label=f'L1: discrete on nodes {k,i}')
                    inserted.add((k,i))
                    inserted.add((i,k))
                if (n,j) not in inserted:
                    cqm.add_constraint(vars.L1[j,n] + vars.L1[n,j] == 1, label=f'L1: discrete on nodes {j,n}')
                    inserted.add((j,n))
                    inserted.add((n,j))    


def _add_constraint_transitive(cqm: ConstrainedQuadraticModel, vars: Variables, G: nx.Graph):
    for index, (i,j) in enumerate(list(vars.L1.keys())[:-1]):
        for (k,n) in list(vars.L1.keys())[index+1:]:
            if k==j and n!=i:
                cqm.add_constraint(vars.L1[i,j]+vars.L1[k,n] - vars.L1[i,n] <= 1, label=f'L1: <= transitive constraint nodes {i,k,n}')
                cqm.add_constraint(vars.L1[i,j]+vars.L1[k,n] - vars.L1[i,n] >= 0, label=f'L1: >= transitive constraint nodes {i,k,n}')
                cqm.add_constraint(vars.L1[n,k]+vars.L1[j,i] - vars.L1[n,i] <= 1, label=f'L1: <= transitive constraint nodes {n,k,i}')
                cqm.add_constraint(vars.L1[n,k]+vars.L1[j,i] - vars.L1[n,i] >= 0, label=f'L1: >= transitive constraint nodes {n,k,i}')   


def _set_objective(cqm: ConstrainedQuadraticModel, vars: BQM_Variables, G: nx.Graph) -> ConstrainedQuadraticModel:
    bqm = BinaryQuadraticModel(dimod.BINARY)
    for index,a1 in enumerate(list(G.edges)[:-1]):
        for a2 in list(G.edges)[index+1:]:
            if a1[0] != a2[0] and a1[1] != a2[1] and a1[0] != a2[1] and a1[1] != a2[0]:
                a = a1[0]
                c = a1[1]
                b = a2[0]
                d = a2[1]
                poly = {(vars.L1[a,b].variables[0],vars.L1[b,c].variables[0],vars.L1[c,d].variables[0]):1,(vars.L1[a,d].variables[0], vars.L1[d,c].variables[0], vars.L1[c,b].variables[0]):1, (vars.L1[c,b].variables[0],vars.L1[b,a].variables[0],vars.L1[a,d].variables[0]):1, (vars.L1[c,d].variables[0],vars.L1[d,a].variables[0],vars.L1[a,b].variables[0]):1,(vars.L1[b,a].variables[0],vars.L1[a,d].variables[0],vars.L1[d,c].variables[0]):1,(vars.L1[b,c].variables[0],vars.L1[c,d].variables[0],vars.L1[d,a].variables[0]):1,(vars.L1[d,a].variables[0],vars.L1[a,b].variables[0],vars.L1[b,c].variables[0]):1,(vars.L1[d,c].variables[0],vars.L1[c,b].variables[0],vars.L1[b,a].variables[0]):1}
                cqm = dimod.make_quadratic_cqm(poly, dimod.BINARY)
                bqm += dimod.make_quadratic(poly, 1.0, dimod.BINARY)
    return cqm.from_bqm(bqm)



def build_cqm(G: nx.Graph) -> ConstrainedQuadraticModel:
    L1 = {}
    for index, (i,k) in enumerate(list(G.edges)[:-1]):
        for (j,n) in list(G.edges)[index+1:]:
            if i != j and k != n and i != n and k != j:
                L1[(i,j)]=None
                L1[(j,i)]=None
                L1[(i,n)]=None
                L1[(n,i)]=None
                L1[(k,j)]=None
                L1[(j,k)]=None
                L1[(k,n)]=None
                L1[(n,k)]=None
                L1[(i,k)]=None
                L1[(k,i)]=None
                L1[(j,n)]=None
                L1[(n,j)]=None

    BQM_variables = BQM_Variables(L1)
    cqm = ConstrainedQuadraticModel()
    _add_variable(cqm, G)
    cqm = _set_objective(cqm, BQM_variables, G)
    _add_constraint_discrete(cqm, BQM_variables, G)
    _add_constraint_transitive(cqm, BQM_variables, G)

    return cqm

if __name__ == '__main__':
    edges = sorted([(1,2),(3,4),(5,6)])
    G = nx.DiGraph(edges)
    cqm = build_cqm(G)
    sampler = LeapHybridCQMSampler()
    #4*sampler.min_time_limit(cqm),
    res = sampler.sample_cqm(cqm,  label="Titto - book emb - hybrid")
    feasible_sampleset = res.filter(lambda row: row.is_feasible)
    pprint.pprint(feasible_sampleset.first)
    

    with open('favourite_book_cqm.txt', 'a') as f:
        f.write('\n' +datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        f.write('\n\n-------------------------\n\n')
        for i in G.edges:
            f.write(', ' + str(i))
        #f.write('\n\n-------------------------\n\n')
        #f.write(str(cqm))
        f.write('\n\n-------------------------\n\n')
        f.write(str(res.first))
        f.write('\n\n-------------------------\n\n')


    bqm, inverse = dimod.cqm_to_bqm(cqm)
    #pprint.pprint(bqm)
    sampler = EmbeddingComposite(DWaveSampler())
    

    sampleset = sampler.sample(bqm, num_reads=2000, label='Titto - book emb - bqm')
    pprint.pprint(sampleset.first)

    with open('favourite_book_bqm.txt', 'a') as f:
        f.write('\n' +datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        f.write('\n\n-------------------------\n\n')
        for i in G.edges:
            f.write(', ' + str(i))
        #f.write('\n\n-------------------------\n\n')
        #f.write(str(bqm))
        f.write('\n\n-------------------------\n\n')
        f.write(str(sampleset.first))
        f.write('\n\n-------------------------\n\n')