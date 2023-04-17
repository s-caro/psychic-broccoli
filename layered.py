import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
from dimod import ConstrainedQuadraticModel, Binary, quicksum, Real, SampleSet, Integer
from dimod import ConstrainedQuadraticModel, Binary, quicksum, cqm_to_bqm, BinaryQuadraticModel, ExactCQMSolver
from dwave.system import LeapHybridCQMSampler, DWaveSampler, EmbeddingComposite
import dimod
import json

import dwave.inspector
import pprint
from datetime import datetime




class Variables:
  def __init__(self, L1: dict, L2: dict):
      self.L1 = {(i,j): Binary(f'x_{i,j}') for (i,j) in L1.keys()}
      self.L2 = {(k,n): Binary(f'x_{k,n}') for (k,n) in L2.keys()}


def _add_variable(cqm: ConstrainedQuadraticModel, g: nx.Graph()):
  for index, (i,k) in enumerate(list(g.edges)[:-1]):
      for (j,n) in list(g.edges)[index+1:]:
        if i != j and k != n:
          cqm.add_variable('BINARY', f'x_{i,j}')
          cqm.add_variable('BINARY', f'x_{j,i}')
          cqm.add_variable('BINARY', f'x_{k,n}')
          cqm.add_variable('BINARY', f'x_{n,k}')


def _add_constraint_discrete(cqm: ConstrainedQuadraticModel, vars: Variables, G: nx.Graph):
  inserted = set()
  for index, (i,k) in enumerate(list(G.edges)[:-1]):
    for (j,n) in list(G.edges)[index+1:]:
      if i != j and k != n and (i,j):
        if (i,j) not in inserted:
          cqm.add_constraint(vars.L1[i,j] + vars.L1[j,i] == 1, label=f'L1: discrete on nodes {i,j}')
          inserted.add((i,j))
          inserted.add((j,i))
        if (k,n) not in inserted:
          cqm.add_constraint(vars.L2[k,n] + vars.L2[n,k] == 1, label=f'L2: discrete on nodes {k,n}')
          inserted.add((k,n))
          inserted.add((n,k))

def _add_constraint_transitive(cqm: ConstrainedQuadraticModel, vars: Variables, G: nx.Graph):
  for index, (i,j) in enumerate(list(vars.L1.keys())[:-1]):
    for (k,n) in list(vars.L1.keys())[index+1:]:
      if k==j and n!=i:
        cqm.add_constraint(vars.L1[i,j]+vars.L1[k,n] - vars.L1[i,n] <= 1, label=f'L1: <= transitive constraint nodes {i,k,n}')
        cqm.add_constraint(vars.L1[i,j]+vars.L1[k,n] - vars.L1[i,n] >= 0, label=f'L1: >= transitive constraint nodes {i,k,n}')
        cqm.add_constraint(vars.L1[n,k]+vars.L1[j,i] - vars.L1[n,i] <= 1, label=f'L1: <= transitive constraint nodes {n,k,i}')
        cqm.add_constraint(vars.L1[n,k]+vars.L1[j,i] - vars.L1[n,i] >= 0, label=f'L1: >= transitive constraint nodes {n,k,i}')
  for index, (i,j) in enumerate(list(vars.L2.keys())[:-1]):
    for (k,n) in list(vars.L2.keys())[index+1:]:
      if k==j and n!=i:
        cqm.add_constraint(vars.L2[i,j]+vars.L2[k,n] - vars.L2[i,n] <= 1, label=f'L2: <= transitive constraint nodes {i,k,n}')
        cqm.add_constraint(vars.L2[i,j]+vars.L2[k,n] - vars.L2[i,n] >= 0, label=f'L2: >= transitive constraint nodes {i,k,n}')
        cqm.add_constraint(vars.L2[n,k]+vars.L2[j,i] - vars.L2[n,i] <= 1, label=f'L1: <= transitive constraint nodes {n,k,i}') 
        cqm.add_constraint(vars.L2[n,k]+vars.L2[j,i] - vars.L2[n,i] >= 0, label=f'L1: >= transitive constraint nodes {n,k,i}')      

def _set_objective(cqm: ConstrainedQuadraticModel, vars: Variables, G: nx.Graph):
  tot = 0
  for index,a1 in enumerate(list(G.edges)[:-1]):
    for a2 in list(G.edges)[index+1:]:
      if a1[0] != a2[0] and a1[1] != a2[1]:
        tot += vars.L1[a1[0],a2[0]]*vars.L2[a2[1],a1[1]] + vars.L1[a2[0],a1[0]]*vars.L2[a1[1],a2[1]]
  cqm.set_objective(tot)


def build_cqm(G: nx.DiGraph) -> ConstrainedQuadraticModel:
  L1 = {}
  L2 = {}
  for index, (i,k) in enumerate(list(G.edges)[:-1]):
    for (j,n) in list(G.edges)[index+1:]:
      if i != j and k != n:
        L1[(i,j)]=None
        L1[(j,i)]=None
        L2[(k,n)]=None
        L2[(n,k)]=None
      elif i!=j:
        L1[(i,j)]=None
        L1[(j,i)]=None
      elif k!=n:
        L2[(k,n)]=None
        L2[(n,k)]=None
  vars = Variables(L1, L2)
  cqm = ConstrainedQuadraticModel()
  _add_variable(cqm, G)
  _add_constraint_discrete(cqm, vars, G)
  _add_constraint_transitive(cqm, vars, G)
  _set_objective(cqm, vars, G)

  return cqm


def ind_couples_edges(G: nx.Graph) -> int:
  indipendent = 0
  for index, (i,k) in enumerate(list(G.edges)[:-1]):
      for (j,n) in list(G.edges)[index+1:]:
        if i != j and k != n:
          indipendent+=1
  return indipendent

if __name__ == '__main__':
  for k in range(20,91,10):
    p = k/100
    for it in range(0,10):
      B = nx.bipartite.random_graph(10,10,p,seed=None, directed=False)
      while not (nx.is_connected(B) and len(B.edges)==k):
        B = nx.bipartite.random_graph(10,10,p,seed=None, directed=False)
      edges = sorted(list(B.edges))
      G = nx.DiGraph(edges)
      indipendent = ind_couples_edges(G)
      print("build model")
      cqm = build_cqm(G)
      print("model build")
      #print(cqm)
      sampler = LeapHybridCQMSampler()
      #4*sampler.min_time_limit(cqm),
      print("start sampling")
      res = sampler.sample_cqm(cqm, label="Titto - layer - hybrid - all runs - indipendent edges")
      feasible_sampleset = res.filter(lambda row: row.is_feasible)
      print("finish sampling")

      print("start exact solver")
      exactSampler = ExactCQMSolver()
      exactRes = exactSampler.sample_cqm(cqm, label = "Titto - layer - exact - all runs - independet edges")
      exact_feasible_sampleset = exactRes.filter(lambda row: row.is_feasible)
      print("finish exact solver")

      finalRes = {}
      finalRes['probability'] = p
      finalRes['number of edges'] = len(G.edges)
      finalRes['indipendent'] = ind_couples_edges(G)
      finalRes['edges'] = list(G.edges)
      finalRes.update(res.info)
      finalRes['sample']= feasible_sampleset.first
      finalRes['exact result']=exact_feasible_sampleset.first
      json_string = json.dumps(finalRes)
      print("dumped result")

      with open('layered_cqm_10_samples_all_info.txt', 'a') as f:
        f.write('\n' +datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        f.write('\n\n-------------------------\n\n')
        f.write(json_string)
        f.write('\n\n-------------------------\n\n')

    #print(feasible_sampleset.first)
    #print("start conversion")
    #bqm, inverse = dimod.cqm_to_bqm(cqm)
    #print("finish conversion")
    #pprint.pprint(bqm)
    #sampler = EmbeddingComposite(DWaveSampler())
    
    #print("start sampling")
    #sampleset = sampler.sample(bqm, num_reads=1000, label='Titto - layer - bqm')
    #pprint.pprint(sampleset.first)
    #print("finish sampling")
    #with open('layered_bqm_ex.txt', 'a') as f:
    #  f.write('\n' +datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    #  f.write('\n\n-------------------------\n\n')
    #  f.write(f'probability: {p}\nnumber of edges: {len(G.edges)}')
    #  f.write('\n\n-------------------------\n\n')
    #  for i in G.edges:
    #      f.write(', ' + str(i))
    #  #f.write('\n\n-------------------------\n\n')
    #  #f.write(str(bqm))
    #  f.write('\n\n-------------------------\n\n')
    #  for j in sampleset.info:
    #    f.write(str(j) + ":\t" + str(sampleset.info[j]) + "\n")
    #  f.write('\n\n-------------------------\n\n')
    #  f.write(str(sampleset.first))
    #  f.write('\n\n-------------------------\n\n')

    print(f"end it {k}")