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

def convert_file_to_list(line): 
  edge_list = []
  for el in line:
    if len(el) == 4:
      if int(el[0]) != int(el[3]):
        edge_list.append((int(el[0]),int(el[3])))
    if len(el) == 6:
      if int(el[0]) != int(el[4]) or int(el[1]) != int(el[5]):
        num1 = int(el[0])*10 + int(el[1])
        num2 = int(el[4])*10 + int(el[5])
        edge_list.append((num1, num2))
    if len(el) == 5:
      if el[1] != ',':
        num1 = int(el[0])*10 + int(el[1])
        edge_list.append((num1, int(el[4])))
      else:
        num2 = int(el[3])*10 + int(el[4])
        edge_list.append((int(el[0]), num2))
  return edge_list

if __name__ == '__main__':
  file_edges = open('evaluation/dataset1_10_10_(5,95).txt','r')
  lines = file_edges.readlines()

  f = open('2-layer-cross-min.csv','a')
  header =['# nodes', '# edges', '#is feasible', '# crosses', 'qpu_access_time', 'charge_time','run_time', 'total time']
  writer = csv.writer(f)
  writer.writerow(header)
  f.close()
  
  for index,line in enumerate(lines[173:]):
    f = open('2-layer-cross-min.csv','a')
    clear_line = line.replace('\n','').replace(')','').replace('(','').split('_')
    edges = sorted(convert_file_to_list(clear_line))

    #edges = [(3,5),(3,6),(1,4),(1,5),(2,5)]
    G = nx.from_edgelist(edges)
    cqm = build_cqm(G)
    sampler = LeapHybridCQMSampler()
    #4*sampler.min_time_limit(cqm),
    res = sampler.sample_cqm(cqm,  label="Titto - cross min - hybrid")
    try:
      feasible_sampleset = res.filter(lambda row: row.is_feasible)
      #pprint.pprint(feasible_sampleset.first)
      data = [str(len(G.nodes)), str(len(G.edges)), 'yes', str(feasible_sampleset.first.energy), str(res.info['qpu_access_time']), str(res.info['charge_time']), str(res.info['run_time'])]
    except (StopIteration, ValueError) as error:
      data = [str(len(G.nodes)), str(len(G.edges)), 'no', '//', str(res.info['qpu_access_time']), str(res.info['charge_time']), str(res.info['run_time'])]
    writer = csv.writer(f)
    writer.writerow(data)
    f.close()
    print(f"done: {index+173} / {len(lines)}")

