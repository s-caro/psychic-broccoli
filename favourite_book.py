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
import csv


class BQM_Variables:
  def __init__(self, L1: dict):
      self.L1 = {(i,j):Binary(f'x_{i,j}') for (i,j) in L1.keys()}

class Z_Variables:
  def __init__(self, ord: dict):
    self.ord = {(a,b,c,d):Binary(f'x_{a,b,c,d}') for (a,b,c,d) in ord.keys()}

class Color_Variables:
  def __init__(self, edges: list, colors: list):
    self.edgeColor = {(i,j,k):Binary(f'x_{i,j,k}') for k in colors for (i,j) in edges}

def _add_variable(cqm: ConstrainedQuadraticModel, g: nx.Graph()):
  for index, a in enumerate(list(G.nodes)[:-1]):
    for b in list(G.nodes)[index+1:]:
      cqm.add_variable('BINARY', f'x_{a,b}')
      cqm.add_variable('BINARY', f'x_{b,a}')

def _add_Z_Variable(cqm: ConstrainedQuadraticModel, ord: dict):
  for (a,b,c,d) in ord.keys():
      cqm.add_variable('BINARY', f'x_{a,b,c,d}')

def _add_Color_variable(cqm:ConstrainedQuadraticModel, edges: list, colors: list):
  for (i,j) in edges:
    for k in colors:
      cqm.add_variable('BINARY', f'x_{i,j,k}')

def _add_color_discrete_constraint(cqm: ConstrainedQuadraticModel, edges: list, colVar: Color_Variables, colors: list):
  for (i,j) in edges:
    cqm.add_constraint(quicksum(colVar.edgeColor[(i,j,k)] for k in colors) == 1, label = f'discrete constraint on edge {i,j}')

def _add_same_edge_color_constraint(cqm: ConstrainedQuadraticModel, edges: nx.Graph.edges, colVar: Color_Variables, colors: list):
  for (i,j) in edges:
    for k in colors:
      cqm.add_constraint(colVar.edgeColor[(i,j,k)] - colVar.edgeColor[(j,i,k)] == 0, label = f'discrete constraint same edge same color edge {i,j} color {k}')

def _add_color_cross_constraint(cqm: ConstrainedQuadraticModel, edges: list, colVar: Color_Variables, zVar: Z_Variables, colors: list):
  for (a,b,c,d) in zVar.ord:
    for k in colors:
      cqm.add_constraint(zVar.ord[(a,b,c,d)]*(colVar.edgeColor[(a,c,k)] + colVar.edgeColor[(b,d,k)]) <= 1, label=f'color cross constraints on cross {a,b,c,d} on color {k}')


def _add_constraint_discrete(cqm: ConstrainedQuadraticModel, vars: BQM_Variables, G: nx.Graph):
  for index, a in enumerate(list(G.nodes)[:-1]):
    for b in list(G.nodes)[index+1:]:
      cqm.add_constraint(vars.L1[a,b] + vars.L1[b,a] == 1, label=f'discrete on nodes {a,b}')

def _add_constraint_transitive(cqm: ConstrainedQuadraticModel, vars: BQM_Variables, G: nx.Graph):
  for index, (i,j) in enumerate(list(vars.L1.keys())[:-1]):
    for (k,n) in list(vars.L1.keys())[index+1:]:
      if k==j and n!=i:
        cqm.add_constraint(vars.L1[i,j]+vars.L1[k,n] - vars.L1[i,n] <= 1, label=f'<= transitive constraint nodes {i,k,n}')
        cqm.add_constraint(vars.L1[i,j]+vars.L1[k,n] - vars.L1[i,n] >= 0, label=f'>= transitive constraint nodes {i,k,n}')
        cqm.add_constraint(vars.L1[n,k]+vars.L1[j,i] - vars.L1[n,i] <= 1, label=f'<= transitive constraint nodes {n,k,i}')
        cqm.add_constraint(vars.L1[n,k]+vars.L1[j,i] - vars.L1[n,i] >= 0, label=f'>= transitive constraint nodes {n,k,i}')   

def _add_z_constraint(cqm: ConstrainedQuadraticModel, vars: BQM_Variables, z_vars: Z_Variables, G: nx.Graph):
  for index,a1 in enumerate(list(G.edges)[:-1]):
    for a2 in list(G.edges)[index+1:]:
      if a1[0] != a2[0] and a1[1] != a2[1] and a1[0] != a2[1] and a1[1] != a2[0]:
        a = a1[0]
        c = a1[1]
        b = a2[0]
        d = a2[1]
        cqm.add_constraint(z_vars.ord[(a,b,c,d)] - vars.L1[a,b] <= 0, label=f'z constraints nodes {a,b} order {a,b,c,d}')
        cqm.add_constraint(z_vars.ord[(a,b,c,d)] - vars.L1[b,c] <= 0, label=f'z constraints nodes {b,c} order {a,b,c,d}')
        cqm.add_constraint(z_vars.ord[(a,b,c,d)] - vars.L1[c,d] <= 0, label=f'z constraints nodes {c,d} order {a,b,c,d}')
        cqm.add_constraint(z_vars.ord[(a,b,c,d)] - vars.L1[a,b] - vars.L1[b,c] - vars.L1[c,d] >= -2, label=f'z constraints all nodes order {a,b,c,d}')

        cqm.add_constraint(z_vars.ord[(b,a,d,c)] - vars.L1[b,a] <= 0, label=f'z constraints nodes {b,a} order {b,a,d,c}')
        cqm.add_constraint(z_vars.ord[(b,a,d,c)] - vars.L1[a,d] <= 0, label=f'z constraints nodes {a,d} order {b,a,d,c}')
        cqm.add_constraint(z_vars.ord[(b,a,d,c)] - vars.L1[d,c] <= 0, label=f'z constraints nodes {d,c} order {b,a,d,c}')
        cqm.add_constraint(z_vars.ord[(b,a,d,c)] - vars.L1[b,a] - vars.L1[a,d] - vars.L1[d,c] >= -2, label=f'z constraints all nodes order {b,a,d,c}')

        cqm.add_constraint(z_vars.ord[(a,d,c,b)] - vars.L1[a,d] <= 0, label=f'z constraints nodes {a,d} order {a,d,c,b}')
        cqm.add_constraint(z_vars.ord[(a,d,c,b)] - vars.L1[d,c] <= 0, label=f'z constraints nodes {d,c} order {a,d,c,b}')
        cqm.add_constraint(z_vars.ord[(a,d,c,b)] - vars.L1[c,b] <= 0, label=f'z constraints nodes {c,b} order {a,d,c,b}')
        cqm.add_constraint(z_vars.ord[(a,d,c,b)] - vars.L1[a,d] - vars.L1[d,c] - vars.L1[c,b] >= -2, label=f'z constraints all nodes order {a,d,c,b}')

        cqm.add_constraint(z_vars.ord[(b,c,d,a)] - vars.L1[b,c] <= 0, label=f'z constraints nodes {b,c} order {b,c,d,a}')
        cqm.add_constraint(z_vars.ord[(b,c,d,a)] - vars.L1[c,d] <= 0, label=f'z constraints nodes {c,d} order {b,c,d,a}')
        cqm.add_constraint(z_vars.ord[(b,c,d,a)] - vars.L1[d,a] <= 0, label=f'z constraints nodes {d,a} order {b,c,d,a}')
        cqm.add_constraint(z_vars.ord[(b,c,d,a)] - vars.L1[b,c] - vars.L1[c,d] - vars.L1[d,a] >= -2, label=f'z constraints all nodes order {b,c,d,a}')

        cqm.add_constraint(z_vars.ord[(c,b,a,d)] - vars.L1[c,b] <= 0, label=f'z constraints nodes {c,b} order {c,b,a,d}')
        cqm.add_constraint(z_vars.ord[(c,b,a,d)] - vars.L1[b,a] <= 0, label=f'z constraints nodes {b,a} order {c,b,a,d}')
        cqm.add_constraint(z_vars.ord[(c,b,a,d)] - vars.L1[a,d] <= 0, label=f'z constraints nodes {a,d} order {c,b,a,d}')
        cqm.add_constraint(z_vars.ord[(c,b,a,d)] - vars.L1[c,b] - vars.L1[b,a] - vars.L1[a,d] >= -2, label=f'z constraints all nodes order {c,b,a,d}')

        cqm.add_constraint(z_vars.ord[(d,a,b,c)] - vars.L1[d,a] <= 0, label=f'z constraints nodes {d,a} order {d,a,b,c}')
        cqm.add_constraint(z_vars.ord[(d,a,b,c)] - vars.L1[a,b] <= 0, label=f'z constraints nodes {a,b} order {d,a,b,c}')
        cqm.add_constraint(z_vars.ord[(d,a,b,c)] - vars.L1[b,c] <= 0, label=f'z constraints nodes {b,c} order {d,a,b,c}')
        cqm.add_constraint(z_vars.ord[(d,a,b,c)] - vars.L1[d,a] - vars.L1[a,b] - vars.L1[b,c] >= -2, label=f'z constraints all nodes order {d,a,b,c}')

        cqm.add_constraint(z_vars.ord[(c,d,a,b)] - vars.L1[c,d] <= 0, label=f'z constraints nodes {c,d} order {c,d,a,b}')
        cqm.add_constraint(z_vars.ord[(c,d,a,b)] - vars.L1[d,a] <= 0, label=f'z constraints nodes {d,a} order {c,d,a,b}')
        cqm.add_constraint(z_vars.ord[(c,d,a,b)] - vars.L1[a,b] <= 0, label=f'z constraints nodes {a,b} order {c,d,a,b}')
        cqm.add_constraint(z_vars.ord[(c,d,a,b)] - vars.L1[c,d] - vars.L1[d,a] - vars.L1[c,d] >= -2, label=f'z constraints all nodes order {c,d,a,b}')

        cqm.add_constraint(z_vars.ord[(d,c,b,a)] - vars.L1[d,c] <= 0, label=f'z constraints nodes {d,c} order {d,c,b,a}')
        cqm.add_constraint(z_vars.ord[(d,c,b,a)] - vars.L1[c,b] <= 0, label=f'z constraints nodes {c,b} order {d,c,b,a}')
        cqm.add_constraint(z_vars.ord[(d,c,b,a)] - vars.L1[b,a] <= 0, label=f'z constraints nodes {b,a} order {d,c,b,a}')
        cqm.add_constraint(z_vars.ord[(d,c,b,a)] - vars.L1[d,c] - vars.L1[c,b] - vars.L1[b,a] >= -2, label=f'z constraints all nodes order {d,c,b,a}')

def _set_objective(cqm: ConstrainedQuadraticModel, z_vars: Z_Variables, G: nx.Graph) -> ConstrainedQuadraticModel:
  tot = 0
  for index,a1 in enumerate(list(G.edges)[:-1]):
    for a2 in list(G.edges)[index+1:]:
      if a1[0] != a2[0] and a1[1] != a2[1] and a1[0] != a2[1] and a1[1] != a2[0]:
        a = a1[0]
        c = a1[1]
        b = a2[0]
        d = a2[1]
        tot += z_vars.ord[(a,b,c,d)] + z_vars.ord[(b,a,d,c)] + z_vars.ord[(a,d,c,b)] + z_vars.ord[(b,c,d,a)] + z_vars.ord[(c,b,a,d)] + z_vars.ord[(d,a,b,c)] + z_vars.ord[(c,d,a,b)] + z_vars.ord[(d,c,b,a)]
  cqm.set_objective(tot)

def build_cqm(G: nx.Graph) -> ConstrainedQuadraticModel:
  colors = ['b','y','r','p']
  L1 = {}
  L2 = {}
  L3 = []
  for (i,j) in G.edges:
    L3.append((i,j))
    L3.append((j,i))
  for index, a in enumerate(list(G.nodes)[:-1]):
    for b in list(G.nodes)[index+1:]:
      L1[(a,b)]=None
      L1[(b,a)]=None
  for index,a1 in enumerate(list(G.edges)[:-1]):
      for a2 in list(G.edges)[index+1:]:
        if a1[0] != a2[0] and a1[1] != a2[1] and a1[0] != a2[1] and a1[1] != a2[0]:
          a = a1[0]
          c = a1[1]
          b = a2[0]
          d = a2[1]
          L2[(a,b,c,d)]=None
          L2[(b,a,d,c)]=None
          L2[(a,d,c,b)]=None
          L2[(b,c,d,a)]=None
          L2[(c,b,a,d)]=None
          L2[(d,a,b,c)]=None
          L2[(c,d,a,b)]=None
          L2[(d,c,b,a)]=None
  BQM_variables = BQM_Variables(L1)
  z_var = Z_Variables(L2)
  color_variables = Color_Variables(L3, colors)
  cqm = ConstrainedQuadraticModel()
  _add_variable(cqm, G)
  _add_Z_Variable(cqm, L2)
  _add_Color_variable(cqm, G.edges, colors)
  _add_constraint_discrete(cqm, BQM_variables, G)
  _add_constraint_transitive(cqm, BQM_variables, G)
  _add_z_constraint(cqm, BQM_variables, z_var, G)
  _add_same_edge_color_constraint(cqm, G.edges, color_variables, colors)
  _add_color_discrete_constraint(cqm, L3, color_variables, colors)
  _add_color_cross_constraint(cqm, L3, color_variables, z_var, colors)
  _set_objective(cqm, z_var, G)

  return cqm

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
  file_edges = open('evaluation/planarConnectedGraph_nodes_5_30_edges_min_max.txt','r')
  lines = file_edges.readlines()

  f = open('4-page-book-embedding_time.csv','a')
  header =['# nodes', '# edges', '#is feasible', '# crosses', 'qpu_access_time', 'charge_time','run_time']
  writer = csv.writer(f)
  writer.writerow(header)
  f.close()
  
  for index,line in enumerate(lines[173:]):
    f = open('4-page-book-embedding_time.csv','a')
    clear_line = line.replace('\n','').replace(')','').replace('(','').split('_')
    edges = sorted(convert_file_to_list(clear_line))

    #edges = [(3,5),(3,6),(1,4),(1,5),(2,5)]
    G = nx.from_edgelist(edges)
    cqm = build_cqm(G)
    sampler = LeapHybridCQMSampler()
    #4*sampler.min_time_limit(cqm),
    res = sampler.sample_cqm(cqm,  label="Titto - book emb - hybrid")
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

