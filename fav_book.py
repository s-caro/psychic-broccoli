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
import time

class Variables:
  def __init__(self, L1: dict):
      self.L1 = {(i,j):Binary(f'x_{i,j}') for (i,j) in L1.keys()}


class Color_Variables:
  def __init__(self, edges: list, colors: list):
    self.edgeColor = {(i,j,k):Binary(f'x_{i,j,k}') for k in colors for (i,j) in edges}

def _add_variable(cqm: ConstrainedQuadraticModel, g: nx.Graph()):
  for index, a in enumerate(list(G.nodes)[:-1]):
    for b in list(G.nodes)[index+1:]:
      cqm.add_variable('BINARY', f'x_{a,b}')
      cqm.add_variable('BINARY', f'x_{b,a}')

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



def _add_color_cross_constraint(cqm: ConstrainedQuadraticModel, vars: Variables, colVar: Color_Variables, colors: list, G: nx.Graph):
  for index,a1 in enumerate(list(G.edges)[:-1]):
    for a2 in list(G.edges)[index+1:]:
      if a1[0] != a2[0] and a1[1] != a2[1] and a1[0] != a2[1] and a1[1] != a2[0]:
          for k in colors:
            a = a1[0]
            c = a1[1]
            b = a2[0]
            d = a2[1]

            cqm.add_constraint(colVar.edgeColor[(a,c,k)] + colVar.edgeColor[(b,d,k)] + vars.L1[a,b] + vars.L1[b,c] + vars.L1[c,d]<=4, label = f'color {k} cross constraint edges {a,b,c,d}')
            cqm.add_constraint(colVar.edgeColor[(a,c,k)] + colVar.edgeColor[(b,d,k)] + vars.L1[b,a] + vars.L1[a,d] + vars.L1[d,c]<=4, label = f'color {k} cross constraint edges {b,a,d,c}')
            cqm.add_constraint(colVar.edgeColor[(a,c,k)] + colVar.edgeColor[(b,d,k)] + vars.L1[a,d] + vars.L1[d,c] + vars.L1[c,b]<=4, label = f'color {k} cross constraint edges {a,d,c,b}')
            cqm.add_constraint(colVar.edgeColor[(a,c,k)] + colVar.edgeColor[(b,d,k)] + vars.L1[b,c] + vars.L1[c,d] + vars.L1[d,a]<=4, label = f'color {k} cross constraint edges {b,c,d,a}')
            cqm.add_constraint(colVar.edgeColor[(a,c,k)] + colVar.edgeColor[(b,d,k)] + vars.L1[c,b] + vars.L1[b,a] + vars.L1[a,d]<=4, label = f'color {k} cross constraint edges {c,b,a,d}')
            cqm.add_constraint(colVar.edgeColor[(a,c,k)] + colVar.edgeColor[(b,d,k)] + vars.L1[d,a] + vars.L1[a,b] + vars.L1[b,c]<=4, label = f'color {k} cross constraint edges {d,a,b,c}')
            cqm.add_constraint(colVar.edgeColor[(a,c,k)] + colVar.edgeColor[(b,d,k)] + vars.L1[c,d] + vars.L1[d,a] + vars.L1[a,b]<=4, label = f'color {k} cross constraint edges {c,d,a,b}')
            cqm.add_constraint(colVar.edgeColor[(a,c,k)] + colVar.edgeColor[(b,d,k)] + vars.L1[d,c] + vars.L1[c,b] + vars.L1[b,a]<=4, label = f'color {k} cross constraint edges {d,c,b,a}')

def _add_constraint_discrete(cqm: ConstrainedQuadraticModel, vars: Variables, G: nx.Graph):
  for index, a in enumerate(list(G.nodes)[:-1]):
    for b in list(G.nodes)[index+1:]:
      cqm.add_constraint(vars.L1[a,b] + vars.L1[b,a] == 1, label=f'discrete on nodes {a,b}')

def _add_constraint_transitive(cqm: ConstrainedQuadraticModel, vars: Variables, G: nx.Graph):
  for index, (i,j) in enumerate(list(vars.L1.keys())[:-1]):
    for (k,n) in list(vars.L1.keys())[index+1:]:
      if k==j and n!=i:
        cqm.add_constraint(vars.L1[i,j]+vars.L1[k,n] - vars.L1[i,n] <= 1, label=f'<= transitive constraint nodes {i,k,n}')
        cqm.add_constraint(vars.L1[i,j]+vars.L1[k,n] - vars.L1[i,n] >= 0, label=f'>= transitive constraint nodes {i,k,n}')
        cqm.add_constraint(vars.L1[n,k]+vars.L1[j,i] - vars.L1[n,i] <= 1, label=f'<= transitive constraint nodes {n,k,i}')
        cqm.add_constraint(vars.L1[n,k]+vars.L1[j,i] - vars.L1[n,i] >= 0, label=f'>= transitive constraint nodes {n,k,i}')   

def _add_transitive_implication_constraint(cqm: ConstrainedQuadraticModel, vars: Variables, G: nx.Graph):
  for index, (i,j) in enumerate(list(vars.L1.keys())[:-1]):
    for (k,n) in list(vars.L1.keys())[index+1:]:
      if k==j and n!=i:
        cqm.add_constraint(1-(vars.L1[i,j]*vars.L1[k,n])+vars.L1[i,n] >= 1, label=f'L1: transitive constraint nodes {i,k,n}')
        cqm.add_constraint(1-(vars.L1[n,k]*vars.L1[j,i])+vars.L1[n,i] >= 1, label=f'L1: transitive constraint nodes {n,k,i}')


def build_cqm(G: nx.Graph) -> ConstrainedQuadraticModel:
  colors = ['b','y','r','p']
  L1 = {}
  L3 = []
  for (i,j) in G.edges:
    L3.append((i,j))
    L3.append((j,i))
  for index, a in enumerate(list(G.nodes)[:-1]):
    for b in list(G.nodes)[index+1:]:
      L1[(a,b)]=None
      L1[(b,a)]=None
  variables = Variables(L1)
  color_variables = Color_Variables(L3, colors)
  cqm = ConstrainedQuadraticModel()
  _add_variable(cqm, G)
  _add_Color_variable(cqm, G.edges, colors)
  _add_constraint_discrete(cqm, variables, G)
  _add_transitive_implication_constraint(cqm, variables, G)
  _add_same_edge_color_constraint(cqm, G.edges, color_variables, colors)
  _add_color_discrete_constraint(cqm, L3, color_variables, colors)
  _add_color_cross_constraint(cqm, variables, color_variables, colors, G)

  return cqm

def convert_file_to_list(line): 
  edge_list = []
  for el in line:
    if len(el) == 4:
      if int(el[0]) != int(el[3]):
        edge = sorted((int(el[0]),int(el[3])))
        edge_list.append(edge)
    if len(el) == 6:
      if int(el[0]) != int(el[4]) or int(el[1]) != int(el[5]):
        num1 = int(el[0])*10 + int(el[1])
        num2 = int(el[4])*10 + int(el[5])
        edge = sorted((num1, num2))
        edge_list.append(edge)
    if len(el) == 5:
      if el[1] != ',':
        num1 = int(el[0])*10 + int(el[1])
        edge = sorted((num1, int(el[4])))
        edge_list.append(edge)
      else:
        num2 = int(el[3])*10 + int(el[4])
        edge = sorted((int(el[0]), num2))
        edge_list.append(edge)
  return edge_list

if __name__ == '__main__':
  file_edges = open('evaluation/planarConnectedGraph_nodes_5_30_edges_min_max.txt','r')
  lines = file_edges.readlines()

  f = open('4-page-book-embedding_time_noObj.csv','a')
  header =['# nodes', '# edges', '#is feasible', '# crosses', 'qpu_access_time', 'charge_time','run_time', 'total time','num_Constraints']
  writer = csv.writer(f)
  writer.writerow(header)
  f.close()
  
  for index,line in enumerate(lines[0:]):
    f = open('4-page-book-embedding_time_noObj.csv','a')
    clear_line = line.replace('\n','').replace(')','').replace('(','').split('_')
    print("start reading")
    edges = sorted(convert_file_to_list(clear_line))

    #edges = [(3,5),(3,6),(1,4),(1,5),(2,5)]
    G = nx.from_edgelist(edges)
    print("start building")
    cqm = build_cqm(G)
    sampler = LeapHybridCQMSampler()
    #4*sampler.min_time_limit(cqm),
    print("start sampling")
    start_time = time.time()
   
    res = sampler.sample_cqm(cqm,  label="Titto - book emb - hybrid")
    end = time.time() - start_time
    try:
      feasible_sampleset = res.filter(lambda row: row.is_feasible)
      #pprint.pprint(feasible_sampleset.first)
      data = [len(G.nodes), len(G.edges), 'yes', feasible_sampleset.first.energy, res.info['qpu_access_time'], res.info['charge_time'], res.info['run_time'], end,len(cqm.constraints)]
    except (StopIteration, ValueError) as error:
      data = [len(G.nodes), len(G.edges), 'no', '//', res.info['qpu_access_time'], res.info['charge_time'], res.info['run_time'], end,len(cqm.constraints)]
    writer = csv.writer(f)
    writer.writerow(data)
    f.close()
    print(f"done: {index} / {len(lines)}")

