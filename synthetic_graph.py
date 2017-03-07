#! coding = utf-8
#__author__ = "Xiaoyan"
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt 

import numpy as np 
import pandas as pd
np.random.seed(42) # so we can reproduce results
import random
random.seed(31) # so we can reproduce results

import networkx as nx

num_nodes = 500
num_cascades = 2000 
#g_weight = 1
g_weight = 1 
outputfile = "../data/cascade-synthentic-graph.txt" 
T = 0.7 # the observation window

# output the nodes as part of the input for the program
print "Write to", outputfile
txt = open(outputfile, "w+") 
for node in range(num_nodes):
  txt.write("%d,%d\n" % (node, node))
txt.write("\n")


def loadSBM(N, C, beta = 0.001):
    '''
    @N : number of nodes
    @C : size of communities 
    '''
    print "beta=", beta
    G = nx.Graph(gnc = {}, membership = {})
    for i in range(N):
      comm = int(i / C)
      G.graph["membership"][i] = [comm] 
      if comm in G.graph["gnc"]:
        G.graph["gnc"][ comm ].append( i ) 
      else:
        G.graph["gnc"][ comm ] = [ i ]

    for i in range(N):
      G.add_node(i)

    for comm in range(int((N + C - 1)/C)):
      nodes = G.graph["gnc"][ comm ]
      print nodes
      for l in range(len(nodes)):
        for r in range(l+1, len(nodes)):
          if np.random.rand() < 0.3:
            G.add_edge(nodes[l], nodes[r])

    for l in range(N-1):
      for r in range(l+1,N):
        if not G.has_edge(l,r):
          if np.random.rand() < beta:
            G.add_edge(l,r)

    #G.graph["membership"][i] = [comm] 
    return G

sG = loadSBM(num_nodes, 40)
#sG = nx.fast_gnp_random_graph(num_nodes, p=0.007, directed=True)
#sG = nx.barabasi_albert_graph(num_nodes, m=3)
#sG = nx.complete_graph(num_nodes)

print nx.info(sG)
G = nx.Graph()
casc_list = []

# generate weights
C = np.zeros( (num_nodes, num_nodes) )
for e in sG.edges():
  #G.add_edge(e[0],e[1], weight = 1.5 + random.random())
  G.add_edge(e[0],e[1], weight = g_weight / 2.0 + random.random() * g_weight / 2.0) 
  #G.add_edge(e[0],e[1], weight = 5)
  #G.add_edge(e[0],e[1], weight = 5)
  C[e[0],e[1]] = G[e[0]][e[1]]["weight"] 
  C[e[1],e[0]] = G[e[0]][e[1]]["weight"] 

C.dump("../data/C.dat")
#nx.write_gpickle(G, "../data/synthetic_graph.gpickle")

# plot data
#C = np.dot(A, B.transpose())
#M = np.squeeze(np.asarray(C))
#plt.hist(M, bins=10)
#plt.savefig("a.pdf")
#exit(1)

sizes = []
for casc in range(num_cascades):
  print casc
  infected_node = random.randint(0, num_nodes-1)
  current_time = 0
  earliest_time = pd.DataFrame({"init": [float("inf") for _ in range(num_nodes)]})
  Infected = {}
  Visited = set()
  output = [(infected_node, 0.0)]
  #print earliest_time 

  for loop in range(num_nodes-1):
    u = infected_node 
    Infected.pop(u, None)
    Visited.add(u)
    #print "update from", u, "at", current_time
    #print Infected.keys()

    # infection from node u ------> every node v
    for v in G.neighbors(u): 
      if v in Visited: continue
      if v in Infected:
        Infected[v] = min(Infected[v],\
                current_time + np.random.exponential( 1.0 / G[u][v]["weight"] ) )
      else:
        Infected[v] = current_time + np.random.exponential( 1.0 / G[u][v]["weight"] )

    # the earliest infection time 
    infected_node = min(Infected, key=Infected.get)
    # jump to the next infection time 
    current_time = Infected[infected_node]

    #print "infect", infected_node, "at time", current_time 
    output.append((infected_node, current_time))

    if (current_time > T):
      sizes.append(len(output))
      print "size", len(output)
      casc_list.append(output)
      break

  txt.write("%d;" % casc)
  txt.write(",".join(["%d,%.5f" % (node, time) for node, time in output]) + "\n") 


import pickle 
with open('../data/casc_list.pickle', 'w') as outfile:
    pickle.dump(casc_list, outfile)
print "write to",'../data/casc_list.pickle' 

print "Cascade Size"
for percent in [50,75,90,95]: 
  print percent, "\% :", np.percentile(sizes, percent) # return 50th percentile, e.g median.

print nx.info(sG)
txt.close()
