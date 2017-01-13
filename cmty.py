#!/usr/bin/env python
import networkx as nx
import math
import csv
import random as rand
import sys
import time
import matplotlib.pyplot as plt

_DEBUG_ = False
# command to run 
# python cmty.py graph.txt -g -d -b -m -c
# -g means draw graph after built
# -d means debug mode
# -b means plot betweenes values
# -m means plot modularity values
# -c means draw community structure

#this method just reads the graph structure from the file
def buildG(G, file_, delimiter_):
	#construct the weighted version of the contact graph from cgraph.dat file
	#reader = csv.reader(open("/home/kazem/Data/UCI/karate.txt"), delimiter=" ")
	reader = csv.reader(open(file_), delimiter=delimiter_)
	for line in reader:
		if len(line) > 2:
			if float(line[2]) != 0.0:
				#line format: u,v,w
				G.add_edge(int(line[0]),int(line[1]),weight=float(line[2]))
		else:
			#line format: u,v
			G.add_edge(int(line[0]),int(line[1]),weight=1.0)

#keep removing edges from Graph until one of the connected components of Graph splits into two
#compute the edge betweenness
def CmtyGirvanNewmanStep(G, bwt_values):
	global _DEBUG_
	if _DEBUG_:
		print "Calling CmtyGirvanNewmanStep"
	init_ncomp = nx.number_connected_components(G)    #no of components
	ncomp = init_ncomp
	while ncomp <= init_ncomp:
		bw = nx.edge_betweenness_centrality(G)    #edge betweenness for G
		#find the edge with max centrality
		if len(bw.values()) == 0:
				break
		max_ = max(bw.values())
		bwt_values.append(max_)
		#find the edge with the highest centrality and remove all of them if there is more than one!
		for k, v in bw.iteritems():
				if float(v) == max_:
						G.remove_edge(k[0],k[1])    #remove the central edge
		ncomp = nx.number_connected_components(G)    #recalculate the no of components

#compute the modularity of current split
def _GirvanNewmanGetModularity(G, deg_, m_):
	global _DEBUG_
	New_A = nx.adj_matrix(G)
	New_deg = {}
	New_deg = UpdateDeg(New_A, G.nodes())
	#Let's compute the Q
	comps = nx.connected_components(G)    #list of components    
	# print 'No of communities in decomposed G: %d' % nx.number_connected_components(G)
	Mod = 0    #Modularity of a given partitionning
	for c in comps:
		EWC = 0    #no of edges within a community
		RE = 0    #no of random edges
		for u in c:
			EWC += New_deg[u]
			RE += deg_[u]        #count the probability of a random edge
		Mod += ( float(EWC) - float(RE*RE)/float(2*m_) )
	Mod = Mod/float(2*m_)
	if _DEBUG_:
		print "Modularity: %f" % Mod
	return Mod

def UpdateDeg(A, nodes):
	deg_dict = {}
	n = len(nodes)  #len(A) ---> some ppl get issues when trying len() on sparse matrixes!
	B = A.sum(axis = 1)
	for i in range(n):
		deg_dict[nodes[i]] = B[i, 0]
	return deg_dict

#run GirvanNewman algorithm and find the best community split by maximizing modularity measure
def runGirvanNewman(G, Orig_deg, m_, modularities, bwt_values):
	#let's find the best split of the graph
	BestQ = 0.0
	Q = 0.0
	Bestcomps_t = []
	best_graph = nx.Graph()
	while True:    
		CmtyGirvanNewmanStep(G, bwt_values)
		Q = _GirvanNewmanGetModularity(G, Orig_deg, m_);
		modularities.append(Q)
		#print "Modularity of decomposed G: %f" % Q
		if Q > BestQ:
				BestQ = Q
				#Bestcomps = nx.connected_components(G)    #Best Split
				best_graph = G.copy()
		if G.number_of_edges() == 0:
				break
	return best_graph

def main(argv):
	global _DEBUG_
	draw_graph = False
	draw_community = False
	plot_betweenness = False
	plot_modularity = False
	short_cut_modularity = False
	_DEBUG_ = False

	if len(argv) < 2:
		sys.stderr.write("Usage: %s <input graph>\n" % (argv[0],))
		return 1
	graph_fn = argv[1]
	
	for a in argv:
		if a == '-d':
			_DEBUG_ = True
		elif a == '-g':
			draw_graph = True
		elif a == '-b':
			plot_betweenness = True
		elif a == '-m':
			plot_modularity = True
		elif a == '-c':
			draw_community = True

	G = nx.Graph()  #let's create the graph first
	
	start = time.time()
	buildG(G, graph_fn, ',')
	end = time.time()
	print 'graph builded in %f' % (end - start), ' seconds'

	if _DEBUG_:
		print 'G nodes:', G.nodes()
		print 'G no of nodes:', G.number_of_nodes()

	if draw_graph:
		nx.draw(G)
		plt.show()

	n = G.number_of_nodes()    #|V|
	A = nx.adj_matrix(G)    #adjacenct matrix

	m_ = 0.0    #the weighted version for number of edges
	for i in range(0,n):
		for j in range(0,n):
			m_ += A[i,j]
	m_ = m_/2.0
	if _DEBUG_:
		print "m: %f" % m_

	#calculate the weighted degree for each node
	Orig_deg = {}
	Orig_deg = UpdateDeg(A, G.nodes())

	modularities = []
	bwt_values = []
	#run Newman alg
	start = time.time()
	best_graph = runGirvanNewman(G, Orig_deg, m_, modularities, bwt_values)
	end = time.time()
	print 'algo run in %f' % (end - start), ' seconds'
	
	if draw_community:
		nx.draw(best_graph)
		plt.show()

	if plot_betweenness:
		plt.plot(bwt_values)
		plt.ylabel('betweenness value')
		plt.xlabel('iteration')
		plt.show()
	
	if plot_modularity:
		plt.plot(modularities)
		plt.ylabel('modularity value')
		plt.xlabel('iteration')
		plt.show()

if __name__ == "__main__":
		sys.exit(main(sys.argv))
