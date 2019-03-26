import matplotlib.pyplot as plt
import networkx as nx
import json
import math
import numpy as np

def index(h,aa):
    return h*21+aa

aminoacids = "ACDEFGHIKLMNPQRSTVWY"
colors = ["#777775", "#fedd00", "#ef3340", "#ef3340", "#000000", "#fedd00", "#0087c7", "#333334", "#0087c7", "#333334", "#333334", "#65428a", "#fedd00", "#65428a", "#0087c7", "#0087c7", "#333334", "#777775", "#000000", "#000000", "#ffffff"]    
textcolors = ["#ffffff", "#000000", "#000000", "#000000", "#000000", "#000000", "#000000", "#000000", "#000000", "#000000", "#000000", "#000000", "#000000", "#000000", "#000000", "#000000", "#000000", "#000000", "#000000", "#000000", "#000000"]    
hiddenstates = [11,1,4,14,19]

with open("ratenetwork.json") as json_file:  
    modelparams = json.load(json_file)
    
hiddenfreqs = modelparams["hiddenfreqs"]

plt.figure(figsize=(10,10))
G = nx.Graph()
labels = {}
fontcolours = {}
for h in hiddenstates:
    G.add_node(index(h,20))
    G.nodes[index(h,20)]['radius'] = 0.0
    G.nodes[index(h,20)]['textcolor'] = "#ffffff"
    labels[index(h,20)] = ""
    for aa in range(20):
        node = G.add_node(index(h,aa))        
        radius =  math.sqrt(hiddenfreqs[h]*modelparams["aafreqs_h%d" % (h+1)][aa])
        G.nodes[index(h,aa)]['radius'] = radius
        G.nodes[index(h,aa)]['textcolor'] = textcolors[aa]
        fontcolours[index(h,aa)] = textcolors[aa]
        if radius < 0.02:
            labels[index(h,aa)] = ""
            G.nodes[index(h,aa)]['radius'] = 0.0
        else:
            labels[index(h,aa)] = str(aminoacids[aa])
        

edgelist = []
for h1  in hiddenstates:
    for aa1 in range(20):
        G.add_edge(index(h1,20), index(h1,aa1), weight=400.0)
        for aa2 in range(aa1+1,20):
            freq1 = modelparams["aafreqs_h%d" % (h1+1)][aa1]
            freq2 = modelparams["aafreqs_h%d" % (h1+1)][aa2]
            edgeweight = modelparams["aa_exchangeablities"][aa1][aa2]*(freq1+freq2)*5.0
            if  G.nodes[index(h1,aa1)]['radius'] > 0.0 and G.nodes[index(h1,aa2)]['radius'] > 0.0:
                G.add_edge(index(h1,aa1), index(h1,aa2), weight=edgeweight)
                edgelist.append((index(h1,aa1),index(h1,aa2)))


h_edgelist = []
for h1 in hiddenstates:
    for h2 in hiddenstates:
        if h1 < h2:
             G.add_edge(index(h1,20), index(h2,20), weight=modelparams["transitionrates"][h1][h2]*(hiddenfreqs[h1]+hiddenfreqs[h2])*20000.0)
             h_edgelist.append((index(h1,20), index(h2,20)))
        
"""    
h_edgelist = []
for h1 in hiddenstates:
    for h2 in hiddenstates:
        if h1 < h2:
            for aa1 in range(20):
                weight = modelparams["transitionrates"][h1][h2]*(hiddenfreqs[h1]+hiddenfreqs[h2])
                G.add_edge(index(h1,aa1), index(h2,aa1), weight=weight*200.0)
                h_edgelist.append((index(h1,aa1), index(h2,aa1)))
"""                

pos = nx.spring_layout(G)
for h in hiddenstates:
    vs = np.zeros(2)
    for aa in range(20):
        ind = index(h,aa)
        n = 0.0
        if ind in pos:
            vs += pos[ind]
            n += 1.0
    #pos[index(h,20)] = vs/n
            
        

#for node in G.nodes:
    #print("node: ",labels.get(node,""))
nx.draw_networkx_labels(G, pos, labels=labels, font_size=16, font_color="#eeeeee", font_weight="bold")

nx.draw_networkx_nodes(G, pos,
                       nodelist=G.nodes,
                       node_color=[colors[index%21] for index in G.nodes],
                       node_size=[float(G.nodes[node]['radius'])*20000.0 for node in G.nodes],
                       alpha=0.9)
nx.draw_networkx_edges(G, pos, edgelist=edgelist, width=[G[edge[0]][edge[1]]['weight'] for edge in edgelist], alpha=0.5, edge_color='b')
#nx.draw_networkx_edges(G, pos, edgelist=h_edgelist, width=[G[edge[0]][edge[1]]['weight']*2.0 for edge in edgelist], alpha=0.5, edge_color='g')

plt.axis('off')
#plt.show()
plt.savefig("ratenetwork.svg")