import networkx as nx
from networkx.readwrite import json_graph
from scipy.io import loadmat
graphData = loadmat('MulNet.mat')

#this function returns a directed graph
def txt2Graph(FielName):
	Graphtype=nx.DiGraph()
	g= nx.read_edgelist(FielName,create_using=Graphtype,nodetype=int)
	for i in range(1, 6408):
		if i not in g.nodes:
			g.add_node(i)

	H = nx.DiGraph()
	H.add_nodes_from(sorted(g.nodes(data=True)))
	print("hereeee:  ", H.number_of_nodes())
	return graphe_TO_json(H,FielName) 

def graphe_TO_json(g,FielName):
    
    data =  json_graph.node_link_data(g,{"link": "links", "source": "source", "target": "target","weight":"weight"})
    
    data['nodes'] = [ {"id": i,"state":"non_infected","Protector":"false","opinion":"normal","beta":0,"omega":0,"delta":0,"jug":0,"Infetime":0,"AccpR":0,"SendR":0,"Accp_NegR":0,"value":0,"blocked":'false',"degree":list(graphData['NodeDegree'][i-1]),"neighbors":[n for n in g.successors(i)]} for i in range(1,len(data['nodes'])+1) ]
    data['links'] = [ {"source":u,"target":v,"weight":(g.degree[u]+g.degree[v])/2} for u,v in g.edges ]
    return data
if __name__ == '__main__':
	#ytGraph =txt2Graph("yt.txt")
	twGraph =txt2Graph("tw.txt")
	#ffGraph =txt2Graph("ff.txt")
	

	print(twGraph['nodes'][0]['degree'])