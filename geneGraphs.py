import networkx as nx
from networkx.readwrite import json_graph
from scipy.io import loadmat
graphData = loadmat('MulNet.mat')

#this function returns a directed graph
def txt2Graph(FielName):
	H = nx.path_graph(6408)
	g=nx.DiGraph()
	g.add_nodes_from(H)
	g.remove_node(0)
	Graphtype=nx.DiGraph()
	g0= nx.read_edgelist(FielName,create_using=Graphtype,nodetype=int)
	g.add_edges_from(g0.edges())

	#print("hereeee:  ", len(list(g.successors(2787))))
	return graphe_TO_json(g,FielName) 

def graphe_TO_json(g,FielName):
    
    data =  json_graph.node_link_data(g,{"link": "links", "source": "source", "target": "target","weight":"weight"})
    
    data['nodes'] = [ {"id": i,"state":"non_infected","Protector":"false","opinion":"normal","beta":0,"omega":0,"delta":0,"jug":0,"Infetime":0,"AccpR":0,"SendR":0,"Accp_NegR":0,"value":0,"blocked":'false',"degree":list(graphData['NodeDegree'][i-1]),"neighbors":[n for n in g.successors(i)]} for i in range(1, len(data['nodes'])+1) ]
    data['links'] = [ {"source":u,"target":v,"weight":(g.degree[u]+g.degree[v])/2} for u,v in g.edges ]
    return data
if __name__ == '__main__':
	#ytGraph =txt2Graph("yt.txt")
	twGraph =txt2Graph("tw.txt")
	#ffGraph =txt2Graph("ff.txt")


	print(twGraph['nodes'][2786]['degree'])
	print(graphData['NodeDegree'][2786])	