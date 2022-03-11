import numpy as np
from networkx.readwrite import json_graph
import networkx as nx
import matplotlib.pyplot as plt
import random
import time
#from matplotlib.ticker import NullFormatter  
import multiprocessing 
from multiprocessing import Manager
import math
import geneGraphs
def Neighbour_finder(g,new_active):
    
    targets = []
    for node in new_active:
        targets += g.neighbors(node)
    return(targets)
def MyUnion(lst1, lst2):
    final_list = list(set(lst1) | set(lst2))
    return final_list
def HISBmodel (Graph,graphTypeIndice,allNodesInfected,GlobalRumorPopularity,Seed_Set,Opinion_Set,Statistical,paramater,K,timeOfTheRumorDetection,method):
    
    #Opinion:normal/denying/supporting
    #State:non_infected/infected/spreaders 
    #Statistical:{'NonInfected':NbrOFnodes,'Infected':**,'Spreaders':**,OpinionDenying':**,'OpinionSupporting':**,'RumorPopularity':**}
    bl=0
    ListInfectedNodes=Seed_Set[:]
    Opinion_Set=Opinion_Set[:]
    time=0.125
    Probability=1
    i=0
    #Initialis Parameters----------------------------
    #-------------------------
    Nbr_Spreaders=len(ListInfectedNodes)  #why: Nbr_Spreaders=len(ListInfectedNodes)
    Nbr_nonInfected=6407
    Nbr_Infected=len(ListInfectedNodes)
    OpinionDenying=0
    OpinionSupporting=0
    RumorPopularity=0
    InitParameters(Graph,paramater)
    ''' if(L_protector!=None):
       for each  in L_protector:
        Graph.nodes[each]['jug']=1'''
    layerIndice = graphTypeIndice%3  #which is in reality the id of the current simulation, it help us to know if it's a yt, ff or a tw simulation
  
    for each  in ListInfectedNodes:
      if each!=0:
        Graph.nodes[each]['Infetime']=0.125 
        Graph.nodes[each]['state']='spreaders'
        Graph.nodes[each]['AccpR']+=1
        RumorPopularity+=Graph.nodes[each]['degree'][layerIndice]
        Nbr_Infected+=1
        Nbr_nonInfected-=1
        if (Opinion_Set[i]=='denying'):

            Graph.nodes[each]['opinion']='denying'
            Graph.nodes[each]['Accp_NegR']+=1
            OpinionDenying+=1
        else:
          Graph.nodes[each]['opinion']='supporting'
          OpinionSupporting+=1
        i+=1
    for each in allNodesInfected:
        GlobalRumorPopularity+=Graph.nodes[each]['degree'][layerIndice]
        
    
    #------------------------------
    Statistical.append({'NonInfected':Nbr_nonInfected,'Infected':Nbr_Infected,'graphTypeIndice':graphTypeIndice, 'AllInfected':len(allNodesInfected), 'Spreaders':Nbr_Spreaders,'OpinionDenying':OpinionDenying,'OpinionSupporting':OpinionSupporting,'RumorPopularity':RumorPopularity,'GlobalRumorPopularity':GlobalRumorPopularity,'graph':0})
    #----------------------

    #if the list is empty we stop the propagation
    
    while ListInfectedNodes:
      #add infection coming from the other networks
      for node in allNodesInfected:
        if not node in ListInfectedNodes:
          ListInfectedNodes.append(node)
      RumorPopularity = 0
      Nbr_Spreaders = 0
      L=len(ListInfectedNodes)
      for X in reversed(range(0,L)):
        
        id = ListInfectedNodes[X]
        
        #relative time of rumor spreading: Verify is the current infected node will lose his intrest on the rumor or note
        RelativeTime = time - Graph.nodes[id]['Infetime'] 
        if (np.exp(-RelativeTime * Graph.nodes[id]['beta']) < 0.15) :   #give a time to the infection
          ListInfectedNodes.pop(X)
          Graph.nodes[id]['state'] = "infected"
          
              

        else:
            #atrraction of nodes
            ActualAttraction = np.exp(-RelativeTime * Graph.nodes[id]['beta']) * np.abs(np.sin((RelativeTime * Graph.nodes[id]['omega'] )+ Graph.nodes[id]['delta']))
            
            RumorPopularity += ActualAttraction * Graph.nodes[id]['degree'][layerIndice]
            GlobalRumorPopularity += ActualAttraction * Graph.nodes[id]['degree'][layerIndice]
            #rumor spreading
            
            c=np.random.random_sample() # a random between 0 and 1
            
            sumDegree = sum(Graph.nodes[id]['degree'])
            probSpreadingInThisLayer =  Graph.nodes[id]['degree'][layerIndice]/sumDegree
            if (c<=ActualAttraction * probSpreadingInThisLayer):
            #if (c<=ActualAttraction):  # yes the current node will speard the rumor
                Nbr_Spreaders+=1
                Graph.nodes[id]['state']='spreaders'
                
                #Calculating if any nodes of those neighbours can be activated, if yes add them to new_ones.
                success = np.random.uniform(0,1,len(Graph.nodes[id]['neighbors'])) < Probability #choic alpha nodes  ### this return list: false true flase .... 
                # success == [ True  True  True False  True .... True False False  True False]                
                new_ones = list(np.extract(success, sorted(Graph.nodes[id]['neighbors'])))  # to avoid the case of sending to all the neighbors
                
                Graph.nodes[id]['SendR']+=len(new_ones)
                
                #Sending Rumor
                for each in new_ones:
                    #Accepted Rumor Probability 
                    ProbToAccRumor = Graph.nodes[id]['degree'][layerIndice]/ (Graph.nodes[id]['degree'][layerIndice] + Graph.nodes[each]['degree'][layerIndice])*0.3
                    if (Graph.nodes[each]['blocked'] =='false'):
                        if(np.random.random_sample()<=ProbToAccRumor):
                        
                            Graph.nodes[each]['AccpR']+=1
    
                            if (Graph.nodes[each]['Infetime']==0 ):
                                Nbr_Infected+=1
                                Nbr_nonInfected-=1
                                Graph.nodes[each]['Infetime'] =time
                                Graph.nodes[each]['opinion'] =Graph.nodes[id]['opinion']        # !!! why
                                
                                if each != 0:
                                  ListInfectedNodes.append(each)
                                if (Graph.nodes[each]['opinion']=="denying"):
                                    #negativ opinion
                                    Graph.nodes[each]['Accp_NegR']+=1
                                    OpinionDenying+=1
                                else:
                                     OpinionSupporting+=1
                            elif (Graph.nodes[id]['opinion']=="denying"):
                                Graph.nodes[each]['Accp_NegR']+=1
                    
                   # else:
                        #print('le noeud',each,'is blocked')
                        #updateOpinion(id)'''
                if (Graph.nodes[id]['opinion']=="denying"):
                    OpinionDenying-=1
                else:
                    OpinionSupporting-=1
                Graph.nodes[id]['opinion']= updateOpinion(jug=Graph.nodes[id]['jug'],Accpet_NegR=Graph.nodes[id]['Accp_NegR'],Nbr_OF_R=Graph.nodes[id]['AccpR'],Role=Graph.nodes[id]['Protector'])
                if (Graph.nodes[id]['opinion']=="supporting"):
                    
                    OpinionSupporting+=1
                else:
                    
                    OpinionDenying+=1       
      
      #save each step to send it to viewing later
      Statistical.append({'NonInfected':Nbr_nonInfected,'Infected':Nbr_Infected, 'AllInfected':len(allNodesInfected), 'Spreaders':Nbr_Spreaders,'OpinionDenying':OpinionDenying,'OpinionSupporting':OpinionSupporting,'RumorPopularity':RumorPopularity,'GlobalRumorPopularity':GlobalRumorPopularity,'graph':0})
      
      allNodesInfected = allNodesInfected.union(ListInfectedNodes)
      
      if time >=timeOfTheRumorDetection*0.125 and bl<K  and method != 'NP' :     #############################
          print(method," At time:", time, "blocked nodes Nbr:", bl)
          #Nodes=len(Graph.nodes)
          
          p=K-bl
          if (method=='BNLS'):
              Random_Blocking_nodes(Graph, p)
              bl=len(blocked(Graph))
          if (method=='BNLSM'):
              Degree_MAX_Blocking_nodes(Graph,p)
              bl=len(blocked(Graph))
          if (method=='BNLSCen'):
              Centrality_Blocking_nodes(Graph,p)
              bl=len(blocked(Graph))
          if (method=='BNLSBeta'):
              Beta_Blocking_nodes(Graph,p)
              bl=len(blocked(Graph))
          if (method=='BNLSBetaD'):
              BetaD_Blocking_nodes(Graph,p)
              bl=len(blocked(Graph))
                                                            #############################
              
          elif method=='TCS':
              Random_TRuth_comp(Graph, p)
              Liste_protector=Protector(Graph)
              taille=len(Liste_protector)
             # print(taille)
              for i in range(taille):
                  Graph.nodes[Liste_protector[i]]['Infetime'] =time
                  Graph.nodes[Liste_protector[i]]['opinion']=="denying"
                  Graph.nodes[Liste_protector[i]]['state']='spreaders'
                  Graph.nodes[Liste_protector[i]]['AccpR']+=1
                  Graph.nodes[Liste_protector[i]]['Accp_NegR']+=1
                  ListInfectedNodes.append(Liste_protector[i])
              bl=len(Liste_protector)
          elif method=='TCSM':
             MaxDegree_TRuth_comp(Graph, p)
             Liste_protector=Protector(Graph)
             taille=len(Liste_protector)
             for i in range(taille):
                 ListInfectedNodes.append(Liste_protector[i])
                 Graph.nodes[Liste_protector[i]]['Infetime'] =time
                 Graph.nodes[Liste_protector[i]]['opinion']=="denying"
                 Graph.nodes[Liste_protector[i]]['state']='spreaders'
                 Graph.nodes[Liste_protector[i]]['AccpR']+=1
                 Graph.nodes[Liste_protector[i]]['Accp_NegR']+=1
             bl=len(Liste_protector)
          elif method=='TCSCen':
             Centrality_TRuth_comp(Graph, p)
             Liste_protector=Protector(Graph)
             taille=len(Liste_protector)
             for i in range(taille):
                 ListInfectedNodes.append(Liste_protector[i])
                 Graph.nodes[Liste_protector[i]]['Infetime'] =time
                 Graph.nodes[Liste_protector[i]]['opinion']=="denying"

                 Graph.nodes[Liste_protector[i]]['state']='spreaders'
                 Graph.nodes[Liste_protector[i]]['AccpR']+=1
                 Graph.nodes[Liste_protector[i]]['Accp_NegR']+=1
             bl=len(Liste_protector)
          elif method=='TCSBeta':
             Beta_TRuth_comp(Graph, p)
             Liste_protector=Protector(Graph)
             taille=len(Liste_protector)
             for i in range(taille):
                 ListInfectedNodes.append(Liste_protector[i])
                 Graph.nodes[Liste_protector[i]]['Infetime'] =time
                 Graph.nodes[Liste_protector[i]]['opinion']=="denying"
                 Graph.nodes[Liste_protector[i]]['state']='spreaders'
                 Graph.nodes[Liste_protector[i]]['AccpR']+=1
                 Graph.nodes[Liste_protector[i]]['Accp_NegR']+=1
             bl=len(Liste_protector)
          elif method=='TCSBetaD':
             BetaD_TRuth_comp(Graph, p)
             Liste_protector=Protector(Graph)
             taille=len(Liste_protector)
             for i in range(taille):
                 ListInfectedNodes.append(Liste_protector[i])
                 Graph.nodes[Liste_protector[i]]['Infetime'] =time
                 Graph.nodes[Liste_protector[i]]['opinion']=="denying"
                 Graph.nodes[Liste_protector[i]]['state']='spreaders'
                 Graph.nodes[Liste_protector[i]]['AccpR']+=1
                 Graph.nodes[Liste_protector[i]]['Accp_NegR']+=1
             bl=len(Liste_protector)
          print(method," At time:", time, "blocked nodes Nbr:", bl)
            
              
      time += 0.25;   
def InitParameters(Graph,parameters):
    #Individual back ground knowledge:Beta
    #Forgetting and remembering factore:Omega
    #Hesitating factore:Deleta
    #Subjective judjement:Jug
  
    for node in Graph.nodes:
       Graph.nodes[node]['omega']=Inclusive(parameters[0]['omega_min'],parameters[0]['omega_max'])
       Graph.nodes[node]['beta']=Inclusive(parameters[0]['beta_min'],parameters[0]['beta_max'])
       Graph.nodes[node]['delta']=Inclusive(parameters[0]['delta_min'],parameters[0]['delta_max'])
       Graph.nodes[node]['jug']=Inclusive(parameters[0]['Jug_min'],parameters[0]['Jug_max'])
def Inclusive(min,max):
   
   b= ((np.random.random_sample()*(max - min )) + min)
    
   return b
def updateOpinion(jug,Accpet_NegR,Nbr_OF_R,Role): 
    if(Role=='True'):
       return 'denying' 
   
    opinion=jug
    if Accpet_NegR != 0:
        opinion*=(Accpet_NegR / Nbr_OF_R)
   
    
    if(np.random.random_sample()<= opinion):
        return 'denying'
    else:
        return 'supporting'
def graphe_TO_json(g):
    
    data =  json_graph.node_link_data(g,{"link": "links", "source": "source", "target": "target","weight":"weight"})
    data['nodes'] = [ {"id": i,"state":"non_infected","Protector":"false","opinion":"normal","beta":0,"omega":0,"delta":0,"jug":0,"Infetime":0,"AccpR":0,"SendR":0,"Accp_NegR":0,"value":0,"blocked":'false',"degree":g.degree[i],"neighbors":[n for n in g.neighbors(i)]} for i in range(len(data['nodes'])) ]
    data['links'] = [ {"source":u,"target":v,"weight":(g.degree[u]+g.degree[v])/2} for u,v in g.edges ]
    return data
def geneList_Infectede(Listinfected,Listopinion,N,percentage):
    #10% of Popularity is infected 
    Nbr_OF_ndodesI=int(N*percentage/100)
    L=list(range(1,N+1))
    List=random.sample(L, Nbr_OF_ndodesI)
    opinion=np.random.uniform(0,1,Nbr_OF_ndodesI)
    
    for each in range(Nbr_OF_ndodesI):
        Listinfected.append(List[each])
        if opinion[each]<=0.2:
           Listopinion.append('denying')
        else:
            Listopinion.append('supporting')
            
def parameters(parameter,stepBeta=1,Beta=0.2,stepOmega=5.2,Omega=math.pi/3,stepDelta=0.65,Delta=math.pi/24,stepJug=0.6,Jug=0.1):
    Beta_max=Beta+stepBeta
    Omega_max=Omega +stepOmega
    Delta_max=Delta +stepDelta
    Jug_max=Jug+stepJug
    parameter.append({'beta_min':round(Beta,2),'beta_max':round(Beta_max,2),'omega_min':round(Omega,2),'omega_max':round(Omega_max,2),'delta_min':round(Delta,2),'delta_max':round(Delta_max,2),'Jug_min':round(Jug,2),'Jug_max':round(Jug_max,2)})
def Start(i,Graph,allNodesInfected,GlobalRumorPopularity,parameter,Stat,percentage,K,timeOfTheRumorDetection,method):
    print("The ", int(i/3)+1, "th simulation")
    for each in range(1, len(Graph.nodes)+1):
        Graph.nodes[each]['opinion']="normal"
        Graph.nodes[each]['Infetime']=0 
        Graph.nodes[each]['state']='non_infected'
        Graph.nodes[each]['Protector']='false'
        Graph.nodes[each]['blocked']='false'
        
    Statistical=[]
    ListInfected=[]
    Listopinion=[]
    if i%3==0:
      allNodesInfected=set()

      GlobalRumorPopularity=0
    #X% of Popularity is infected 
    geneList_Infectede(ListInfected,Listopinion,6407,percentage)
    for n in ListInfected:
      allNodesInfected.add(n)
   
    HISBmodel(Graph,i,allNodesInfected,GlobalRumorPopularity,ListInfected,Listopinion,Statistical,parameter,K,timeOfTheRumorDetection,method)  
    Stat.append(Statistical)    
    
    
def globalStat(S,Stat_Global,GlobalStatYt,GlobalStatTw,GlobalStatFf,parameter,method):
    max=0
    Stat=[]
    for each in S:
        
        L=len(each)
        Stat.append(each)
        if(L>max):
            max=L
    for i in range(len(Stat)):
        L=len(Stat[i])
        Nbr_nonInfected=Stat[i][L-1]['NonInfected']
        Nbr_Infected=Stat[i][L-1]['Infected']
        Nbr_Spreaders=Stat[i][L-1]['Spreaders']
        OpinionDenying=Stat[i][L-1]['OpinionDenying']
        OpinionSupporting=Stat[i][L-1]['OpinionSupporting']
        RumorPopularity=Stat[i][L-1]['RumorPopularity']
        GlobalRumorPopularity=Stat[i][L-1]['GlobalRumorPopularity']
        AllInfected=Stat[i][L-1]['AllInfected']

        for j in range(L,max):
            Stat[i].append({'NonInfected':Nbr_nonInfected,'Infected':Nbr_Infected,'AllInfected':AllInfected, 'Spreaders':Nbr_Spreaders,'OpinionDenying':OpinionDenying,'OpinionSupporting':OpinionSupporting,'RumorPopularity':RumorPopularity,'GlobalRumorPopularity':GlobalRumorPopularity,'graph':0})       

    y1=[]
    y2=[]
    y3=[]
    y4=[]
    y5=[]
    y6=[]
    y8=[]
    z=list()
    for i in range(3):
      z.append([list(),list(),list()])
    Len=len(Stat)
  
    for i in range(max):
        
        Infected=0
        Spreaders=0
        RumorPopularity=0
        GlobalRumorPopularity=0
        OpinionDenying=0
        OpinionSupporting=0
        AllInfected=0
        infectedInNetworks=[0]*3
        spreadersINetworks=[0]*3
        RPopularityINetworks=[0]*3
        counter=0
        print("number of simulations: ", len(Stat))
        for each in Stat: 
            counter+=1          
            Infected+=(each[i]['Infected'])
            Spreaders+=(each[i]['Spreaders'])
            RumorPopularity+=(each[i]['RumorPopularity'])
            GlobalRumorPopularity+=(each[i]['GlobalRumorPopularity'])
            OpinionDenying+=(each[i]['OpinionDenying'])
            OpinionSupporting+=(each[i]['OpinionSupporting'])
            AllInfected+=(each[i]['AllInfected'])
            for j in range(3):
              if counter % 3 == j:
                infectedInNetworks[j]+=(each[i]['Infected'])
                spreadersINetworks[j]+=(each[i]['Spreaders'])
                print("thiiiiiiiiiiiiiiiiiiiiiiis spreaders: ", each[i]['Spreaders'])
                RPopularityINetworks[j]+=(each[i]['RumorPopularity'])
        y1.append(Infected/Len)
        y2.append(Spreaders/Len)
        y3.append(RumorPopularity/Len)
        y8.append(GlobalRumorPopularity/Len)
        y4.append(OpinionDenying/Len)
        y5.append(OpinionSupporting/Len)
        y6.append(AllInfected/Len)
        #the first list in Z contrain the lists of infected, spreaders and RumorPopularity of youtube, the second for twitter...
        for j in range(3):
                z[j][0].append(infectedInNetworks[j]/int(Len/3))
                z[j][1].append(spreadersINetworks[j]/int(Len/3))
                z[j][2].append(RPopularityINetworks[j]/int(Len/3))
        
    GlobalStatYt.append({'Infected':z[1][0],'Spreaders':z[1][1],'RumorPopularity':z[1][2],'max':max,'method':method})
    GlobalStatTw.append({'Infected':z[0][0],'Spreaders':z[0][1],'RumorPopularity':z[0][2],'max':max,'method':method})
    GlobalStatFf.append({'Infected':z[2][0],'Spreaders':z[2][1],'RumorPopularity':z[2][2],'max':max,'method':method})
    Stat_Global.append({'Infected':y1,'Spreaders':y2,'RumorPopularity':y3,'GlobalRumorPopularity':y8,'OpinionDenying':y4,'OpinionSupporting':y5, 'AllInfected':y6, 'parameter':parameter,'max':max,'method':method})       
    

def Display(Stat_Global,GlobalStatYt,GlobalStatTw,GlobalStatFf,xx,title_fig,nb):
    print("Stat_Global  : ", Stat_Global)
    Title=['NP']
    max=0
    Stat=[]
    Infected=[]
    para=[]
    for each in Stat_Global:
        L=each['max']   #It's the bigger number of observations of the done simulations
        
        para.append(each['method'])
        metho=str(each['method'])
        if metho.startswith('TCS'):
            Infected.append(each['OpinionSupporting'][L-1]/Nodes)
        else:
            Infected.append(each['Infected'][L-1]/Nodes)
                
        
        Stat.append(each)
        if(L>max):
            max=L
    for each in Stat:
        L=each['max']
        if (L<max):
            Nbr_Infected=each['Infected'][L-1]
            Nbr_Spreaders=each['Spreaders'][L-1]
            OpinionDenying=each['OpinionDenying'][L-1]
            OpinionSupporting=each['OpinionSupporting'][L-1]
            RumorPopularity=each['RumorPopularity'][L-1]
            GlobalRumorPopularity=each['GlobalRumorPopularity'][L-1]
            allNodesInfected=each['AllInfected'][L-1]
            for j in range(L,max):
                each['Infected'].append(Nbr_Infected)
                each['Spreaders'].append(Nbr_Spreaders)
                each['OpinionDenying'].append(OpinionDenying)
                each['OpinionSupporting'].append(OpinionSupporting)
                each['RumorPopularity'].append(RumorPopularity)
                each['GlobalRumorPopularity'].append(GlobalRumorPopularity)
                each['AllInfected'].append(allNodesInfected)

    pro=int(max/50)
    if max%2 == 1:
      max-=1
    for each in Stat:
            for j in reversed(range(max)):
                d=j%pro
                if(d!=0):
                    each['Infected'].pop(j)
                    each['Spreaders'].pop(j)
                    each['OpinionDenying'].pop(j)
                    each['OpinionSupporting'].pop(j)
                    each['RumorPopularity'].pop(j)
                    each['GlobalRumorPopularity'].pop(j)
                    each['AllInfected'].pop(j)

    for each in Stat:
            for j in reversed(range(10)):
                    each['Infected'].pop(20+j)
                    each['Spreaders'].pop(20+j)
                    each['OpinionDenying'].pop(20+j)
                    each['OpinionSupporting'].pop(20+j)
                    each['RumorPopularity'].pop(20+j)
                    each['GlobalRumorPopularity'].pop(20+j)
                    each['AllInfected'].pop(20+j)
    x = range(0,len(Stat[0]['AllInfected']))
    x=np.array(x)*pro
    #for the yt,tw and ff: ------------------------------------------------------------
    #for yt
    maxyt=0
    StatYt=[]
    InfectedYt=[]
    paraYt=[]

    for each in GlobalStatYt:
        L=each['max']   #It's the bigger number of observations of the done simulations
        paraYt.append(each['method'])
        metho=str(each['method'])
        if metho.startswith('TCS'):
            InfectedYt.append(each['OpinionSupporting'][L-1]/Nodes)
        else:
            InfectedYt.append(each['Infected'][L-1]/Nodes)
                
        
        StatYt.append(each)

        if(L>maxyt):
            maxyt=L
    for each in StatYt:
        L=each['max']
        if (L<maxyt):
            Nbr_Infected=each['Infected'][L-1]
            Nbr_Spreaders=each['Spreaders'][L-1]
            RumorPopularity=each['RumorPopularity'][L-1]
            for j in range(L,maxyt):
                each['Infected'].append(Nbr_Infected)
                each['Spreaders'].append(Nbr_Spreaders)
                each['RumorPopularity'].append(RumorPopularity)
    
    proYt=int(maxyt/50)
    if maxyt%2 == 1:
      maxyt-=1
    for each in StatYt:
            for j in reversed(range(maxyt)):
                d=j%proYt
                if(d!=0):
                    each['Infected'].pop(j)
                    each['Spreaders'].pop(j)
                    each['RumorPopularity'].pop(j)

    for each in StatYt:
            for j in reversed(range(10)):
                    each['Infected'].pop(20+j)
                    each['Spreaders'].pop(20+j)
                    each['RumorPopularity'].pop(20+j)
    xYt = range(0,len(StatYt[0]['Infected']))
    xYt=np.array(x)*proYt
    ##########
    #for tw
    maxTw=0
    StatTw=[]
    InfectedTw=[]
    paraTw=[]
    for each in GlobalStatTw:
        L=each['max']   #It's the bigger number of observations of the done simulations
        paraTw.append(each['method'])
        metho=str(each['method'])
        if metho.startswith('TCS'):
            InfectedTw.append(each['OpinionSupporting'][L-1]/Nodes)
        else:
            InfectedTw.append(each['Infected'][L-1]/Nodes)
                
        
        StatTw.append(each)
        if(L>maxTw):
            maxTw=L
    for each in StatTw:
        L=each['max']
        if (L<maxTw):
            Nbr_Infected=each['Infected'][L-1]
            Nbr_Spreaders=each['Spreaders'][L-1]
            RumorPopularity=each['RumorPopularity'][L-1]
            for j in range(L,maxTw):
                each['Infected'].append(Nbr_Infected)
                each['Spreaders'].append(Nbr_Spreaders)
                each['RumorPopularity'].append(RumorPopularity)
    
    proTw=int(maxTw/50)
    if maxTw%2 == 1:
      maxTw-=1
    for each in StatTw:
            for j in reversed(range(maxTw)):
                d=j%proTw
                if(d!=0):
                    each['Infected'].pop(j)
                    each['Spreaders'].pop(j)
                    each['RumorPopularity'].pop(j)

    for each in StatTw:
            for j in reversed(range(10)):
                    each['Infected'].pop(20+j)
                    each['Spreaders'].pop(20+j)
                    each['RumorPopularity'].pop(20+j)
    xTw = range(0,len(StatTw[0]['Infected']))
    xTw=np.array(x)*proTw
    #######
    #for Ff
    maxFf=0
    StatFf=[]
    InfectedFf=[]
    paraFf=[]
    for each in GlobalStatFf:
        L=each['max']   #It's the bigger number of observations of the done simulations
        paraFf.append(each['method'])
        metho=str(each['method'])
        if metho.startswith('TCS'):
            InfectedFf.append(each['OpinionSupporting'][L-1]/Nodes)
        else:
            InfectedFf.append(each['Infected'][L-1]/Nodes)
                
        
        StatFf.append(each)
        if(L>maxFf):
            maxFf=L
    for each in StatFf:
        L=each['max']
        if (L<maxFf):
            Nbr_Infected=each['Infected'][L-1]
            Nbr_Spreaders=each['Spreaders'][L-1]
            RumorPopularity=each['RumorPopularity'][L-1]
            for j in range(L,maxFf):
                each['Infected'].append(Nbr_Infected)
                each['Spreaders'].append(Nbr_Spreaders)
                each['RumorPopularity'].append(RumorPopularity)
    
    proFf=int(maxFf/50)
    if maxFf%2 == 1:
      maxFf-=1
    for each in StatFf:
            for j in reversed(range(maxFf)):
                d=j%proFf
                if(d!=0):
                    each['Infected'].pop(j)
                    each['Spreaders'].pop(j)
                    each['RumorPopularity'].pop(j)

    for each in StatFf:
            for j in reversed(range(10)):
                    each['Infected'].pop(20+j)
                    each['Spreaders'].pop(20+j)
                    each['RumorPopularity'].pop(20+j)
    xFf = range(0,len(StatFf[0]['Infected']))
    xFf=np.array(x)*proFf
    # plot 
    
    type=['x','*','p','8','h','H','.','+','4','1','2','3']
    
   
    ###################################################################################
    #AllInfected
    plt.figure(num=xx)
    plt.subplot()
    #k="{}:{},{}]" 
    k="{}" 
    for infected,j in zip( Stat,range(len(Stat))):
      quotients = [number/Nodes  for number in infected["AllInfected"]]
      plt.plot(x,quotients,marker=type[j],markersize=7,linewidth=1,label=k.format(Title[j]))
    plt.legend(fontsize=12)   

    plt.xlabel('Temps',fontsize=10)
    plt.ylabel('Nombre des individues')
    plt.title("infected in the whole popularity")
    plt.grid(True)
    plt.savefig(title_fig+'infectedInWholePopularity.pdf',dpi=50)
    #Infected yt
    xx+=1
    plt.figure(num=xx)
    plt.subplot()
    #k="{}:{},{}]" 
    k="{}" 
    for infected,j in zip( StatYt,range(len(StatYt))):
            quotients = [number /663  for number in infected["Infected"]]
            plt.plot(x,quotients,marker=type[j],markersize=7,linewidth=1,color="red",label=k.format(Title[j]))
    plt.legend(fontsize=12)   

    plt.xlabel('Temps',fontsize=10)
    plt.ylabel('Nombre des individues')
    plt.title("infected in youtube")
    plt.grid(True)
    plt.savefig(title_fig+'infectedYoutube.pdf',dpi=50)
    #Infected tw
    xx+=1
    plt.figure(num=xx)
    plt.subplot()
    #k="{}:{},{}]" 
    k="{}" 
    for infected,j in zip( StatTw,range(len(StatTw))):
            quotients = [number /5540  for number in infected["Infected"]]
            plt.plot(x,quotients,marker=type[j],markersize=7,linewidth=1,color="blue",label=k.format(Title[j]))
    plt.legend(fontsize=12)   

    plt.xlabel('Temps',fontsize=10)
    plt.ylabel('Nombre des individues')
    plt.title("infected in twitter")
    plt.grid(True)
    plt.savefig(title_fig+'infectedTwitter.pdf',dpi=50)
    #Infected ff
    xx+=1
    plt.figure(num=xx)
    plt.subplot()
    #k="{}:{},{}]" 
    k="{}" 
    for infected,j in zip( StatFf,range(len(StatFf))):
            quotients = [number /5702  for number in infected["Infected"]]
            plt.plot(x,quotients,marker=type[j],markersize=7,linewidth=1,color="green",label=k.format(Title[j]))
    plt.legend(fontsize=12)   

    plt.xlabel('Temps',fontsize=10)
    plt.ylabel('Nombre des individues')
    plt.title("infected in Friendfeed")
    plt.grid(True)
    plt.savefig(title_fig+'infectedInFriendfeed.pdf',dpi=50)
    
    ###################################################################################
     #All Spreaders
    x = range(0,len(Stat[0]['Infected']))
    y= range(0,int(len(Stat[0]['Infected'])/3))
    x=np.array(x)*pro
    y=np.array(y)*pro
    xx+=1
    plt.figure(num=xx)
    plt.subplot()
    #k="{}:{},{}]" 
    k="{}" 
    for infected ,j in zip( Stat,range(len(Stat))):
      quotients =[]
      for i in range(0, len(infected["Spreaders"])-2, 3):
          q = (infected["Spreaders"][i] + infected["Spreaders"][i+1] + infected["Spreaders"][i+2] ) / Nodes
          quotients.append(q)
      plt.plot(y, quotients,marker=type[j],markersize=6,linewidth=1,label=k.format(Title[j]))
    
    plt.legend(fontsize=12)
    plt.grid(True)
    plt.title("Spreaders of the whole popularity")
    plt.xlabel('Temps')
    plt.ylabel('Nombre des individues')
    plt.savefig(title_fig+'SpreadersInWholePopularity.pdf',dpi=20)
    #Spreaders yt
    xx+=1
    plt.figure(num=xx)
    plt.subplot()
    #k="{}:{},{}]" 
    k="{}" 
    for Spreaders,j in zip( StatYt,range(len(StatYt))):
            quotients = [number /663  for number in Spreaders["Spreaders"]]
            plt.plot(x,quotients,marker=type[j],markersize=7,linewidth=1,color="red",label=k.format(Title[j]))
    plt.legend(fontsize=12)   

    plt.xlabel('Temps',fontsize=10)
    plt.ylabel('Nombre des individues')
    plt.title("Spreaders in youtube")
    plt.grid(True)
    plt.savefig(title_fig+'SpreadersYoutube.pdf',dpi=50)
    #Spreaders tw
    xx+=1
    plt.figure(num=xx)
    plt.subplot()
    #k="{}:{},{}]" 
    k="{}" 
    for Spreaders,j in zip( StatTw,range(len(StatTw))):
            quotients = [number /5540  for number in Spreaders["Spreaders"]]
            plt.plot(x,quotients,marker=type[j],markersize=7,linewidth=1,color="blue",label=k.format(Title[j]))
    plt.legend(fontsize=12)   

    plt.xlabel('Temps',fontsize=10)
    plt.ylabel('Nombre des individues')
    plt.title("Spreaders in twitter")
    plt.grid(True)
    plt.savefig(title_fig+'SpreadersTwitter.pdf',dpi=50)
    #Spreaders ff
    xx+=1
    plt.figure(num=xx)
    plt.subplot()
    #k="{}:{},{}]" 
    k="{}" 
    for Spreaders,j in zip( StatFf,range(len(StatFf))):
            quotients = [number /5702  for number in Spreaders["Spreaders"]]
            plt.plot(x,quotients,marker=type[j],markersize=7,linewidth=1,color="green",label=k.format(Title[j]))
    plt.legend(fontsize=12)   

    plt.xlabel('Temps',fontsize=10)
    plt.ylabel('Nombre des individues')
    plt.title("Spreaders in Friendfeed")
    plt.grid(True)
    plt.savefig(title_fig+'SpreadersInFriendfeed.pdf',dpi=50)
    ##################################################################################
    # Global RumorPopularity
    xx+=1
    plt.figure(num=xx)
    plt.subplot()
    #k="{}:{},{}]" 
    k="{}" 
    for infected,j in zip( Stat,range(len(Stat))):
      quotients = [number /Nodes  for number in infected["GlobalRumorPopularity"]]
      plt.plot(x, quotients,marker=type[j],markersize=6,linewidth=1,label=k.format(Title[j]))
    plt.legend(fontsize=12) 
    plt.xlabel('Temps')
    plt.ylabel('Nombre des individues')
    plt.grid(True)
    plt.title("Global Rumor Popularity")
    plt.savefig(title_fig+'RumorPopularityInWholePopularity.pdf',dpi=20)
    #RumorPopularity yt
    xx+=1
    plt.figure(num=xx)
    plt.subplot()
    #k="{}:{},{}]" 
    k="{}" 
    for RumorPopularity,j in zip( StatYt,range(len(StatYt))):
            quotients = [number /663  for number in RumorPopularity["RumorPopularity"]]
            plt.plot(x,quotients,marker=type[j],markersize=7,linewidth=1,color="red",label=k.format(Title[j]))
    plt.legend(fontsize=12)   

    plt.xlabel('Temps',fontsize=10)
    plt.ylabel('Nombre des individues')
    plt.title("RumorPopularity in youtube")
    plt.grid(True)
    plt.savefig(title_fig+'RumorPopularityYoutube.pdf',dpi=50)
    #RumorPopularity tw
    xx+=1
    plt.figure(num=xx)
    plt.subplot()
    #k="{}:{},{}]" 
    k="{}" 
    for RumorPopularity,j in zip( StatTw,range(len(StatTw))):
            quotients = [number /5540 for number in RumorPopularity["RumorPopularity"]]
            plt.plot(x,quotients,marker=type[j],markersize=7,linewidth=1,color="blue",label=k.format(Title[j]))
    plt.legend(fontsize=12)   

    plt.xlabel('Temps',fontsize=10)
    plt.ylabel('Nombre des individues')
    plt.title("RumorPopularity in twitter")
    plt.grid(True)
    plt.savefig(title_fig+'RumorPopularityTwitter.pdf',dpi=50)
    #RumorPopularity ff
    xx+=1
    plt.figure(num=xx)
    plt.subplot()
    #k="{}:{},{}]" 
    k="{}" 
    for RumorPopularity,j in zip( StatFf,range(len(StatFf))):
            quotients = [number /5702  for number in RumorPopularity["RumorPopularity"]]
            plt.plot(x,quotients,marker=type[j],markersize=7,linewidth=1,color="green",label=k.format(Title[j]))
    plt.legend(fontsize=12)   

    plt.xlabel('Temps',fontsize=10)
    plt.ylabel('Nombre des individues')
    plt.title("RumorPopularity in Friendfeed")
    plt.grid(True)
    plt.savefig(title_fig+'RumorPopularityInFriendfeed.pdf',dpi=50)
    ###################################################################################
    # Format the minor tick labels of the y-axis into empty strings with
    # `NullFormatter`, to avoid cumbering the axis with too many labels.
   
def Simulation(index,graph,Stat_Global,percentage):
     Beta=0.2
     with Manager() as manager:
        Stat=manager.list()  
        parameter=[]
        parameters(parameter,Beta=Beta+index/10)
        start_time = time.time()  
        processes=[multiprocessing.Process(target=Start,args=(i,index,graph,parameter,Stat,percentage))for i in range(10)] 
        [process.start() for process in processes] 
        [process.join() for process in processes]
        end_time = time.time() 
        print("Parallel xx time=", end_time - start_time)
        globalStat(Stat,Stat_Global,parameter)
#gene graph
def Random_networks ( N=300 ,P=0.3):
    # Erdős-Rényi graph
    # number of nodes
    # expected number of partners
    
    g = nx.gnp_random_graph(N, P)  
    return graphe_TO_json(g)
def Small_World_networks(N=300,K=10,P=0.3):
    
    #Watts_strogatz graph
    #N=number of nodes
    #K=Each node is joined with its k nearest neighbors in a ring topology(anneau).
    #P=The probability of rewiring each edge(Probabilite de remplace les arretes)
    g= nx.watts_strogatz_graph(N,K,P)
    return graphe_TO_json(g)
def Scale_free_networks (N=300,M=10):
    #Barabasi_albert graph
    #N= Number of nodes
    #M= Number of edges to attach from a new node to existing nodes
    g=nx.barabasi_albert_graph(N,M)
    return graphe_TO_json(g)
    
    return graphe_TO_json(g) 
def search_spreaders(G,sp):
    
    l=len(G.nodes)
    for i in range (l):
        if ( G.nodes[i]['state']=='spreaders'):
          sp.append(i)
                
def neighbor(Spreaders,g):
    neighb=[]
    MaxD=[]
    Cente=[]
    beta=[]
    betaD=[]
    Cent=((nx.degree_centrality(g)))
    
    for i in Spreaders:
        n=g.neighbors(i)
        
        for j in n:
          
          if g.nodes[j]['state'] =='non_infected':
              if j not in neighb :
                  neighb.append(j)
                  Cente.append(Cent[j])
                  MaxD.append(g.nodes[j]['degree'])
                  beta.append(g.nodes[j]['beta'])
                  betaD.append(g.nodes[j]['degree']/g.nodes[j]['beta'])
           
              
   
    return neighb,MaxD,Cente,beta,betaD
def simulation_strategy(x,K,timeOfTheRumorDetection,method,G,numberOfNetworks):
    allNodesInfected = set()
    GlobalRumorPopularity=0
    #for each network in the multiplex
        
    with Manager() as manager:
        Stat_Global=manager.list()
        GlobalStatYt=manager.list()
        GlobalStatTw=manager.list()
        GlobalStatFf=manager.list()
        v=0
        for met in method :
            print(met)
            with Manager() as manager:
                Stat=manager.list()  
                parameter=[]
                parameters(parameter)
                start_time = time.time()  
                processes=[multiprocessing.Process(target=Start,args=(i,G[i],allNodesInfected,GlobalRumorPopularity,parameter,Stat,percentage,K,timeOfTheRumorDetection,met))for i in range(len(G))] 
                
                [process.start() for process in processes] 
                [process.join() for process in processes]
                end_time = time.time() 
                print("Parallel xx time=", end_time - start_time)
                globalStat(Stat,Stat_Global,GlobalStatYt,GlobalStatTw,GlobalStatFf,parameter,met)
            v+=1     
        Display(Stat_Global,GlobalStatYt,GlobalStatTw,GlobalStatFf,x,'np',Nodes)   #NBLS
def Iterative():
    start_time = time.time()  
    StatI=[]

    for i in range(6):
        parameter=[]
        parameters(parameter,Omega=0.2+i/10)
        Stat=[]
        start_time1 = time.time() 
        for j in range(50):
            Start(i,j,g,parameter,Stat,percentage)
        end_time1 = time.time()
        print("Serial xx time=", end_time1 - start_time1)
        globalStat(Stat,StatI,parameter)
    end_time = time.time()
    print("Serial time=", end_time - start_time)
    Display(StatI)
def Random_Blocking_nodes(Graphe,k):
    sp=[]
    search_spreaders(Graphe,sp)
    nb,d,cen,Bet,betaD=neighbor(sp,Graphe)
    size=len(nb)
    if k>size:
      k=size-1
    for i in range(k):
        s=random.randint(0, size-1)
        Graphe.nodes[nb[s]]['blocked']='True'
        nb.pop(s)
        size-=1
       
def Degree_MAX_Blocking_nodes(G,k):
    
    sp=[]
   
    search_spreaders(G,sp)
   
    nb,DNode,cen,Bet,betaD=neighbor(sp,G)
    

    for i in range(k):
            
            ID = DNode.index(max(DNode))
            G.nodes[nb[ID]]['blocked']='True'
            DNode.pop(ID)
            nb.pop(ID)
            
def Centrality_Blocking_nodes(G,k):
    
    sp=[]
   
    search_spreaders(G,sp)
    nb,DNode,cen,Bet,betaD=neighbor(sp,G)
    for i in range(k):
            
            ID = cen.index(max(cen))
            G.nodes[nb[ID]]['blocked']='True'
            cen.pop(ID)
            nb.pop(ID)         
def Beta_Blocking_nodes(G,k):
    
    sp=[]
   
    search_spreaders(G,sp)
   
    nb,DNode,cen,Bet,betaD=neighbor(sp,G)
    

    for i in range(k):
            
            ID = Bet.index(min(Bet))
            G.nodes[nb[ID]]['blocked']='True'
            Bet.pop(ID)
            nb.pop(ID)         
def BetaD_Blocking_nodes(G,k):
    
    sp=[]
   
    search_spreaders(G,sp)
   
    nb,DNode,cen,Bet,betaD=neighbor(sp,G)
    

    for i in range(k):
            
            ID = betaD.index(max(betaD))
            G.nodes[nb[ID]]['blocked']='True'
            betaD.pop(ID)
            nb.pop(ID)   

def Random_TRuth_comp(Graphe,k):
    sp=[]
    search_spreaders(Graphe,sp)
    nb,d,cen,Bet,betaD=neighbor(sp,Graphe)
    size=len(nb)
    if k > size :
       k=size-1
    for i in range(k):
        s=random.randint(0, size-1)
        Graphe.nodes[nb[s]]['Protector']='True'
        Graphe.nodes[nb[s]]['state']='infected'
        nb.pop(s)
        size-=1
def MaxDegree_TRuth_comp(Graphe,K):
    sp=[]
    search_spreaders(Graphe,sp)
    nb,d,cen,Bet,betaD=neighbor(sp,Graphe)
    size=len(nb)
    k=K
    if k > size :
       k=size-1
    for i in range(k):
        s = d.index(max(d))
        Graphe.nodes[nb[s]]['Protector']='True'
        Graphe.nodes[nb[s]]['state']='infected'
        nb.pop(s)
        d.pop(s)
def Centrality_TRuth_comp(Graphe,K):
    sp=[]
    search_spreaders(Graphe,sp)
   
    nb,d,cen,Bet,betaD=neighbor(sp,Graphe)
   
    size=len(nb)
    k=K
    if k > size :
       k=size-1
    for i in range(k):
        s = cen.index(max(cen))
        Graphe.nodes[nb[s]]['Protector']='True'
        Graphe.nodes[nb[s]]['state']='infected'
        nb.pop(s)
        cen.pop(s)
def Beta_TRuth_comp(Graphe,K):
    sp=[]
    search_spreaders(Graphe,sp)
    nb,d,bet,cen,betaD=neighbor(sp,Graphe)
    size=len(nb)
    k=K
    if k > size :
       k=size-1
    for i in range(k):
        s = cen.index(min(cen))
        Graphe.nodes[nb[s]]['Protector']='True'
        Graphe.nodes[nb[s]]['state']='infected'
        nb.pop(s)
        cen.pop(s)
def BetaD_TRuth_comp(Graphe,K):
    sp=[]
    search_spreaders(Graphe,sp)
    nb,d,bet,cen,betaD=neighbor(sp,Graphe)
    size=len(nb)
    k=K
    if k > size :
       k=size-1
    for i in range(k):
        s = betaD.index(max(betaD))
        Graphe.nodes[nb[s]]['Protector']='True'
        Graphe.nodes[nb[s]]['state']='infected'
        nb.pop(s)
        betaD.pop(s)


def blocked(G):
    
    L=[]
    for i in range (len(G.nodes)):
        if(G.nodes[i]['blocked']=='True'):
            L.append(i)           
    return L
def Protector(G):  
    
    L=[]
    for i in range (len(G.nodes)):
        if(G.nodes[i]['Protector']=='True'):
            L.append(i)           
    return L
if __name__ == '__main__':
       # use net.Graph() for undirected graph

# How to read from a file. Note: if your egde weights are int, 
# change float to int.
   
    #Graph's Parametres 
    P=0.3
    K=100
    M=7
    nbb=0
   
    #g=json_graph.node_link_graph(Scale_free_networks(Nodes,M))
    #g=Scale_free_networks(Nodes,M)
    #g=json_graph.node_link_graph(Small_World_networks(Nodes,K,P))
    #g=json_graph.node_link_graph(Random_networks(Nodes,P))
    ytGraph=json_graph.node_link_graph(geneGraphs.txt2Graph("yt.txt"))
    twGraph=json_graph.node_link_graph(geneGraphs.txt2Graph("tw.txt"))
    ffGraph=json_graph.node_link_graph(geneGraphs.txt2Graph("ff.txt"))  

    G={"yt":ytGraph, "tw":twGraph, "ff":ffGraph}
    numberOfNetworks=3
    #G={"tw":twGraph} 
    NumOFsumi=10
    G=[]
    for i in range(NumOFsumi):
      G.append(ytGraph)
      G.append(twGraph)
      G.append(ffGraph)
    Nodes=6407
    allNodesInfected={}
    static="Nodes :{},Edegs:{}."
    percentage=5 #1% of popularity" is infected 
    
    beta=0.2
    omega=0
    juge=0.1
    delta=0
    K=int(Nodes*0.2)
    print(K)
    timeOfTheRumorDetection=1
    
    
  

    
    simulation_strategy(1,  K, timeOfTheRumorDetection, ["NP"],G,numberOfNetworks)
    plt.show()
