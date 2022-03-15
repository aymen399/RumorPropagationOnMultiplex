from functools import partial
from itertools import repeat
from multiprocessing import Pool, freeze_support

def func(a, b):
    print( a + b)

def main():
    a_args = [1,2,3]
    second_arg = 1
    with Pool() as pool:
        L = pool.starmap(func, [(1, 1), (2, 1), (3, 1)])
       
        #print(N)

if __name__=="__main__":
    freeze_support()
    main()

#todo: parallelisme of 3 simulations with diffent argument(the graph) use ....) for i,g 
''' 
dict1 = {"tryIndice": "Bitcoin", 2: "Ethereum"} contient l'indice i et le graphe g de les trois different propagations;ff,tw et yt
for key, value in dict1.items():
    print(f"Key {key} has value {value}")
'''

#Todo: then parallize the afromentioned si mafter putting it in a function as start() let's be start2()
#Todo: then add the pk of the layer of the multiplex and add the list of all nodes in dict to indicate if they are infected or not
# to update the llst of infected nodes in each network if the infected nodes update is coming from another network: xxx
# add xxx in model2.py before line 150.