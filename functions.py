#%%
import matplotlib.pyplot as plt
import numpy as np
from math import factorial
from pynverse import inversefunc
import numpy as np
from scipy.optimize import fsolve
from collections import Counter
import random as rand
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

from scipy.stats import *

#%%
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:.2f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()
def dyn_vax(beta,gamma,omega,graph,radius,tmax,block = True):
    #Note that with omega = 0, the model is the SIR model
    S = []
    I = []
    R = []
    V = []
    S.append(len(graph.nodes)-1)
    I.append(1)
    R.append(0)
    V.append(0)
    list_nodes = list(graph.nodes)
    if len(list_nodes) == 0:
        return S,I,R,V
    #choose a random node to be infected
    patient0 = rand.choice(list_nodes)
    infected = [patient0]
    #choose a random node to be recovered
    recovered = []
    vaccinated = []
    #choose a random node to be susceptible
    susceptible = [i for i in list_nodes if i not in infected]
    t=0
    while t < tmax:
        #print("time: "+str(t))
        #add new infected nodes
        inf1 = infected.copy() #to avoid changing the list while iterating
        #shuffle the list of infected nodes 
        rand.shuffle(inf1)
        newly_infected = []
        for node in inf1:
            tot_neigh = nx.single_source_shortest_path_length(graph,node,radius)
            neigh = set(key for key, value in tot_neigh.items() if value <= radius) if block else set(key for key, value in tot_neigh.items() if value == radius)

            #remove node from neigh
            #neigh.remove(node)
            for neighbor in neigh:
                if neighbor in susceptible:
                    k = rand.random()
                    if k < omega:
                        #print("k: "+str(k)+": new vaccinated node: "+str(neighbor))
                        vaccinated.append(neighbor)
                        susceptible.remove(neighbor)
            for neighbor in graph.neighbors(node):
                if neighbor in susceptible:
                    k = rand.random()
                    if k < beta:
                        #print("k: "+str(k)+": new infected node: "+str(neighbor))
                        infected.append(neighbor)
                        newly_infected.append(neighbor)
                        susceptible.remove(neighbor)
        for node in infected:
            if node not in newly_infected:
            #print("node: "+str(node))
                if rand.random() < gamma:
                    #print("node: "+str(node)+" is recovered")
                    recovered.append(node)
                    infected.remove(node)
                #print("dict_recovering_nodes[node]: "+str(dict_recovering_nodes[node])
        #update lists
        #print("t",t)
        t += 1
        #print("S: "+str(len(susceptible)))
        S.append(len(susceptible))
        I.append(len(infected))
        R.append(len(recovered))
        V.append(len(vaccinated))
        #print("Susceptible: "+str(len(susceptible)))
        #print("Infected: "+str(len(infected)))
       # print("Recovered: "+str(len(recovered)))
        #print sum
        #print("tot",sum([len(susceptible),len(infected),len(recovered),len(vaccinated)]))
    return S,I,R,V

#%%
def get_k(graph):
    list_degrees = []
    for node in graph.nodes:
        list_degrees.append(graph.degree(node))
    k = sum(list_degrees)/len(graph.nodes)
    return k

def get_k_s(graph):
    list_degrees = []
    for node in graph.nodes:
        list_degrees.append(graph.degree(node))
    k_s = sum([el**(2) for el in list_degrees])/len(graph.nodes)
    return k_s


def g_1(x,prob_list,mean_degree):
    sol = 0
    for i in range(0,len(prob_list)-1):
        if prob_list[i] == 0:
            continue
        else:
            sol += ((i+1)*prob_list[i+1]*x**i)/mean_degree
    return sol

def g_0(x,prob_list):
    sol = 0
    for i in range(len(prob_list)-1):
        sol += (prob_list[i]*(x**(i)))
    return sol

def newton_method(f,df,x0 = 0.5, epsilon=0.0001):
    x = x0
    while abs(f(x)) > epsilon:
        x = x - f(x)/df(x)
    return x


def factorial(n):
    if n == 0:
        return 1
    else:
        return n*factorial(n-1)

#%%
#### get probability distributions
def prob_list_poisson(lambda1,n):
    prob_vec = []
    k = 0
    while k<n:
        prob_vec.append((poisson.pmf(k,lambda1)))
        k+=1
    prob_vec = [i/sum(prob_vec) for i in prob_vec]
    return prob_vec

def prob_list_powerlaw(a,min,n):
    prob_vec = [0]*min
    for k in range(min,int(np.sqrt(n))):
        prob_vec.append(k**(-a))
    prob_vec = [i/sum(prob_vec) for i in prob_vec]
    return prob_vec




### get degree lists

def discrete_inverse_trans(prob_vec):
    U=uniform.rvs(size=1)
    if U<=prob_vec[0]:
        return 1
    else:
        for i in range(1,len(prob_vec)+1):
            if sum(prob_vec[0:i])<U and sum(prob_vec[0:i+1])>U:
                return i

def discrete_samples(prob_vec,n):
    sample=[]
    for i in range(0,n):
        k = discrete_inverse_trans(prob_vec)
        if k == None:
            k = 1
        sample.append(k)

    for i in range(len(sample)):
        if sample[i] == None:
            print(sample[i],i)

    #check if sum is even   
    if sum(sample)%2 != 0:
        sample[0] = sample[0]+1
    return sample

##### get configuration model 

def config_model(degree_list):
    #check if sum is even
    if sum(degree_list)%2 != 0:
        degree_list[0] = degree_list[0]+1
    graph = nx.configuration_model(degree_list)
    graph1= graph.copy()
    for (u,v) in graph1.edges():
        if u == v:
            graph.remove_edge(u,v)
    #remove parallel edges
    graph = nx.Graph(graph)
    return graph
    

##### percolation

def remove_edges_uniform(graph,p):
    lista_random = []
    for i in range(len(graph.edges)):
        lista_random.append(rand.random())
    graph1 = graph.copy()
    i = 0
    for edge in graph.edges:
        if lista_random[i] > p:
            graph1.remove_edge(edge[0],edge[1])
        i = i+1
    return graph1


def remove_nodes_degree(graph,p):
    graph1 = graph.copy()
    i = 0
    for node in graph.nodes:
        if graph.degree(node) <p:
            graph1.remove_node(node)
        i = i+1
    return graph1

def remove_nodes_uniform(graph,p):
    lista_random = []
    for i in range(len(graph.nodes)):
        lista_random.append(rand.random())
    graph1 = graph.copy()
    i = 0
    for node in graph.nodes:
        if lista_random[i] > p:
            graph1.remove_node(node)
        i = i+1
    return graph1

##### Get size of biggest component after percolation 

def get_size_biggest_component(g):
    if len(g.nodes) == 0:
        return 0
    else:
        return len(sorted(nx.connected_components(g), key = len, reverse=True)[0])
    
def f_0(z,prob_vec,max_phi):
    sol = 0
    for k in range(max_phi):
        sol = sol + prob_vec[k]*(z**k)
    return sol

def f_1(z,prob_vec,max_phi,mean_deg):
    sol = 0
    for k in range(max_phi):
        if k != 0:
            sol = sol + k*prob_vec[k]*(z**(k-1))
    sol = sol/mean_deg
    return sol

def get_S_synth_site(degree_list, n):
    size_biggest_component = [0]*n
    p = np.linspace(0,1,n)
    k = 0
    k_s = 0
    for j in range(n):
        graph = config_model(degree_list)
        for i in range(len(p)):
            graph1 = remove_nodes_uniform(graph,p[i])
            size_biggest_component[i] = size_biggest_component[i] + get_size_biggest_component(graph1)
        k +=get_k(graph1)
        k_s += get_k_s(graph1)
        print(j)
    size_biggest_component = [elem/(n*len(degree_list)) for elem in size_biggest_component]
    k = k/n
    k_s = k_s/n
    return size_biggest_component, k, k_s

def get_S_theory_site(prob_vec):
    p = np.linspace(0,1,10)
    degree_list = discrete_samples(prob_vec,1000)
    k1 = sum([degree_list[i] for i in range(len(degree_list))])/len(degree_list)
    k_s1 = sum([degree_list[i]**(2) for i in range(len(degree_list))])/len(degree_list)
    S = []
    for i in range(len(p)):
        x_intersect = sp.optimize.fsolve(lambda x: 1-p[i]+p[i]*g_1(x,prob_vec,np.mean(degree_list)) - x,0.05)
        S.append(p[i]*(1-g_0(x_intersect,prob_vec)))
    return S, k1, k_s1


def get_S_synth_bond(degree_list, n):
    size_biggest_component = [0]*n
    p = np.linspace(0,1,n)
    k = 0
    k_s = 0
    graph = config_model(degree_list)
    for j in range(n):
        for i in range(len(p)):
            graph1 = remove_edges_uniform(graph,p[i])
            size_biggest_component[i] = size_biggest_component[i] + get_size_biggest_component(graph1)
        k +=get_k(graph1)
        k_s += get_k_s(graph1)
        print(j)
    size_biggest_component = [elem/(n*len(degree_list)) for elem in size_biggest_component]
    k = k/n
    k_s = k_s/n
    return size_biggest_component, k, k_s

def get_S_theory_bond(prob_vec):
    p = np.linspace(0,1,100)
    degree_list = discrete_samples(prob_vec,len(prob_vec))
    k1 = sum([degree_list[i] for i in range(len(degree_list))])/len(degree_list)
    k_s1 = sum([degree_list[i]**(2) for i in range(len(degree_list))])/len(degree_list)
    S = []
    for i in range(len(p)):
        x_intersect = sp.optimize.fsolve(lambda x: 1-p[i]+p[i]*g_1(x,prob_vec,np.mean(degree_list)) - x,0.05)
        S.append(abs((1-g_0(x_intersect,prob_vec))))
    return S, k1, k_s1


def get_S_theory_non_uniform(degree_list,prob_vec,mean_deg):
    S = []
    for max_phi in range(max(degree_list)):
        x_intersect = sp.optimize.fsolve(lambda x: 1-f_1(1,prob_vec,max_phi,mean_deg)+f_1(x,prob_vec,max_phi,mean_deg) - x,0.05)
        S.append(f_0(1,prob_vec,max_phi)-f_0(x_intersect,prob_vec,max_phi))
    return S

def get_S_synth_non_uniform(graph):
    components = []
    for max_phi in range(max([deg[1] for deg in graph.degree])):
        graph1 = graph.copy()
        for node in graph.nodes:
            if graph.degree[node] >= max_phi:
                graph1.remove_node(node)
        components.append((get_size_biggest_component(graph1))/len(graph.nodes))
    return components
# %%
def prob_list(degree_list,n): #get a probability distribution from a given degree list
    prob_vec = []
    for i in range(n):
        prob_vec.append(degree_list.count(i)/len(degree_list))
    return prob_vec
# %%
