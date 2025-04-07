
### DYNAMIC BLOCK AND RING VACCINATION (THEORY AND SIMULATION) ###


#%%
from functions import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#%%
#define a progressbar function to monitor the progress of the simulation
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


## FUNCTION FOR SIMULATION OF DYNAMIC VACCINATION 
def dyn_vax(beta,gamma,omega,graph,radius,tmax,block = True):
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
    recovered = []
    vaccinated = []
    susceptible = [i for i in list_nodes if i not in infected]
    t=0
    while t < tmax:
        inf1 = infected.copy() #to avoid changing the list while iterating
        #shuffle the list of infected nodes 
        rand.shuffle(inf1)
        newly_infected = []
        for node in inf1:
            tot_neigh = nx.single_source_shortest_path_length(graph,node,radius)
            neigh = set(key for key, value in tot_neigh.items() if value <= radius) if block else set(key for key, value in tot_neigh.items() if value == radius)
            for neighbor in neigh:
                if neighbor in susceptible:
                    k = rand.random()
                    if k < omega:
                        vaccinated.append(neighbor)
                        susceptible.remove(neighbor)
            for neighbor in graph.neighbors(node):
                if neighbor in susceptible:
                    k = rand.random()
                    if k < beta:
                        infected.append(neighbor)
                        newly_infected.append(neighbor)
                        susceptible.remove(neighbor)
        for node in infected:
            if node not in newly_infected:
                if rand.random() < gamma:
                    recovered.append(node)
                    infected.remove(node)
        t += 1
        S.append(len(susceptible))
        I.append(len(infected))
        R.append(len(recovered))
        V.append(len(vaccinated))
    return S,I,R,V


## FUNCTION FOR SIMULATION OF DYNAMIC VACCINATION (VARIATION WITH VACCINE EFFICACY <100%)
def dyn_vax_lambda(beta,gamma,omega,graph,radius,tmax,block = True):
    S = []
    I = []
    R = []
    V = []
    lambdino = 0.5
    S.append(len(graph.nodes)-1)
    I.append(1)
    R.append(0)
    V.append(0)
    list_nodes = list(graph.nodes)
    if len(list_nodes) == 0:
        return S,I,R,V
    patient0 = rand.choice(list_nodes)
    infected = [patient0]
    recovered = []
    vaccinated = []
    imposters = []
    imposters_infected = []
    newly_imposters_infected = []
    susceptible = [i for i in list_nodes if i not in infected]
    t=0
    while t < tmax:
        inf1 = infected.copy()
        rand.shuffle(inf1)
        newly_infected = []
        for node in inf1:
            tot_neigh = nx.single_source_shortest_path_length(graph,node,radius)
            neigh = set(key for key, value in tot_neigh.items() if value <= radius) if block else set(key for key, value in tot_neigh.items() if value == radius)

            for neighbor in neigh:
                if neighbor in susceptible:
                    k = rand.random()
                    if k < omega:
                        l = rand.random()
                        if l<lambdino:
                            vaccinated.append(neighbor)
                            susceptible.remove(neighbor)
                        else:
                            imposters.append(neighbor)
                            susceptible.remove(neighbor)
                            
            for neighbor in graph.neighbors(node):
                if neighbor in susceptible:
                    k = rand.random()
                    if k < beta:
                        infected.append(neighbor)
                        newly_infected.append(neighbor)
                        susceptible.remove(neighbor)
                elif neighbor in imposters:
                    k = rand.random()
                    if k < beta:
                        imposters_infected.append(neighbor)
                        imposters.remove(neighbor)
                        
        imp_inf1 = imposters_infected.copy()
        for node in imp_inf1:
            tot_neigh = nx.single_source_shortest_path_length(graph,node,radius)
            neigh = set(key for key, value in tot_neigh.items() if value <= radius) if block else set(key for key, value in tot_neigh.items() if value == radius)
            for neighbor in graph.neighbors(node):
                if neighbor in susceptible:
                    k = rand.random()
                    if k < beta:
                        infected.append(neighbor)
                        newly_infected.append(neighbor)
                        susceptible.remove(neighbor)
                elif neighbor in imposters:
                    
                    k = rand.random()
                    if k < beta:
                        imposters_infected.append(neighbor)
                        imposters.remove(neighbor)
        
        for node in infected:
            if node not in newly_infected:
                if rand.random() < gamma:
                    recovered.append(node)
                    infected.remove(node)

        t += 1
        S.append(len(susceptible))
        I.append(len(infected))
        R.append(len(recovered))
        V.append(len(vaccinated)+len(imposters)+len(imposters_infected))

    return S,I,R,V

#%%#Initialization of the parameters: Simulation of dynamic Vaccination
Block = True
path_sim = ""
prob_vec = prob_list_powerlaw(2.5,3,400)
#prob_vec = prob_list_poisson(4,20)
n_nodes = #insert number of nodes in the initial graph
niter = #number of stochastic simulations
Beta = np.linspace(0,1,20) #set 
gamma = 1/3
Omega = np.linspace(0,1,20)
radius = #choose vaccination radius
t_max = 50

#%%
#SIMULATION
df = pd.DataFrame(columns = ["Omega", "Beta", "Recovered", "Infected", "Susceptible", "Vaccinated"])
#df = pd.read_csv("/Users/matteoscianna/Desktop/Paper Tesi/data/dynamic_vax/simulation/final/dynamic_vax_gamma_1-3.csv")
for omega in Omega:
    #if omega == Omega[-1]:
    #    Beta = np.linspace(0,1,20)[:]
    #else:
    #    Beta = np.linspace(0,1,20)
    for beta in Beta:
        iterat = 0
        S = [0]*t_max
        I = [0]*t_max
        R = [0]*t_max
        V = [0]*t_max
        while iterat < niter:
            printProgressBar(iterat,niter,prefix = "Progress",suffix = "Complete",length = 50)
            degree_vec_in = discrete_samples(prob_vec,n_nodes)
            G = config_model(degree_vec_in)
            sir = dyn_vax(beta,gamma,omega,G,radius,t_max,block = Block) #block = True for block vaccination
            S = [S[i]+sir[0][i] for i in range(0,len(S))]
            I = [I[i]+sir[1][i] for i in range(0,len(I))]
            R = [R[i]+sir[2][i] for i in range(0,len(R))]
            V = [V[i]+sir[3][i] for i in range(0,len(V))]
            iterat += 1
        
        #average
        S = [i/(niter*n_nodes) for i in S]
        I = [i/(niter*n_nodes) for i in I]
        R = [i/(niter*n_nodes) for i in R]
        V = [i/(niter*n_nodes) for i in V]
        print("Results for omega :"+str(omega)+", beta: "+str(beta)+ "::: +\n"+ "Recovered: "+str(R[-1])+"\n"+"Infected: "+str(I[-1])+"\n"+"Susceptible: "+str(S[-1])+"\n"+"Vaccinated: "+str(V[-1]))
        df.loc[len(df)] = [omega,beta,R[-1],I[-1],S[-1],V[-1]]
#%%
###############  THEORY  #########################
#try with vax - r=2
path_th = ""
prob_vec = prob_list_powerlaw(2.5,3,400)
Block = True
n_nodes = #insert number of nodes in the initial graph
niter = #number of stochastic simulations
Beta = np.linspace(0,1,20) #set 
gamma = 1/3
Omega = np.linspace(0,1,20)
radius = #choose vaccination radius
t_max = 50
#%%
for omega in Omega:
    for beta in Beta:
        #theory
        MEAN_S = np.zeros(51)
        MEAN_S[0]  = MEAN_S[0]  = ((n_nodes-1)/n_nodes)*niter
        MEAN_I = np.zeros(51)
        MEAN_R = np.zeros(51)
        MEAN_V = np.zeros(51)
        S = [0]*t_max
        I = [0]*t_max
        R = [0]*t_max
        V = [0]*t_max
        iterat = 0
        while iterat < niter:
            printProgressBar(iterat,niter,prefix = "Progress",suffix = "Complete",length = 50)
            #initialization of the graph
            degree_vec_in = discrete_samples(prob_vec,n_nodes)
            G = config_model(degree_vec_in)
            ADJ_DICT_TRUE = {}
            for i in range(len(G.nodes)):
                tot_neigh = nx.single_source_shortest_path_length(G,i,radius)
                nth_neigh = set(key for key, value in tot_neigh.items() if value > 0 and value <= radius) if Block else set(key for key, value in tot_neigh.items() if value == radius)
                ADJ_DICT_TRUE[i] = {}
                ADJ_DICT_TRUE[i]["first_neigh"] = {neigh:{"theta-t-1":1,"theta-t":1,"phi-t-1":1/n_nodes,"phi-t":1/n_nodes,"delta-t-1":1,"delta-t":1} for neigh in G.neighbors(i)}
                ADJ_DICT_TRUE[i]["nth_neigh"] = {neigh:{"delta-t-1":1,"delta-t":1,"phi-t-1":1/n_nodes,"phi-t":1/n_nodes} for neigh in nth_neigh}
            t = 0
            MAT_S = np.zeros((n_nodes,t_max+1))
            MAT_S[:,0] = (n_nodes-1)/n_nodes
            MAT_R = np.zeros((n_nodes,t_max+1))
            MAT_I = np.zeros((n_nodes,t_max+1))
            MAT_I[:,0] = 1/n_nodes
            MAT_V = np.zeros((n_nodes,t_max+1))
            ADJ_DICT = ADJ_DICT_TRUE.copy()

            while t < t_max:
                printProgressBar(t,t_max,prefix = "Time",suffix = "Complete",length = 50)
                t+=1
                for i in ADJ_DICT.keys():
                    for k in ADJ_DICT[i]["first_neigh"].keys():
                        ADJ_DICT[i]["first_neigh"][k]["theta-t-1"] = ADJ_DICT[i]["first_neigh"][k]["theta-t"]    
                        ADJ_DICT[i]["first_neigh"][k]["theta-t"] -= beta*(1-omega)*ADJ_DICT[i]["first_neigh"][k]["phi-t"]
                        ADJ_DICT[i]["first_neigh"][k]["delta-t-1"] = ADJ_DICT[i]["first_neigh"][k]["delta-t"]
                        ADJ_DICT[i]["first_neigh"][k]["delta-t"] -= omega*ADJ_DICT[i]["first_neigh"][k]["phi-t"]
                    
                    for k in ADJ_DICT[i]["nth_neigh"].keys():
                        ADJ_DICT[i]["nth_neigh"][k]["delta-t-1"] = ADJ_DICT[i]["nth_neigh"][k]["delta-t"]
                        ADJ_DICT[i]["nth_neigh"][k]["delta-t"] -= omega*ADJ_DICT[i]["nth_neigh"][k]["phi-t"]

                for i in ADJ_DICT.keys():
                    for k in ADJ_DICT[i]["first_neigh"].keys():        
                        save= ADJ_DICT[i]["first_neigh"][k]["phi-t"]
                        ADJ_DICT[i]["first_neigh"][k]["phi-t"] = ADJ_DICT[i]["first_neigh"][k]["phi-t-1"]*(1-gamma)*(1-beta)*(1-omega) + (MAT_S[i,0])*((np.prod([ADJ_DICT[k]["first_neigh"][j]["theta-t-1"] for j in ADJ_DICT[k]["first_neigh"].keys() if j!= i])-np.prod([ADJ_DICT[k]["first_neigh"][j]["theta-t"] for j in ADJ_DICT[k]["first_neigh"].keys() if j!= i]))*(np.prod([ADJ_DICT[k]["first_neigh"][j]["delta-t"] for j in ADJ_DICT[k]["first_neigh"].keys() if j!= i])))*(np.prod([ADJ_DICT[k]["nth_neigh"][j]["delta-t"] for j in ADJ_DICT[k]["nth_neigh"].keys() if j!= i]))
                        ADJ_DICT[i]["first_neigh"][k]["phi-t-1"] = save

                for i in ADJ_DICT.keys():
                    for k in ADJ_DICT[i]["nth_neigh"].keys():        
                        save= ADJ_DICT[i]["nth_neigh"][k]["phi-t"]
                        ADJ_DICT[i]["nth_neigh"][k]["phi-t"] = ADJ_DICT[i]["nth_neigh"][k]["phi-t-1"]*(1-gamma)*(1-beta)*(1-omega) + (MAT_S[i,0])*(np.prod([ADJ_DICT[k]["first_neigh"][j]["theta-t-1"] for j in ADJ_DICT[k]["first_neigh"].keys() if j!= i])-np.prod([ADJ_DICT[k]["first_neigh"][j]["theta-t"] for j in ADJ_DICT[k]["first_neigh"].keys() if j!= i]))*(np.prod([ADJ_DICT[k]["nth_neigh"][j]["delta-t"] for j in ADJ_DICT[k]["nth_neigh"].keys() if j!= i]))
                        ADJ_DICT[i]["nth_neigh"][k]["phi-t-1"] = save

                for i in ADJ_DICT.keys():
                        MAT_S[i,t] = MAT_S[i,0]*(np.prod([ADJ_DICT[i]["first_neigh"][k]["theta-t"] for k in ADJ_DICT[i]["first_neigh"].keys()]))*(np.prod([ADJ_DICT[i]["nth_neigh"][k]["delta-t"] for k in ADJ_DICT[i]["nth_neigh"].keys()]))
                        MAT_R[i,t] = MAT_R[i,t-1] + gamma*MAT_I[i,t-1]
                        #MAT_V[i,t] = probability V at t-1 + probability S at t-1 * probability V at t
                        MAT_I[i,t] = MAT_I[i,t-1]*(1-gamma) + (MAT_S[i,0])*(np.prod([ADJ_DICT[i]["first_neigh"][j]["theta-t-1"] for j in ADJ_DICT[i]["first_neigh"].keys() if j!= i])-np.prod([ADJ_DICT[i]["first_neigh"][j]["theta-t"] for j in ADJ_DICT[i]["first_neigh"].keys() if j!= i]))*(np.prod([ADJ_DICT[i]["nth_neigh"][k]["delta-t"] for k in ADJ_DICT[i]["nth_neigh"].keys()]))
                        MAT_V[i,t] = 1-MAT_S[i,t]-MAT_I[i,t]-MAT_R[i,t]

                
            
                #print("t: "+str(t))
                MEAN_S[t] += np.mean(MAT_S[:,t])
                MEAN_I[t] += np.mean(MAT_I[:,t])
                MEAN_R[t] += np.mean(MAT_R[:,t])
                MEAN_V[t] += np.mean(MAT_V[:,t])
            iterat += 1
        print("Results for omega :"+str(omega)+", beta: "+str(beta)+ "::: +\n"+ "Recovered: "+str(MEAN_R[-1])+"\n"+"Infected: "+str(MEAN_I[-1])+"\n"+"Susceptible: "+str(MEAN_S[-1])+"\n"+"Vaccinated: "+str(MEAN_V[-1]))

    

#%%
Block = True
beta = 1.0
omega = 0.15
gamma = 1
radius = 2
#Here we plot the evolution of the SIR-V model with dynamic vaccination
if Block:
    R = np.load(path_sim+"/block_vax/SF_2.5/EVOL_SIR-V/R_beta="+str(beta)+",omega="+str(omega)+",gamma="+str(gamma)+",r="+str(radius)+".npy")
    I = np.load(path_sim+"/block_vax/SF_2.5/EVOL_SIR-V/I_beta="+str(beta)+",omega="+str(omega)+",gamma="+str(gamma)+",r="+str(radius)+".npy")
    S = np.load(path_sim+"/block_vax/SF_2.5/EVOL_SIR-V/S_beta="+str(beta)+",omega="+str(omega)+",gamma="+str(gamma)+",r="+str(radius)+".npy")
    V = np.load(path_sim+"/block_vax/SF_2.5/EVOL_SIR-V/V_beta="+str(beta)+",omega="+str(omega)+",gamma="+str(gamma)+",r="+str(radius)+".npy")
else:

    R = np.load(path_sim+"/ring_vax/SF_2.5/EVOL_SIR-V/R_beta="+str(beta)+",omega="+str(omega)+",gamma="+str(gamma)+",r="+str(radius)+".npy")
    I = np.load(path_sim+"/ring_vax/SF_2.5/EVOL_SIR-V/I_beta="+str(beta)+",omega="+str(omega)+",gamma="+str(gamma)+",r="+str(radius)+".npy")
    S = np.load(path_sim+"/ring_vax/SF_2.5/EVOL_SIR-V/S_beta="+str(beta)+",omega="+str(omega)+",gamma="+str(gamma)+",r="+str(radius)+".npy")
    V = np.load(path_sim+"/ring_vax/SF_2.5/EVOL_SIR-V/V_beta="+str(beta)+",omega="+str(omega)+",gamma="+str(gamma)+",r="+str(radius)+".npy")

t = np.arange(0,50,1)
plt.xticks(fontsize = 0)
plt.yticks(fontsize = 0)
plt.plot(t,R,label = "R",color = "grey",linestyle = "--")
plt.plot(t,I,label = "I",color = "red",linestyle = "--")
plt.plot(t,S,label = "S",color = "green",linestyle = "--")
plt.plot(t,V,label = "V",color = "blue",linestyle = "--")
plt.scatter(t,R,edgecolor = "grey",color = "w")
plt.scatter(t,I,edgecolor = "red",color = "w")
plt.scatter(t,S,edgecolor = "green",color = "w")
plt.scatter(t,V,edgecolor = "blue",color = "w")
#legend on two columns
plt.legend(loc = "upper right",fontsize = 20,ncol = 2)
plt.xlabel(r"$t$",fontsize = 25)
plt.ylabel("Fraction of nodes",fontsize = 25)

# %%
Ring = True
omega = 0.05
Beta = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0]
evol_r1 = []
evol_r2 = []
evol_r3 = []
for beta in Beta:
    if Ring:
        evol_r1_th.append(np.load(path_sim+"/ring_vax/SF_2.5/EVOL_SIR-V/V_beta="+str(beta)+",omega="+str(omega)+",gamma="+str(gamma)+",r=1.npy")[-1])
        evol_r2_th.append(np.load(path_sim+"/ring_vax/SF_2.5/EVOL_SIR-V/V_beta="+str(beta)+",omega="+str(omega)+",gamma="+str(gamma)+",r=2.npy")[-1])
        evol_r3_th.append(np.load(path_sim+"/ring_vax/SF_2.5/EVOL_SIR-V/V_beta="+str(beta)+",omega="+str(omega)+",gamma="+str(gamma)+",r=3.npy")[-1])
    else:
        evol_r1_sim.append(np.load(path_sim+"/block_vax/SF_2.5/EVOL_SIR-V/V_beta="+str(beta)+",omega="+str(omega)+",gamma="+str(gamma)+",r=1.npy")[-1])
        evol_r2_sim.append(np.load(path_sim+"/block_vax/SF_2.5/EVOL_SIR-V/V_beta="+str(beta)+",omega="+str(omega)+",gamma="+str(gamma)+",r=2.npy")[-1])
        evol_r3_sim.append(np.load(path_sim+"/block_vax/SF_2.5/EVOL_SIR-V/V_beta="+str(beta)+",omega="+str(omega)+",gamma="+str(gamma)+",r=3.npy")[-1])

R_SIR_th =np.load("evolution_R_SIR_th.csv") # evolution classic SIR model
R_SIR_sim =np.load("evolution_R_SIR_sim.csv") # evolution classic SIR model

#inferno colormap
cmap = plt.get_cmap('viridis')
#same but colored according to the colormap
colors = [cmap(0.7),cmap(0.5),cmap(0.2)]
#with smooth
plt.plot(Beta,evol_r1_th,label = r"$r=1$",color = colors[0],alpha = 1)
plt.plot(Beta,evol_r2_th,label = r"$r=2$",color = colors[1],alpha = 1)
plt.plot(Beta,evol_r3_th,label = r"$r=3$",color = colors[2],alpha = 1)
plt.scatter(Beta,evol_r1_sim,edgecolor = colors[0],color = "w")
plt.scatter(Beta,evol_r2_sim,edgecolor = colors[1],color = "w")
plt.scatter(Beta,evol_r3_sim,,edgecolor = colors[2],color = "w")
plt.scatter(np.linspace(0,1,len(R_SIR_sim[::5])),R_SIR_sim[::5],label = "SIR",color = "w",marker = "o",edgecolors="black")
plt.plot(np.linspace(0,1,len(R_SIR_th)),R_SIR_th,color = "black")

plt.legend(fontsize = 20)

plt.xlabel(r"$\beta$",fontsize = 25)
plt.ylabel(r"$R$",fontsize = 25)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
#grid
plt.grid(linewidth=0.2)

# %%
#Plot of the Heatmap

Block = False
Omega = [0.05,0.15,0.2,0.25,0.3,0.35]
Beta = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0]
evol_r = np.zeros((len(Omega),len(Beta)))
radius = 2
for i in range(len(Omega)):
    for j in range(len(Beta)):
        evol_r[i,j] = np.load(path_sim+"/block_vax/SF_2.5/EVOL_SIR-V/V_beta="+str(Beta[j])+",omega="+str(Omega[i])+",gamma="+str(gamma)+",r="+str(radius)+".npy")[-1] if Block else np.load(path_sim+"/ring_vax/SF_2.5/EVOL_SIR-V/R_beta="+str(Beta[j])+",omega="+str(Omega[i])+",gamma="+str(gamma)+",r="+str(radius)+".npy")[-1]

cmap = plt.get_cmap('viridis')
plt.figure(figsize=(20,20))
plt.imshow(evol_r, cmap=cmap, interpolation='nearest',vmin=0,vmax=1)
plt.yticks(np.arange(len(Omega)),Omega,fontsize = 20)
plt.xticks(np.arange(len(Beta)),Beta,fontsize = 20,rotation = 45)
plt.xlabel(r"$\beta$",fontsize = 25)
# %%
