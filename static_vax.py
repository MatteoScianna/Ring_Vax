#%%
from functions import *
from pynverse import inversefunc
from sympy import *
from sympy import symbols
from scipy.stats import *
from scipy.optimize import fsolve
import numpy as np

#%%
def static_vax(n_nodes,perc_vax,prob_vec,radius,Ring):


    remaining_nodes = n_nodes
    #before: n = n_nodes fixed. Hence, where there is n should be n_nodes
    alpha = 1
    t = 0
    p = np.zeros((n_nodes+1, len(prob_vec)))
    d = np.zeros((n_nodes+1, len(prob_vec)))
    q = np.zeros((n_nodes+1, len(prob_vec)))
    k_avg = np.zeros(n_nodes+1)
    r_avg = np.zeros(n_nodes+1)
    S = []
    # Set initial values
    p[0] = prob_vec
    q[-1] = p[0]
    q[0] = p[0]
    d[1] = p[0]
    k_avg[0] = np.sum(np.arange(len(prob_vec))*p[0])
    r_avg[1] = np.sum(np.arange(len(prob_vec))*p[0])
    removed_nodes = 0
    size_gc = 1
    while remaining_nodes > n_nodes*(1-perc_vax) and size_gc > 0.005:
        
        t+=1

    #################### teleportation #####################


        for k in range(len(prob_vec)): #for each degree

            p[t][k] = (1 / (n_nodes-t)) * (p[t-1, k] * (n_nodes-t+1) - d[t, k] +  #update p
                                            r_avg[t] * (q[t-1, k] - (q[t-1, k-1] if k > 0 else 0)))
            
        k_avg[t] = np.sum(np.arange(len(prob_vec)) * p[t]) #consequently, update k_avg

        for k in range(len(prob_vec)): #and update q, d and r_avg
            q[t][k] = ((k+1)*p[t,k+1]) / k_avg[t] if k < len(prob_vec)-1 else 0
            d[t+1, k] = p[t, k]
        r_avg[t+1] = np.sum(np.arange(len(prob_vec)) * d[t+1])

        remaining_nodes = n_nodes-t

        #solve the equation to find the size of the giant component after the teleportation 

        x_intersect = fsolve(lambda x: g_1(x,p[t],k_avg[t])-x,0.05)
        size_gc = (1-g_0(x_intersect,p[t]))[0]
        S.append((remaining_nodes/n_nodes,(remaining_nodes/n_nodes)*size_gc))
        print("time: ",t,"remaining nodes: ",remaining_nodes,"size GC: ",size_gc)

        copy_prob_dist = p[t].copy() #useful because t will be updated but we need the current p[t] for the next step

        


    ############## vaccination ####################

        #update metrics
        mean_deg = sum([i*p[t][i] for i in range(len(p[t]))])
        sq_mean_deg = sum([i**2*p[t][i] for i in range(len(p[t]))])
        removed_nodes = sum([((((sq_mean_deg-mean_deg)/mean_deg)**(i-1))*mean_deg) for i in range(1,radius+1)])
        reintroduced_nodes = sum([((((sq_mean_deg-mean_deg)/mean_deg)**(i-1))*mean_deg) for i in range(radius)])
        remaining_nodes = (n_nodes-t)-removed_nodes
        perc_remaining_nodes = remaining_nodes/(n_nodes-t)


        #solve equation for the size of the giant component after vaccination (beforre eventual removal)
        g_01 = (lambda x: sum([p[t][i]*(x**(i)) for i in range(len(p[t]))]))
        f = inversefunc(g_01,y_values = perc_remaining_nodes)
        #print("new f",f)
        x = symbols('x')
        g_12 = (sum([p[t][i]*(x**(i)) for i in range(len(p[t]))]))
        der = diff(g_12,x)
        gen_fun =(1/g_0(f,p[t]))*(g_0(f+(der.subs(x,f)/der.subs(x,1))*(x-1),p[t]))
        degree_distribution_itermediate = []
        for j in range(0,(len(p[t]))):
            df = diff(gen_fun, x, j)
            deg1 = (1/factorial(j))*df.subs(x,0).evalf()
            degree_distribution_itermediate.append(deg1)
        t+=removed_nodes
        t = int(round(t))


        ## Intermediate step only for ring vaccination: reintroduction of nodes ##
        if Ring:
            prob_intermediate = [float(k) for k in degree_distribution_itermediate]
            list_degrees = [(n_nodes-t)*i for i in prob_intermediate]

            #add number of 0 degree nodes equal to the number of removed nodes belonging to internal shells
            list_degrees[0]+= sum([(((sq_mean_deg-mean_deg)/mean_deg)**(i-1))*mean_deg for i in range(radius)])
            list_degrees = [list_degrees[i]/sum(list_degrees) for i in range(len(list_degrees))]
            
            #back in time! we need to reintroduce the nodes that were removed
            t -= int(round(reintroduced_nodes))
            p[t] = list_degrees



        else: # if not ring vaccination we can directly update the degree distribution
            p[t] = [float(k) for k in degree_distribution_itermediate]

        k_avg[t] = np.sum(np.arange(len(prob_vec)) * p[t])


        #update q, d and r_avg
        for k in range(len(prob_vec)):
            q[t][k] = ((k+1)*p[t,k+1]) / k_avg[t] if k < len(prob_vec)-1 else 0
            d[t+1, k] = d[t, 0] * 1*p[t, k] + (1 - d[t, 0]) * (alpha *1*p[t, k] + (1-alpha) * q[t-1, k])
        r_avg[t+1] = np.sum(np.arange(len(prob_vec)) * d[t+1])



    ### Now we have the final degree distribution after vaccination, we can calculate the final size of the giant component ###



    #find last row of p with non zero elements (before we defined p with zeros!!! hence we need to find the last row with non zero elements)
    for i in range(len(p)-1,0,-1):
        if sum(p[i])>0:
            last_row = i
            break

    #this is the probability degree distribution of the network after vaccination 
    prob_vec_final = p[last_row]


    #solve the equation to find the size of the giant component after the whole process

    mean= sum([prob_vec_final[i]*i for i in range(len(prob_vec_final))])
    size_gc_final = []
    step = np.linspace(0,1,100)
    for i in range(len(step)):
        x_intersect = fsolve(lambda x: 1-step[i]+step[i]*g_1(x,prob_vec_final,mean) - x,0.05)
        size_gc_final.append(float(list(abs((1-g_0(x_intersect,prob_vec_final)))*(1-perc_vax))[0]))
    return size_gc_final, prob_vec_final


#%%

##### INITIALIZATION THEORY #######


n_nodes = 10000 #nodes we want to have in the graph 
radius = 2 #radius of the ring
net = "SF2.5"
Ring = True
#prob_vec = prob_list_poisson(6,20)
#prob_vec = [0]*20
#rob_vec[6] = 1
prob_vec = prob_list_powerlaw(2.5,3,400)
#degree_vec_in = discrete_samples(prob_vec,n_nodes)
Perc_vax = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
P = np.zeros((len(Perc_vax)+1,len(prob_vec)))
P[0] = prob_vec
prob_vec_final = prob_vec
path_th = ""
#%%
#Now iterate for all the percentages of vaccinated nodes considered
for perc_vax in Perc_vax:
    GC,prob_vec_final = static_vax(n_nodes,perc_vax,prob_vec_final,radius,Ring)

# %%

####### INITIALIZATION SIMULATION #######

path_sim = ""
Perc_vax = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
niter = 1000
iteraz = 0
Ring = True
alpha = 1
radius = 2
net = "SF2.5"
prob_vec = prob_list_powerlaw(2.5,3,400)
p = np.zeros((len(Perc_vax)+1,len(prob_vec)))
s = np.zeros((len(Perc_vax)+1,100))
nnodes = 10000
degree_vec_in = discrete_samples(prob_vec,nnodes)
graph = config_model(degree_vec_in)
p[0] = prob_vec

#%%

####### STATIC VAX - SIMULATION ########

while iteraz < niter:
    print(iteraz)
    j = 0
    graph_copy = graph.copy()
    for perc_vax in Perc_vax:
        j+=1
        print(j)
        
        while len(graph_copy.nodes()) > nnodes*(1-perc_vax) and len(max(nx.connected_components(graph_copy), key=len))/nnodes>0.005:
            node = rand.choice(list(graph_copy.nodes()))
            tot_neigh = nx.single_source_shortest_path_length(graph_copy,node,radius)
            neigh = set(key for key, value in tot_neigh.items() if value == radius) if Ring else set(key for key, value in tot_neigh.items() if value <= radius)
            for nei in neigh:
                if rand.random() < alpha:
                    graph_copy.remove_node(nei)
            if Ring:
                neigh_copia = tot_neigh.copy()
                for node in neigh_copia:
                    lista_edges_to_remove = []
                    if node not in neigh:
                        #remove all edges of node
                        lista_edges_to_remove.append([(node,neigh) for neigh in graph_copy.neighbors(node)])
                    for i in lista_edges_to_remove:
                        graph_copy.remove_edges_from(i)




        deg_dist = [graph_copy.degree(i) for i in graph_copy.nodes()]
        prob_vec_final = prob_list(deg_dist,len(prob_vec))
        p[j] += prob_vec_final

        mean= sum([prob_vec_final[i]*i for i in range(len(prob_vec_final))])
        S2 = []
        step = np.linspace(0,1,100)
        for i in range(len(step)):
            x_intersect = fsolve(lambda x: 1-step[i]+step[i]*g_1(x,prob_vec_final,mean) - x,0.05)
            S2.append(float(list(abs((1-g_0(x_intersect,prob_vec_final)))*(1-perc_vax))[0]))
        s[j] += S2
    iteraz += 1

p = [i/iteraz for i in p]
p[0] = prob_vec
s = [i/iteraz for i in s]
    

# %%
#LOADING DATA TO PLOT: DEGREE DISTRIBUTION EVOLUTION
net = "ER6"
Perc_vax = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
radius =3
Ring = True
p_th = np.load(path_th+'/ring_vax/'+net+'/evol_prob_vec/prob_vec_r'+str(radius)+'.npy',allow_pickle=True) if Ring else np.load(path_th+'/block_vax/'+net+'/evol_prob_vec/prob_vec_r'+str(radius)+'.npy',allow_pickle=True)
p_sim = np.load(path_sim+'/ring_vax/'+net+'/evol_prob_vec/prob_vec_r'+str(radius)+'.npy',allow_pickle=True) if Ring else np.load(path_sim+'/block_vax/'+net+'/evol_prob_vec/prob_vec_r'+str(radius)+'.npy',allow_pickle=True)
#%%
#PLOT EVOLUTION OF DEGREE DISTRIBUTION - THEORY VS SIMULATION

import matplotlib.pyplot as plt
import matplotlib.cm as cm
colors = cm.cividis(np.linspace(0, 1, len(Perc_vax)))
for i in range(0,len(p_th)-1):
    plt.plot(p_th[i],color = colors[i],alpha = 0.5,label = r"$\phi$ = "+str(Perc_vax[i]))
    plt.plot(p_sim[i],linestyle = 'dashed',color = colors[i],alpha = 0.5)
plt.plot([],label = "Theory",color = 'black')
plt.plot([],label = "Simulation",linestyle = 'dashed',color = 'black')
plt.xlabel(r'k',fontsize = 20)
plt.ylim(0,0.6)
plt.xlim(-0.5,15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 0)

plt.grid(True,which = 'both',linestyle = '--',linewidth = 0.2)
plt.tight_layout()

# %%
#LOAD DATA: EVOLUTION OF GIANT COMPONENT
s_sim = np.load(path_sim+'/ring_vax/'+net+'/evol_impact/p_S_r'+str(radius)+'.npy',allow_pickle=True) if Ring else np.load(path_sim+'/block_vax/'+net+'/evol_impact/p_S_r'+str(radius)+'.npy',allow_pickle=True)
S = np.zeros((len(Perc_vax),100))
for perc_vax in Perc_vax:
    S1 = np.load(path_th+'/ring_vax/'+net+'/evol_impact/p_S_perc_vax'+str(perc_vax)+'r'+str(radius)+'.npy',allow_pickle=True) if Ring else np.load(path_th+'/block_vax/'+net+'/evol_impact/p_S_perc_vax'+str(perc_vax)+'r'+str(radius)+'.npy',allow_pickle=True)
    S[Perc_vax.index(perc_vax)] = S1[1]
#%%
#PLOT EVOLUTION OF IMPACT - THEORY VS SIMULATION

S_SF_block = np.load_csv(path_th+"/block/S_SF.csv")
s_sim_SF_block = np.load_csv(path_th+"S_SF_sim.csv")
S_SF_ring = np.load_csv(path_th+"/ring/S_SF.csv")
s_sim_SF_ring = np.load_csv(path_th+"/ring/S_SF_sim.csv")


#%%
import matplotlib.pyplot as plt
import seaborn as sns  


fig, ax = plt.subplots(1,2, figsize=(45,20), dpi = 200)
sns.set_style("whitegrid")
plt.rcParams['axes.grid'] = True

perc_phi = [0.1,0.3,0.5]
colors = cm.Blues([0.3,0.7,2])
#####################
####### COLUMN 1
####################

for i in range(len(S_SF_block)):
    ax[0].plot(np.linspace(0,1,len(S_SF_block[i])),S_SF_block[i],color = colors[i],alpha = 1,label = r"$\phi$ = "+str(perc_phi[i]),linewidth = 3)
    ax[0].scatter(np.linspace(0,1,len(s_sim_SF_block[i]))[::3],s_sim_SF_block[i][::3],edgecolor = colors[i],alpha = 1,facecolor = 'none',s = 200)



#ax[0,0].set_title('a ( 'r'$sec^{1/2}$ ) : Threshold',fontsize = 30, pad = 15)
ax[0].set_xlabel(r'$\beta$',fontsize = 40)
ax[0].set_ylabel('R', labelpad = 15, fontsize = 40)
ax[0].tick_params(axis ='x', labelsize = 28)
ax[0].tick_params(axis ='y', labelsize = 28)
#ax[0].grid(axis = 'x')
#legend
ax[0].legend(fontsize = 40)


#reduce space between subplots
plt.subplots_adjust(wspace = 0.1)

#####################

###### COLUMN 2

for i in range(len(S_SF_ring)):
    ax[1].plot(np.linspace(0,1,len(S_SF_ring[i])),S_SF_ring[i],color = colors[i],alpha = 1,linewidth = 3)
    ax[1].scatter(np.linspace(0,1,len(s_sim_SF_ring[i]))[::3],s_sim_SF_ring[i][::3],label = r"$\phi$ = "+str(perc_phi[i]),edgecolor = colors[i],alpha = 1,facecolor = 'none',marker = 'o',s = 200)

#ax[0,1].set_title('z : Bias', fontsize = 30, pad = 15)
ax[1].set_xlabel(r'$\beta$',fontsize = 40)
ax[1].set_ylabel(' ')
#ax[1].set_xticks([])
#ax[1].set_ylim(0,1)
#ax[1].s(labelsize = 0)

ax[1].tick_params(axis = 'y', labelsize = 0)
ax[1].tick_params(axis = 'x', labelsize = 28)
#ax[1].get_legend().remove()

#####################

#####################
