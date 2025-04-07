#%%
import json
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
from geovoronoi import voronoi_regions_from_coords, points_to_coords
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from shapely.geometry import box
import contextily as cx
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
# %%
path = ""
puglia_prov = gpd.read_file(path+"/puglia_prov.geojson")
puglia_mun = gpd.read_file(path+"/puglia_mun.geojson")
uliveti = gpd.read_file(path+"/uliveti.geojson")
# %%
gallipoli = puglia_mun[puglia_mun['COMUNE'] == "Lecce"]
gallipoli = gallipoli.to_crs(epsg=3395)    # convert to World Mercator CRS
#%%
#substitute multipolygon with centroid
uliveti1 = uliveti.copy()
uliveti1 = uliveti1.to_crs(epsg=3395)
uliveti1['geometry'] = uliveti1.centroid
#find uliveto in gallipoli
#%%
uliveti_gallipoli = uliveti1[uliveti1.within(gallipoli.iloc[0].geometry)]
#%%
area = puglia_prov[puglia_prov['DEN_PROV'] == "Lecce"]
#area big if in Lecce, Brindisi or Taranto
area_big = puglia_prov[puglia_prov['DEN_PROV'].isin(["Lecce","Brindisi","Taranto"])]
#merge
area = area.to_crs(epsg=3395)    # convert to World Mercator CRS
area_big = area_big.to_crs(epsg=3395)
area_shape = area.iloc[0].geometry
#merge together the polygon of the three provinces in area big
area_shape_big = area_big.unary_union 
uliveti1 = uliveti1[uliveti1.within(area_shape)]

#%%
#compute voronoi regions
coords = points_to_coords(uliveti1.geometry)
region_polys, region_pts = voronoi_regions_from_coords(coords, area_shape.buffer(0.1))

#%%
#create dict with degree of each region
dict_degree = {}
j=0
for key,value in region_polys.items():
    lista = []
    for key1,value1 in region_polys.items():
        if value1.touches(value):
            lista.append(key1)
    dict_degree[key] = lista
    j+=1
    print(j)

#%%
#plot degree distribution
import numpy as np
import matplotlib.pyplot as plt
dict_degree = {int(k): v for k, v in dict_degree.items()}

prob_vec = [0]*(max([len(v) for k,v in dict_degree.items()])+1)
for k,v in dict_degree.items():
    prob_vec[len(v)] += 1
prob_vec = [i/len(dict_degree) for i in prob_vec]

plt.plot(prob_vec)
plt.show()

# %%
#PLOT GEO MAP WITH VORONOI REGIONS COLORED ACCORDING TO THE DEGREE


cmap = plt.cm.viridis
cmap = cmap(np.linspace(0, 1, 14))
import matplotlib.pyplot as plt
#plot area and voronoi, with centroids
fig, ax = plt.subplots(figsize=(10,10))   
area_big.plot(ax=ax, facecolor="none", edgecolor='grey')
n_xtics = len(ax.get_xticks())
n_ytics = len(ax.get_yticks())
min_lat = min(area_big.to_crs(epsg=4326).bounds["miny"])
max_lat = max(area_big.to_crs(epsg=4326).bounds["maxy"])
min_lon = min(area_big.to_crs(epsg=4326).bounds["minx"])
max_lon = max(area_big.to_crs(epsg=4326).bounds["maxx"])
lat = np.linspace(min_lat,max_lat,n_xtics)
lon = np.linspace(min_lon,max_lon,n_ytics)
ax.set_xticklabels([str(round(i,1)) for i in lon],fontsize = 13)
ax.set_yticklabels([str(round(i,1)) for i in lat],fontsize = 13)
#grid
ax.grid(linestyle = "--",alpha = 0.5)
#display ticks and labels
plt.xticks(fontsize = 13)
plt.yticks(fontsize = 13)
sea = gpd.GeoSeries(
    [
        box(*box(*area_shape_big.bounds).buffer(0.1).bounds).difference(
            area_shape_big
        )
    ]
)


xlim = ax.get_xlim()
ylim = ax.get_ylim()

plt.xlim((xlim[0]+100000,xlim[1]-10000))
plt.ylim((ylim[0]+7300,ylim[1]-50000))

sea.plot(ax=ax, facecolor="lightblue", edgecolor='lightblue')








for key,value in region_polys.items():
    if value.geom_type != 'Polygon':
        pass
    else:
        c = cmap[len(dict_degree[str(key)])-1]
        
        ax.plot(*value.exterior.xy, color=c, linewidth=0.5)
        ax.fill(*value.exterior.xy, color=c, alpha=0.5) 

        
    
        ax.plot(*coords[key], 'o', color= c, markersize=1)
cmap = ListedColormap(cmap)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=14))
sm._A = []


plt.show()

# %%
dict_degree = {int(k): v for k, v in dict_degree.items()}
prob_vec = [0]*(max([len(v) for k,v in dict_degree.items()])+1)
for k,v in dict_degree.items():
    prob_vec[len(v)] += 1
prob_vec = [i/len(dict_degree) for i in prob_vec]
G = nx.Graph(dict_degree)
#%%
#PLOT THE MAP COLORED ACCORDING TO THE STATE OF NODES

state = pd.read_csv(path+"SIRV-spatial.csv")
omega = 0.2
beta = 1
I = []
for i in I:
    stato = state[i]
    fig, ax = plt.subplots(figsize=(10,10))   
    #ax.axis("off")
    #change epsg  
    #plot only the area of Lecce and a bit of Brindisi and Taranto (neighbouring provinces)
    area_big.plot(ax=ax, facecolor="none", edgecolor='grey')
    n_xtics = len(ax.get_xticks())
    n_ytics = len(ax.get_yticks())
    min_lat = min(area_big.to_crs(epsg=4326).bounds["miny"])
    max_lat = max(area_big.to_crs(epsg=4326).bounds["maxy"])
    min_lon = min(area_big.to_crs(epsg=4326).bounds["minx"])
    max_lon = max(area_big.to_crs(epsg=4326).bounds["maxx"])
    lat = np.linspace(min_lat,max_lat,n_xtics)
    lon = np.linspace(min_lon,max_lon,n_ytics)
    ax.set_xticklabels([str(round(i,1)) for i in lon],fontsize = 13)
    ax.set_yticklabels([str(round(i,1)) for i in lat],fontsize = 13)
    #grid
    ax.grid(linestyle = "--",alpha = 0.5)
    #display ticks and labels
    plt.xticks(fontsize = 25)
    plt.yticks(fontsize = 25)
    #add labels
    plt.xlabel("Longitude (°)",fontsize = 25) 
    plt.ylabel("Latitude (°)",fontsize = 25) if i == 45 else plt.ylabel("")
    sea = gpd.GeoSeries(
        [
            box(*box(*area_shape_big.bounds).buffer(0.1).bounds).difference(
                area_shape_big
            )
        ]
    )


    #get currend x and y limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    plt.xlim((xlim[0]+100000,xlim[1]-10000))
    plt.ylim((ylim[0]+7300,ylim[1]-50000))

    sea.plot(ax=ax, facecolor="lightblue", edgecolor='lightblue')


    for key,value in region_polys.items(): #voronoi regions dictionary
        if value.geom_type != 'Polygon':
            pass
        else:
            if stato[key] == 1:
                c = 'blue'
            elif stato[key] == 2:
                c = 'red'
            elif stato[key] == 3:
                c = "gray"
            else:    
                c = 'green'
            ax.plot(*value.exterior.xy, color=c, linewidth=0.5)
            ax.fill(*value.exterior.xy, color=c, alpha=0.5) #plot external part of polygon

            
        
            ax.plot(*coords[region_pts[key][0]], 'o', color= c, markersize=1)
            #
    #legend
    custom_lines = [Line2D([0], [0], color='green', lw=4),
                    Line2D([0], [0], color='blue', lw=4),
                    Line2D([0], [0], color='red', lw=4),
                    Line2D([0], [0], color='gray', lw=4)]



    





    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    


    sea.plot(ax=ax, facecolor="lightblue", edgecolor='lightblue')

#%%

#load prob_state
Ring = False
prob_state_block = np.load(path+"prob_state_dynblock_omega0.25beta0.7r2_.npy")
prob_state_ring = np.load(path+"prob_state_dynring_omega0.25beta0.7r2_.npy")
niter = 300
#now I want to link each node with the probability of being recovered at the end of the simulation.
prob_state_block = [i/niter for i in prob_state_block]
prob_state_ring = [i/niter for i in prob_state_ring]

#plot the probability of being in state 3 for each node

ax1 = plt.subplot(121)
#change size of ax1
ax2 = plt.subplot(122)
#change total figure size
fig = plt.gcf()

#make the two plots bigger
fig.set_size_inches(20,10)

#reduce space between plots
plt.subplots_adjust(wspace = 0,hspace = 0)




#cmap infero_r from 0 to 1
cmap = plt.cm.inferno_r


 
#plot only the area of Lecce and a bit of Brindisi and Taranto (neighbouring provinces)
area_big.plot(ax=ax2, facecolor="none", edgecolor='grey')
n_xtics = len(ax2.get_xticks())
n_ytics = len(ax2.get_yticks())
min_lat = min(area_big.to_crs(epsg=4326).bounds["miny"])
max_lat = max(area_big.to_crs(epsg=4326).bounds["maxy"])
min_lon = min(area_big.to_crs(epsg=4326).bounds["minx"])
max_lon = max(area_big.to_crs(epsg=4326).bounds["maxx"])
lat = np.linspace(min_lat,max_lat,n_xtics)
lon = np.linspace(min_lon,max_lon,n_ytics)
ax2.set_xticklabels([str(round(i,1)) for i in lon],fontsize = 13)
ax2.set_yticklabels([str(round(i,1)) for i in lat],fontsize = 13)
#grid
ax2.grid(linestyle = "--",alpha = 0.5)
#display ticks and labels
plt.xticks(fontsize = 13)
plt.yticks(fontsize = 13)
sea = gpd.GeoSeries(
    [
        box(*box(*area_shape_big.bounds).buffer(0.1).bounds).difference(
            area_shape_big
        )
    ]
)


#get currend x and y limits
xlim = ax2.get_xlim()
ylim = ax2.get_ylim()

ax2.set_xlim((xlim[0]+100000,xlim[1]-10000))
ax2.set_ylim((ylim[0]+7300,ylim[1]-50000))


sea.plot(ax=ax2, facecolor="lightblue", edgecolor='lightblue')


for key,value in region_polys.items():
    #if value not multipolygon
    if value.geom_type != 'Polygon':
        pass
    else:
        c = cmap(prob_state_block[key])
        #plot borders of polygon and fill with color
        ax2.plot(*value.exterior.xy, color=c, linewidth=0.5)
        ax2.fill(*value.exterior.xy, color=c, alpha=0.5) #plot external part of polygon

        #plot internal part of polygon
        ax2.plot(*coords[region_pts[key][0]], 'o', color= c, markersize=1)
        #







area_big.plot(ax=ax1, facecolor="none", edgecolor='grey')
n_xtics = len(ax1.get_xticks())
n_ytics = len(ax1.get_yticks())
min_lat = min(area_big.to_crs(epsg=4326).bounds["miny"])
max_lat = max(area_big.to_crs(epsg=4326).bounds["maxy"])
min_lon = min(area_big.to_crs(epsg=4326).bounds["minx"])
max_lon = max(area_big.to_crs(epsg=4326).bounds["maxx"])
lat = np.linspace(min_lat,max_lat,n_xtics)
lon = np.linspace(min_lon,max_lon,n_ytics)
ax1.set_xticklabels([str(round(i,1)) for i in lon],fontsize =20)
ax1.set_yticklabels([str(round(i,1)) for i in lat],fontsize = 20)
#labels for x and y
#y label only for first plot
ax1.set_xlabel("Longitude (°)",fontsize = 25)
ax1.set_ylabel("Latitude (°)",fontsize = 25)

ax2.set_xlabel("Longitude (°)",fontsize = 25)
#grid
ax1.grid(linestyle = "--",alpha = 0.5)
#display ticks and labels
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
sea = gpd.GeoSeries(
    [
        box(*box(*area_shape_big.bounds).buffer(0.1).bounds).difference(
            area_shape_big
        )
    ]
)


#get currend x and y limits
xlim = ax1.get_xlim()
ylim = ax1.get_ylim()

ax1.set_xlim((xlim[0]+100000,xlim[1]-10000))
ax1.set_ylim((ylim[0]+7300,ylim[1]-50000))

sea.plot(ax=ax1, facecolor="lightblue", edgecolor='lightblue')


for key,value in region_polys.items():
    #if value not multipolygon
    if value.geom_type != 'Polygon':
        pass
    else:
        c = cmap(prob_state_ring[key])
        #plot borders of polygon and fill with color
        ax1.plot(*value.exterior.xy, color=c, linewidth=0.5)
        ax1.fill(*value.exterior.xy, color=c, alpha=0.5) #plot external part of polygon

        #plot internal part of polygon
        ax1.plot(*coords[region_pts[key][0]], 'o', color= c, markersize=1)
        #

#add colorbar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1))
sm._A = []

#cbar common to both plots
cbar = plt.colorbar(sm, ax=[ax1,ax2],orientation = 'horizontal')

#reduce length of colorbar



cbar.set_label('Probability of being Recovered',fontsize = 25)
#ticks
cbar.ax.tick_params(labelsize=20)

# %%
####### VARIATION OF THE SPATIAL NETWORK: EDGES DEPENDING ON DISTANCE BETWEEN TREES
from scipy.spatial.distance import squareform, pdist
from math import sin, cos, sqrt, atan2, radians
import pandas as pd

def ProjCoordsTangentSpace(lat,lon,lat0,lon0):
    '''
    Description:
    Projects in the tangent space of the earth in (lat0,lon0) 
    Return: 
    The projected coordinates of the lat,lon  '''
    PI = np.pi
    c_lat= 0.6*100000*(1.85533-0.006222*np.sin(lat0*PI/180))
    c_lon= c_lat*np.cos(lat0*PI/180)
    
    x = c_lon*(lon-lon0)
    y = c_lat*(lat-lat0)
    if isinstance(x,np.ndarray) or isinstance(x,np.float64):
        pass
    else:        
        x = x.to_numpy()
    if isinstance(y,np.ndarray) or isinstance(y,np.float64):
        pass
    else:
        y = y.to_numpy()
    return x,y
# %%
#apply ProjCoordsTangentSpace to uliveti1
uliveti1['x'], uliveti1['y'] = ProjCoordsTangentSpace(uliveti1.geometry.y, uliveti1.geometry.x, uliveti1.geometry.y.mean(), uliveti1.geometry.x.mean())
# %%
#calculate distance matrix between the new coordinates
dist_matrix = distance_matrix([[i["x"],i["y"]] for i in uliveti1.to_dict("records")], [[i["x"],i["y"]] for i in uliveti1.to_dict("records")])
# %%
deg_dist = []
distances = [200,300,400,500,600,700,800,900,1000]
list_giant_comp = []
dict_dict_degrees = {}
#for each distance, build a network with the points that are at distance less than the distance
#calculate degree distribution and add it to deg_dist
for d in distances:
    print(d)
    dist_matrix_1 = dist_matrix < d
    #set True to 1 
    dist_matrix_1 = dist_matrix_1.astype(int)
    G = nx.from_pandas_adjacency(pd.DataFrame(dist_matrix_1))
    #get size of giant component
    giant = max(nx.connected_components(G), key=len)
    list_giant_comp.append(len(giant))
    print("edges added")
    dict_degree = dict(G.degree())
    dict_degree = {int(k): v for k, v in dict_degree.items()}
    #get a prob_vec from dict_degree
    #count the number of nodes with degree k
    prob_vec = [0 for i in range(max(dict_degree.values())+1)]
    for k,v in dict_degree.items():
        prob_vec[v] += 1
    #normalize
    prob_vec = [i/len(dict_degree) for i in prob_vec]
    #deg_dist.append(prob_vec)
    dict_dict_degrees[str(d)] = dict_degree

#%%
plt.plot(distances, [i/len(uliveti1) for i in list_giant_comp])
plt.xlabel("Threshold (meters)", fontsize = 20)
plt.ylabel("Giant component size",fontsize = 20)
#white grid
plt.grid(color = 'black', linestyle = '--', linewidth = 0.3,alpha = 0.4)