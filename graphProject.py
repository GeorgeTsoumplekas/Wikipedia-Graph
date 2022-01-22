import networkx as nx
import wikipedia as wp
from timeit import default_timer as timer
import matplotlib.pyplot as plt
import math

"""
Function that simply returns the input the user gave for the word he wants to look for
"""
def get_query_word():
    word = input("Enter the subject you want to look for: ")
    
    return word  


"""
Function creating the directed graph.
The graph is constructed in such a way that links with similar names that redirect
the user to the same article appear as one node with the same name. (done with wp.search)
The depth level of the graph is not predefined but when a number of nodes has already
been added to the graph we stop continuing searching in the next depth level.
"""
def create_graph():
    G = nx.DiGraph()
    total = 0
    
    #list of pages in each depth level
    page_list = list()
    
    search_word = list()
    
    word = get_query_word()
    search_word.append((wp.search(word,1))[0])
    
    while len(search_word[0]) == 0:
        print("Can't find a wikipedia page for this query. Please try again!")
        print("Note: check for typos and try using lower case only letters")
        print("(Except for names and the first letter of the first word)")
        search_word.clear()
        word = get_query_word()   
        search_word.append((wp.search(word,1))[0])
        
    start = timer()
    
    print("Title is: " + search_word[0])
    
    G.add_node(search_word[0])
    total += 1
    
    #list keeping the wikipedia pages we have in each depth level
    page_list.append(wp.page(search_word,auto_suggest=False))
    
    print("Number of links in initial search page is: {}".format(len(page_list[0].links)))
    print("Please wait, this might take a while..")
    
    #add additional depth levels to the graph setting a limit to stop adding more
    #nodes after we have exceeded 1000    
    while total < 1000:
        if total != 1:
            page_list.clear()
            for s_word in search_word:
                try:
                    pg = wp.page(s_word,auto_suggest=False)
                    #remove 'International Standard Book Number' from the graph on purpose
                    #since it exists in almost every page and it is not related to any subject
                    if pg.title == 'International Standard Book Number':
                        continue
                    page_list.append(pg)
                except wp.DisambiguationError as e:
                    pg = wp.page(e.options[0],auto_suggest=False)
                    if pg.title == 'International Standard Book Number':
                        continue
                    page_list.append(pg)
                except wp.PageError:
                    continue
        search_word.clear()
        for page in page_list:
            for lnk in page.links:
                poss_word = wp.search(lnk,1)
                if poss_word == []:
                    continue
                else:
                    search_word.append(poss_word[0])
                
                #if no article can be found for this link just skip it
                if len(search_word[-1]) == 0:
                    search_word.pop()
                    continue
                else:
                    G.add_edge(page.title,search_word[-1])
                    total += 1
    
    #remove it again, just to be sure
    blocked_node = "International Standard Book Number"
    if blocked_node in G:
        G.remove_node(blocked_node)
    
    end = timer()
    print("The graph has been created")
    print("Elapsed time: {} s".format(end-start))
    
    return G


"""
Function that finds the number of nodes and edges in the graph
"""
def get_nodes_edges(G):
    print("\nTotal number of nodes is: {}".format(nx.number_of_nodes(G)))
    print("Total number of edges: {}".format(nx.number_of_edges(G)))


"""
Function that finds the nodes adjacent to a specific node of the graph
""" 
def get_node_neighbor(G):
    node_name = input("Name of the node you want to find its neighbors: ")
    while node_name not in G:
        print("Node not found. Try again.")
        node_name = input("Name of the node you want to find its neighbors: ")
        
    neighbors = nx.neighbors(G,node_name)
    print("For node {} its neighbors are:".format(node_name))
    for n in neighbors:
        print(n)


"""
Function that finds the in and out degree of a specific node
"""
def get_node_degree(G):
    node_name = input("Name of the node you want to find its degree: ")
    while node_name not in G:
        print("Node not found. Try again.")
        node_name = input("Name of the node you want to find its degree: ")
        
    in_degree = G.in_degree(node_name)
    print("For node {} in-degree is: {}".format(node_name,in_degree))
    
    out_degree = G.out_degree(node_name)
    print("For node {} out-degree is: {}".format(node_name,out_degree))
    
 
"""
Function that finds the node with the biggest degree in the graph
If the given graph is directed it looks for both the node with the biggest in-degree
and the node with the biggest out-degree
In case of many nodes with the same degree it returns the first of them found
"""
def get_biggest_degree(G):
    if G.is_directed():
        big_node_in = None
        big_node_out = None
        big_degree_in = -1
        big_degree_out = -1
        
        for node in G.nodes():
            if G.in_degree(node) > big_degree_in:
                big_degree_in = G.in_degree(node)
                big_node_in = node
            if G.out_degree(node) > big_degree_out:
                big_degree_out = G.out_degree(node)
                big_node_out = node
        
        print("\nBiggest in-degree: {} in node ".format(big_degree_in) + str(big_node_in))
        print("Biggest out-degree: {} in node ".format(big_degree_out) + str(big_node_out))
        
    else:
        big_node = None
        big_degree = -1
        
        for node in G.nodes():
            if G.degree(node) > big_degree:
                big_degree = G.degree(node)
                big_node = node
        
        print("Biggest degree: {} in node ".format(big_degree) + str(big_node))


"""
Function that finds the node with the smallest degree in the graph
If the given graph is directed it looks for both the node with the smallest in-degree
and the node with the smallest out-degree
In case of many nodes with the same degree it returns the first of them found
"""
def get_smallest_degree(G):
    if G.is_directed():
        small_node_in = None
        small_node_out = None
        small_degree_in = 1000000
        small_degree_out = 1000000
        
        for node in G.nodes():
            if G.in_degree(node) < small_degree_in:
                small_degree_in = G.in_degree(node)
                small_node_in = node
            if G.out_degree(node) < small_degree_out:
                small_degree_out = G.out_degree(node)
                small_node_out = node
                
        print("\nSmallest in-degree: {} in node ".format(small_degree_in) + str(small_node_in))
        print("Smallest out-degree: {} in node ".format(small_degree_out) + str(small_node_out))
        
    else:
        small_node = None
        small_degree = 1000000
        
        for node in G.nodes():
            if G.degree(node) < small_degree:
                small_degree = G.degree(node)
                small_node = node
        
        print("Smallest degree: {} in node ".format(small_degree) + str(small_node))    
 

"""
Function creating the degree distribution plots of the given graph
If the graph is directed it creates both the in-degree distribution plot
and the out-degree distribution plot
"""       
def degree_distribution(G,graph_name):
    if G.is_directed():
        in_degrees = [G.in_degree(n) for n in G.nodes()]
        out_degrees = [G.out_degree(n) for n in G.nodes()]
        
        plt.hist(in_degrees,label='in-degree')
        plt.title("Distribution of in degrees in " + graph_name)
        plt.xlabel("in-degree")
        plt.ylabel("# of nodes")
        plt.legend()
        plt.show()
        
        plt.hist(out_degrees, label='out-degree')
        plt.title("Distribution of out degrees in " + graph_name)
        plt.xlabel("out-degree")
        plt.ylabel("# of nodes")
        plt.legend()
        plt.show()
        
    else:
        degrees = [G.degree(n) for n in G.nodes()]
        plt.hist(degrees)
        plt.title("Distribution of degrees in " + graph_name)
        plt.xlabel("degree")
        plt.ylabel("# of nodes")
        plt.show()
 

"""  
Function that finds the average degree of the nodes of the graph
If the graph is directed it finds both average in-degree and average out-degree
"""  
def average_degree(G):
    nnodes = G.number_of_nodes()
    
    if nx.is_directed(G):        
        avg_degree_in = sum(d for n, d in G.in_degree()) / float(nnodes)
        avg_degree_out = sum(d for n, d in G.out_degree()) / float(nnodes)
        
        print("\nAverage in-degree is {}".format(avg_degree_in))
        print("Average out-degree is {}".format(avg_degree_out))
        
    else:
        avg_degree = sum(d for n, d in G.degree()) / float(nnodes)
        
        print("Average degree is {}".format(avg_degree))
        
    
"""
Function that calculates the density of a *directed* graph
"""
def get_density(G):
    nedges = float(G.number_of_edges())
    nnodes = float(G.number_of_nodes())
    density = nedges/(nnodes*(nnodes-1)) 
    print("\nDensity = {}".format(density))
    

"""
Function that finds the successors of a specific node in the given directed graph
"""
def get_successors(G):
    node_name = input("Name of the node you want to find its successors: ")
    while node_name not in G:
        print("Node not found. Try again.")
        node_name = input("Name of the node you want to find its successors: ")

    succs = G.successors(node_name)
    
    print("For node " + node_name + " its succesors are:")
    for s in succs:
        print(s)

"""
Function that finds the predecessors of a specific node in the given directed graph
"""
def get_predecessors(G):
    node_name = input("Name of the node you want to find its predecessors: ")
    while node_name not in G:
        print("Node not found. Try again.")
        node_name = input("Name of the node you want to find its predecessors: ")

    preds = G.predecessors(node_name)
    
    print("For node " + node_name + " its predecessors are:")
    for s in preds:
        print(s)       
        
"""
Function that checks if the given directed graph is strongly and/or weakly connected
If not it finds the number of its strongly and/or weakly connected components
"""        
def check_components(G):
    if nx.is_strongly_connected(G):
        print("\nThe graph is strongly connected")
    else:
        print("The graph is not strongly connected")
        scc = nx.number_strongly_connected_components(G)
        print("The graph has {} strongly connected components".format(scc))
    
    if nx.is_weakly_connected(G):
        print("The graph is weakly connected")
    else:
        print("The graph is not weakly connected")
        wcc = nx.number_weakly_connected_components(G)
        print("The number of weakly connected components is {}".format(wcc))
    

"""
Function that calculate different centrality metrics for a specific node in the given directed graph
More specifically it calculates the in-degree centrality, out-degree centrality 
and closeness centrality of the node
"""
def node_centrality(G):
    node_name = input("Name of the node you want to find its centrality: ")
    while node_name not in G:
        print("Node not found. Try again.")
        node_name = input("Name of the node you want to find its centrality: ")
        
    in_centr = nx.in_degree_centrality(G)
    out_centr = nx.out_degree_centrality(G)
    
    print("For node " + node_name + ":")
    print("The in-degree centrality is: {}".format(in_centr.get(node_name)))
    print("The out-degree centrality is: {}".format(out_centr.get(node_name)))
    print("The closeness centrality is: {}".format(nx.closeness_centrality(G,node_name)))


"""
Function that calculates the eccentricity of a specific node in the given directed graph
Special care taken if the graph is not strongly connected
"""
def get_eccentricity(G):
    node_name = input("Name of the node you want to find its eccentricity: ")
    while node_name not in G:
        print("Node not found. Try again.")
        node_name = input("Name of the node you want to find its eccentricity: ")
      
    if nx.is_strongly_connected(G):
        ecc = nx.eccentricity(G,node_name)
        print("The eccentricity of node " + node_name + " is {}".format(ecc))
    else:
        print("Cannot calculate eccentricity because the graph is not strongly connected.")
    
"""
Function that calculates various distance metrics of the graph such as radius,
diameter,center and average shortest path length.
If the given graph is directed and is not strongly connected we use its undirected
equivalent (for which we are sure that it is connected due to the way the directed wikipedia
graph is created). In either case, we have to check first if the graph is connected before we proceed
"""  
def get_graph_distance_metrics(G):
    H = None
    
    if G.is_directed():
        #if the graph is not strongly connected work with its undirected equivalent
        if nx.is_strongly_connected(G):
            H = G
        else:
            H = G.to_undirected()
    else:
        H = G
    
    #we can only calculate these metrics in a connected graph
    if nx.is_connected(H):
        
        radius = nx.radius(H)
        print("Radius = {}".format(radius))
        
        diam = nx.diameter(H)
        print("Diameter = {}".format(diam))
        
        center = nx.center(H)
        print("Center consists of:")
        for c in center:
            print(c)
            
        avg_path_length = nx.average_shortest_path_length(H)    
        print("Length of average shortest path: {}".format(avg_path_length))
    
    else:
        print("Cannot really calculate the graph's distance metrics because the graph is not connected")
    

"""
Function that calculates different clustering metrics of a given graph. More 
specifically, it calculates the number of triangles, the graph's transitivity and
its average clustering coefficient.
If the given graph is directed, in order to find the triangles in it we have to work
with its undirected equivalent
"""
def get_graph_clustering_metrics(G):   
    H = G.to_undirected()   
     
    #calculate the trangles in the equivalent undirected graph
    triangles_dict = nx.triangles(H)
    triangles = 0
    for x in triangles_dict.values():
        triangles += x
    triangles = triangles/3
        
    print("\nTotal triangles: {}".format(triangles))
    
    trans = nx.transitivity(G)
    print("Transitivity: {}".format(trans))
    
    avg_clustering = nx.average_clustering(G)
    print("Average clustering coefficient: {}".format(avg_clustering))
    

"""
Function that calculates the clustering coefficient of a specific node in the given graph.
"""    
def get_clustering_coeff(G):
    node_name = input("Name of the node you want to find its clustering coefficient: ")
    while node_name not in G:
        print("Node not found. Try again.")
        node_name = input("Name of the node you want to find its clustering coefficient: ")
        
    clust_coeff = nx.clustering(G,node_name)
    print("For node " + node_name + " the clustering coefficient is: {}".format(clust_coeff))
    

"""
Function that finds the nodes with the biggest degree centrality in a directed graph.
More specifically, it finds both the node with the biggest in-degree and the node with
the biggest out-degree.
If there are more than one nodes with the same maximum degree find all of them.
"""    
def degree_centrality(G):
    in_degree = list(G.in_degree())
    out_degree = list(G.out_degree())
    
    max_in_degree = -1
    max_in_node = list()
    
    for x in in_degree:
        if x[1] == max_in_degree:
            max_in_node.append(x[0])
        elif x[1] > max_in_degree:
            max_in_node.clear()
            max_in_degree = x[1]
            max_in_node.append(x[0])
            
    max_out_degree = -1
    max_out_node = list()
    
    for x in out_degree:
        if x[1] == max_out_degree:
            max_out_node.append(x[0])
        elif x[1] > max_out_degree:
            max_out_node.clear()
            max_out_degree = x[1]
            max_out_node.append(x[0])
    
    print("\nImportance of nodes based on their in-degree")
    print("The nodes with the biggest in-degree:")
    for x in max_in_node:
        print(x)
    print("with in-degree = {}".format(max_in_degree))
    
    print("\nImportance of nodes based on their out-degree")
    print("The nodes with the biggest out-degree:")
    for x in max_out_node:
        print(x)
    print("with out-degree = {}".format(max_out_degree))


"""
Function that finds the nodes with the biggest eigenvector centrality in a directed graph.
If there are more than one nodes with the same maximum eigenvector centrality find all of them.
"""    
def eigenvector_centrality(G):
    eigenvect_centr = nx.eigenvector_centrality(G)
    max_centr = -1
    max_node = list()
    
    for x,y in eigenvect_centr.items():
        if y == max_centr:
            max_node.append(x)
        elif y > max_centr:
            max_centr = y
            max_node.clear()
            max_node.append(x)
            
    print("\nImportance of nodes based on their eigenvenctor centrality")
    print("The nodes with the biggest eigenvector centrality:")
    for x in max_node:
        print(x)
    print("with eigenvector centrality value = {}".format(max_centr))


"""
Function that finds the nodes with the biggest katz centrality in a directed graph.
If there are more than one nodes with the same maximum katz centrality find all of them.
"""             
def katz_centrality(G):
    l_max = max(nx.adjacency_spectrum(G))   
    a = 1/(2*l_max) #epilexthike authaireta gia na isxuei oti a < 1/lmax
    
    katz_centr = nx.katz_centrality(G,alpha=a)
    max_centr = -1
    max_node = list()
    
    for x,y in katz_centr.items():
        if y == max_centr:
            max_node.append(x)
        elif y > max_centr:
            max_centr = y
            max_node.clear()
            max_node.append(x)
            
    print("\nImportance of nodes based on Katz centrality")
    print("The nodes with the biggest Katz centrality:")
    for x in max_node:
        print(x)
    print("with Katz centrality value = {}".format(max_centr))
  

"""
Function that finds the nodes with the biggest pagerank centrality in a directed graph.
If there are more than one nodes with the same maximum pagerank centrality find all of them.
"""    
def pagerank_centrality(G):
    pagerank_centr = nx.pagerank(G)
    max_centr = -1
    max_node = list()
    
    for x,y in pagerank_centr.items():
        if y == max_centr:
            max_node.append(x)
        elif y > max_centr:
            max_centr = y
            max_node.clear()
            max_node.append(x)
            
    print("\nImportance of nodes based on PageRank centrality")
    print("The nodes with the biggest PageRank centrality:")
    for x in max_node:
        print(x)
    print("with PageRank centrality value = {}".format(max_centr))
 

"""
Function that finds the hubs and the authorities in a directed graph.
If there are more than one nodes with the same hub or authority value find all of them.
"""   
def hits_centrality(G):
    hits_centr = nx.hits(G,max_iter=1000)
    
    hubs = hits_centr[0]
    authorities = hits_centr[1]
    
    #find hubs
    max_hub_val = -1
    max_hub_node = list()
    for x,y in hubs.items():
        if y == max_hub_val:
            max_hub_node.append(x)
        elif y > max_hub_val:
            max_hub_val = y
            max_hub_node.clear()
            max_hub_node.append(x)
    
    #find authorities
    max_auth_val = -1
    max_auth_node = list()
    for x,y in authorities.items():
        if y == max_auth_val:
            max_auth_node.append(x)
        elif y > max_auth_val:
            max_auth_val = y
            max_auth_node.clear()
            max_auth_node.append(x)
    
    print("\nImportance of nodes based on HITS")
    
    print("Most important hub(s):")
    for x in max_hub_node:
        print(x)
    print("with value {}".format(max_hub_val))
    
    print("Most important authorities:")
    for x in max_auth_node:
        print(x)
    print("with value {}".format(max_auth_val))  
    

"""
Function that finds the nodes with the biggest closeness centrality in a directed graph.
If there are more than one nodes with the same maximum closeness centrality find all of them.
"""    
def distance_centrality(G):
    dist_centr = nx.closeness_centrality(G)
    max_centr = -1
    max_node = list()
    
    for x,y in dist_centr.items():
        if y == max_centr:
            max_node.append(x)
        elif y > max_centr:
            max_centr = y
            max_node.clear()
            max_node.append(x)
            
    print("\nImportance of nodes based on closeness centrality")
    print("The nodes with the biggest closeness centrality:")
    for x in max_node:
        print(x)
    print("with closeness centrality value = {}".format(max_centr))
    

"""
Function that finds the nodes with the biggest betweenness centrality in a directed graph.
If there are more than one nodes with the same maximum betweenness centrality find all of them.
"""    
def position_centrality(G):
    pos_centr = nx.betweenness_centrality(G)
    max_centr = -1
    max_node = list()
    
    for x,y in pos_centr.items():
        if y == max_centr:
            max_node.append(x)
        elif y > max_centr:
            max_centr = y
            max_node.clear()
            max_node.append(x)
            
    print("\nImportance of nodes based on betweeness centrality")
    print("The nodes with the biggest betweeness centrality:")
    for x in max_node:
        print(x)
    print("with betweeness centrality value = {}".format(max_centr))  


"""
Function that creates a random G(n,l) graph. The created graph has the same
number of nodes (n) and edges (l) as the given graph.
"""    
def create_gnl_graph(G):
    l = G.number_of_edges()
    n = G.number_of_nodes()
    
    gnl_graph = nx.gnm_random_graph(n,l,directed=True)
    
    return gnl_graph


"""
Function that creates a graph based on the Watts-Strogatz model.
More specifically, the number of nodes is equal to that of the given graph,
each node is set to be connected with k=10 nearest neighbors in the ring topology
and the probability of rewiring each  edge is set to p = 0.01. 
The function also tries to create a connected graph so that it can be more easily used,
for example when calculating distance metrics of the graph.
"""
def create_small_world_graph(G):
    nnodes = G.number_of_nodes()
    k = 10
    p = 0.01
    ws_graph = nx.connected_watts_strogatz_graph(nnodes,k,p)
    
    return ws_graph
 
   
"""
Function that compares two graphs, the first one being the wikipedia graph we created
and the second one being either a random G(n,l) or a Watts-Strogatz model graph.
The two graphs are compared in terms of node degrees and distance and clustering metrics.
"""    
def compare_networks(G,G_new,graph_type):
    
    #compare average in and out degree
    print("\nFor the graph from wikipedia:")
    average_degree(G)
    print("For the " + graph_type + ":")
    average_degree(G_new)
    
    #compare biggest in and out degree
    print("\nFor the graph from wikipedia:")
    get_biggest_degree(G)
    print("\nFor the " + graph_type + ":")
    get_biggest_degree(G_new)
    
    #compare smallest in and out degree
    print("\nFor the graph from wikipedia:")
    get_smallest_degree(G)
    print("\nFor the " + graph_type + ":")
    get_smallest_degree(G_new)
    
    #compare degree distribution
    degree_distribution(G, "graph from wikipedia")
    degree_distribution(G_new, graph_type)
    
    #compare distance_metrics
    print("\nFor the graph from wikipedia:")
    get_graph_distance_metrics(G)
    print("\nFor the " + graph_type + ":")
    get_graph_distance_metrics(G_new)
    
    #for the random graph also calculate the theoretically expected
    #length of the average shortest path given by lnN/ln<k>
    if graph_type == "random graph":
        nnodes = G_new.number_of_nodes()  
        avg_degree = 0.0
        if G_new.is_directed():
            avg_degree = sum(d for n, d in G_new.in_degree())
        else:
            avg_degree = sum(d for n, d in G_new.degree())
        theor_avg_path_length = math.log(nnodes)/math.log(avg_degree)
        print("Theoretical length of average shortest path: {}".format(theor_avg_path_length))
    
    #compare clustering metrics
    print("\nFor the graph from wikipedia:")
    get_graph_clustering_metrics(G)
    print("\nFor the " + graph_type + ":")
    get_graph_clustering_metrics(G_new)
   

#create the graphs
nt = create_graph()
gnl_graph = create_gnl_graph(nt)
ws_graph = create_small_world_graph(nt)

#create .graphml files for the wikipedia graph, the random graph and the Watts-Strogatz model graph
nx.write_graphml(nt,"wikipedia_graph.graphml")
nx.write_graphml(gnl_graph,"random_graph.graphml")
nx.write_graphml(ws_graph,"small_world_graph.graphml")

flag = True

while flag:
    print("\nWhat do you want to do next?")
    print("Press 1 if you want to examine a specific node of the graph")
    print("Press 2 if you want to examine the properties of the whole graph")
    print("Press 3 if you want to compare the wikipedia graph with a random G(n,l) graph")
    print("Press 4 if you want to compare the wikipedia graph with a Watts-Strogatz model graph")
    print("Press 0 to exit")
    case = input("Enter a number depending on what you want to do:")
    
    if int(case) == 1:
        get_node_neighbor(nt)
        get_node_degree(nt)
        get_successors(nt)
        get_predecessors(nt)
        node_centrality(nt)
        get_eccentricity(nt)
        get_clustering_coeff(nt)
     
    elif int(case) == 2:
        get_nodes_edges(nt)
        get_biggest_degree(nt)
        get_smallest_degree(nt)
        degree_distribution(nt,"wikipedia graph")
        average_degree(nt)
        get_density(nt)
        check_components(nt)
        get_graph_distance_metrics(nt)
        get_graph_clustering_metrics(nt)
        degree_centrality(nt)
        eigenvector_centrality(nt)
        katz_centrality(nt)
        pagerank_centrality(nt)
        hits_centrality(nt)
        distance_centrality(nt)
        position_centrality(nt)

    elif int(case) == 3:
        compare_networks(nt,gnl_graph,"random graph")
        
    elif int(case) == 4:
        compare_networks(nt,ws_graph,"Watts-Strogatz model")

    elif int(case) == 0:
        flag = False
        
    else:
        print("Invalid input. Please try again.")


