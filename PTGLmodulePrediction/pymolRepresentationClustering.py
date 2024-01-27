########### functions ###########


def check_file_writable(fp):
    """Checks if the given filepath is writable"""
    if os.path.exists(fp):
        if os.path.isfile(fp):
            return os.access(fp, os.W_OK)
        else:
            return False
    # target does not exist, check perms on parent dir
    parent_dir = os.path.dirname(fp)
    if not parent_dir: parent_dir = '.'
    return os.access(parent_dir, os.W_OK)


########### simple classes ###########

class PTGLComplexGraphClusterDictionary:
	"""
	A class to prepare a dictionary, which holds the information about a PTGL complex graph and the corresponding clustering.
	Structure:
		key -> numbers of the nodes (here nodes represent chains of the structure):
		[chain ID: str, cluster membership: int, [chain ID contacts (edges): str], [edge weights corresponding to edges: str]]

	Attributes
	----------
	cluster_dictionary : dict, dictionary containing the PTGL graph information and the corresponding clustering

	Methods
	-------
	prepareDictionary():
		Parses the graph file input and the graph partitioning and prepares a dictionary for further usage.
	"""
	
	def __init__(self, graph_file_input: str, graph_partition: [int], vertexListOfGraph: list = [int], edge_weight_type: str = "absoluteWeight"):
		"""
		
		Initialization of PTGLComplexGraphClusterDictionary
		
		:param input_structure: str, the coordinate file of the structure (mmCIF or PDB)
		:param graph_file_input: str, the graph input file (only GML input supported)
		:param graph_partition: [[int]], list of vectors, which assigns the number of the community to a chain. 
                                        It is the result of a partition from iGraph.
        	:param vertexListOfGraph: list, list with the vertex IDs of the graph. This is
        							 necessary if the graph does not have all nodes in
        							 the partition, e.g. if a subgraph was used for 
        							 partitioning. For iGraph this list can be generated
        							 with [int(vertex['id']) for vertex in G.vs()]
		:param edge_weight_type: str, the edge weight, which should be read
			default: "absoluteWeight"
		
		"""
		
		self.cluster_dictionary: dict[str: str,int,int,int,[str],[str]] = {}
		self.graph_file: str = graph_file_input
		self.graph_partition: list = graph_partition
		self.edge_weight_type: str = edge_weight_type
		self.vsList: list = vertexListOfGraph
		self.isAaGraph: bool = False
		
		self.prepareDictionary()
		
		
	def prepareDictionary(self):
		"""
		
		Prepares a dictionary, which stores PTGL complex graph information and the corresponding partitioning.
		
		"""
		
		hit_node: bool = False
		hit_edge: bool = False
		in_block: bool = False
		node_id: int
		label: str
		fromResidue: int = None
		toResidue: int = None
		edge_source: str 
		edge_target: str
		edge_weight: str = None # In the case it is not available
		add_to_dictionary = True
		
		
		with open(self.graph_file) as f:
			for line in f.read().split("\n"):
				if not self.isAaGraph:
					tmp = line.split()
					if tmp and tmp[0] == "isAaGraph":
						self.isAaGraph = int(tmp[1])

				if not in_block:
					# in graph header
					if "node [" in line:
						hit_node = True
						in_block = True
					elif "edge [" in line:
						hit_edge = True
						in_block = True
						# once all the nodes are read, we go through the edges and have to add two lists to the dictionary first
						if add_to_dictionary:
							for item in self.cluster_dictionary:
								self.cluster_dictionary[item].append([])
								self.cluster_dictionary[item].append([])
							add_to_dictionary = False
					else:
						pass
				else:
					# first go through all nodes
					if hit_node:
						if line.split()[0] == "id":
							node_id = str(line.split()[-1])
						elif line.split()[0] == "label":
							label = "".join(line.split()[1:]).replace("\"", "") # In case of quotes inside the string, maybe some better solution required.
						elif line.split()[0] == "pdbResStart":
							fromResidue = str(line.split()[-2].split('-')[1])
						elif line.split()[0] == "pdbResEnd":
							toResidue = str(line.split()[-2].split('-')[1])
						elif line.split()[0] == "]":
							in_block = False
							hit_node = False
							node_in_cluster = False
							if self.isAaGraph:
								fromResidue =str(label.split('-')[2])
								toResidue = str(label.split('-')[2])
							for number, item in enumerate(self.graph_partition):
								for entry in item:
									if node_id == str(self.vsList[entry]):
										self.cluster_dictionary[node_id] = [label, number, fromResidue, toResidue]
										node_in_cluster = True
							if not node_in_cluster:
								self.cluster_dictionary[node_id] = [label, -1, fromResidue, toResidue]
							
							fromResidue = None
							toResidue = None
							pass
						else: pass

					# then go through all edges            
					elif hit_edge:
						if line.split()[0] == "source":
							edge_source = str(line.split()[-1])
						elif line.split()[0] == "target":
							edge_target = str(line.split()[-1])
						elif line.split()[0] == f"{self.edge_weight_type}":
							edge_weight = str(line.split()[-1])
						elif line.split()[0] == "]":
							in_block = False
							hite_edge = False
							self.cluster_dictionary[edge_source][4].append(self.cluster_dictionary[edge_target][0])
							self.cluster_dictionary[edge_source][5].append(edge_weight)
							self.cluster_dictionary[edge_target][4].append(self.cluster_dictionary[edge_source][0])
							self.cluster_dictionary[edge_target][5].append(edge_weight)
							pass
						else: pass

					else:
				    		pass
        


########### vamos ###########

def preparePML(ClusterDictionary: dict, clusterColors: list, input_structure: str, isAaGraph: bool):
	""" 
	:param clusterColor: list, the list with the color values of the clusters. For iGraph
								this list can be obtained by the following command using
								the length of the partition:
								[list(ig.color_name_to_rgb(color)) for color in ig.drawing.colors.ClusterColoringPalette(len(partition))]
	"""

	# prepare variables for the PML script
	number_of_clusters: int = 0
	clustered_chains: dict = dict()    
	pseudoatoms: str = ""
	group_distances: str = ""
	distance_calculations: str = ""
	pseudoatoms_for_label: str = ""
	close_contact_groups: str = ""
	SSEgraph: bool = False
		

	for key, nested in ClusterDictionary.items():
		if not SSEgraph and not isAaGraph:
			if nested[2] is not None and nested[3] is not None:
				SSEgraph = True
		if nested[1] > number_of_clusters:
			number_of_clusters = nested[1]
		if str(nested[1]) in clustered_chains:
			clustered_chains[str(nested[1])].append(nested[0])
		else:
			clustered_chains[str(nested[1])] = [nested[0]]
		if isAaGraph:
			pseudoatoms += f"pseudoatom ps{nested[0][1:]}, chain {nested[0].split('-')[1]} and residue {nested[2]}\n"
		elif SSEgraph:
			pseudoatoms += f"pseudoatom ps{nested[0]}, chain {nested[0].split('-')[0]} and residue {nested[2]}-{nested[3]}\n"
		else:
			pseudoatoms += f"pseudoatom ps{nested[0]}, chain {nested[0]}\n"

		for position, item in enumerate(nested[4]):
			distance_calculations += f"dist dist_{nested[0]}-{item}, ps{nested[0]}////PS1, ps{item}////PS1\nhide label, dist_{nested[0]}-{item}\n"
			pseudoatoms_for_label += f"pseudoatom ps_{nested[0]}-{item}, ps{nested[0]} or ps{item}, label={nested[5][position]}\n"

		group_distances += f"group chain-{nested[0]}_contacts, dist_{nested[0]}-*\n"
		close_contact_groups += f"group chain-{nested[0]}_contacts, close\n"

	# Prepare coloring of chains according to cluster membership
	# First get number of clusters

	filePML: str
	select_clusters: str = ""
	color_cluster: str = ""

	for cluster in range(number_of_clusters+1):
		if isAaGraph:
			#select_clusters += "select cluster{}, {}\n".format(cluster, "".join(str(" or " + "chain " + e.split('-')[1] + " and residue " + ClusterDictionary[e.split('-')[0][1:]][2]) for e in clustered_chains[str(cluster)])[3:])
			select_clusters += "select cluster{}, {}\n".format(cluster, "".join(str(" or " + "chain " + e.split('-')[1] + " and residue " + ClusterDictionary[e.split('-')[0]][2]) for e in clustered_chains[str(cluster)])[3:])
		elif SSEgraph:
			select_clusters += "select cluster{}, {}\n".format(cluster, "".join(str(" or " + "chain " + e.split('-')[0] + " and residue " + ClusterDictionary[e.split('-')[1]][2] + "-" + ClusterDictionary[e.split('-')[1]][3]) for e in clustered_chains[str(cluster)])[3:])
		else:
			select_clusters += "select cluster{}, chain {}\n".format(cluster, "".join(str("+" + e) for e in clustered_chains[str(cluster)])[1:])
		color_cluster += f"set_color newColor{cluster}, {clusterColors[cluster]} \ncolor newColor{cluster}, cluster{cluster}\n"
	if any([key == "-1" for key in clustered_chains.keys()]):
		if isAaGraph:
			#select_clusters += "select unclustered, chain {}\n".format("".join(str(" or " + "chain " + e.split('-')[1] + " and residue " + ClusterDictionary[e.split('-')[0][1:]][2]) for e in clustered_chains["-1"])[3:])
			select_clusters += "select unclustered, chain {}\n".format("".join(str(" or " + "chain " + e.split('-')[1] + " and residue " + ClusterDictionary[e.split('-')[0]][2]) for e in clustered_chains["-1"])[3:])
		elif SSEgraph:
			select_clusters += "select unclustered, chain {}\n".format("".join(str(" or " + "chain " + e.split('-')[0] + " and residue " + ClusterDictionary[e.split('-')[1]][2] + "-" + ClusterDictionary[e.split('-')[1]][3]) for e in clustered_chains["-1"])[3:])
		else:
			select_clusters += "select unclustered, chain {}\n".format("".join(str("+" + e) for e in clustered_chains["-1"])[1:])
		color_cluster += f"color white, unclustered\n"

	# Prepare the PML script for PyMol
	filePML = f"""
### PML file to represent the clusters in the structure.
# Python is by default case insensitive for chain IDs but for large clusters with upper and lower case letters we have to care!
set ignore_case, off

# Load the structure
load {input_structure}, complex

# First color everything white
color white, all
# Select and color the found clusters
{select_clusters}
{color_cluster}
# Create pseudoatoms at the geometric center of the chains
{pseudoatoms}
group pseudoatoms, ps*

# Create all distance objects originating from each chain (reverse dublicates are needed!)
# Hide the distance label, generate a pseudoatom between both pseudoatoms for which the distance was calculated
# and label it according to the edge weight
{distance_calculations}
# Group all distance measurements originating from one chain to all other connected ones
{group_distances}
# Disallow addition/removal from the group (allow with 'open')
{close_contact_groups}
# Set pseudoatoms along the distance measurement line and label them with the corresponding edge weights
{pseudoatoms_for_label}
group label_pseudoatoms, ps_*

# Set some parameters
set dash_gap, 0.7
set dash_radius, 0.2
set dash_color, black
set label_position,(1,1,1) # label offset
set label_color, black
bg_color white
set cartoon_highlight_color, grey40

## If coordinates of the pseudoatoms should be written
#xyz = cmd.get_coords("pseudoatoms", 1)
#python
#import numpy as np
#with open('pseudoatomCoords.npy', 'wb') as f:
#    np.save(f, xyz)
#python end

	"""

	return filePML
	
def prepareClusteringRepresentationPyMOL(protein_structure, graph_file_path, graph_partition, vertexListOfGraph, clusterColors):
	# prepare dictionary with PTGL complex graph and partitioning information
	ClusterDictionary = PTGLComplexGraphClusterDictionary(graph_file_path, graph_partition, vertexListOfGraph)
	return preparePML(ClusterDictionary.cluster_dictionary, clusterColors, protein_structure, ClusterDictionary.isAaGraph)

	
	
	

