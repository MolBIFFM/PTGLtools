/* File: GraphService.cpp
 * Author: ben
 * 
 * Created on April 29, 2015 3:49 PM
 */


#include "GraphService.h"
#include <unordered_map>
#include <utility>

using namespace boost;

/* This Service Class is meant for graphs which were parsed from a GML file.
 * GML files are of the following form.
 * 
 * x    is any arbitrary integer
 * s    is a string
 * att  represents an arbitrary number of other attributes which are written in
 *      the same form as id, label, comment and directed
 * 
 * graph [
 * comment s
 * id x
 * label s
 * directed x (x=0 => undirected; x=1 => directed)
 * att
 * ]
 * 
 * node [
 * comment s
 * id x
 * label s
 * att
 * ]
 * 
 * For edges, the values for source and target refer to the id of the
 * corresponding nodes
 * 
 * edge [
 * source x
 * target x
 * comment s
 * label s
 * att
 * ]
 * 
 */

GraphService::GraphService() {
    Graph g_tmp(1);
    g = g_tmp;
    gc = GraphletCounts(g);
    graphlet_identifier = "";
    graphlet_patterns = std::vector<std::string>();
}

GraphService::GraphService(const Graph& graph) {
    g = graph;
    gc = GraphletCounts(g);
    graphlet_identifier = "";
    graphlet_patterns = std::vector<std::string>();
}


/* Returns the graph */
Graph GraphService::getGraph() const {
  return g;  
}

/* Return the keys i.e. the property names of a given graph
 * @param graph g
 * @return vector<string> properties */
std::vector<std::string> GraphService::getGraphProperties() {
    std::unordered_map<std::string, std::string> propertyMap = g[graph_bundle].properties;
    
    std::vector<std::string> keys = std::vector<std::string>();
    keys.reserve(propertyMap.size());

    for(auto kv : propertyMap) {
        keys.push_back(kv.first);
    } 
    return keys;   
};


    
    
/* Return the value of a given property as a string.
 * @param1 graph g - the graph
 * @param2 string prop - the property
 * @return string property */
std::string GraphService::getPropertyValue(std::string prop) {
    
    std::string propValue = g[graph_bundle].properties[prop];
    return propValue;
}; 
    

    
/* Return the number of vertices
 * @param Graph g
 * @return int number of edges */    
int GraphService::getNumVertices() {
    return num_vertices(g);
};

/* Return the number of edges
 * @param Graph g
 * @return int number of edges */
int GraphService::getNumEdges() {
    return num_edges(g);
};



/* Return a vector, containing the vertex indices for a given graph
 * @return vector<int> vertices */
std::vector<int> GraphService::getVertices() {


    std::vector<int> vertexVector;
    vertexVector.reserve(num_vertices(g));   
    typedef property_map<Graph,vertex_index_t>::type IndexMap;
    IndexMap index = get(vertex_index, g);
    std::pair<VertexIterator, VertexIterator> vertexPair;
    for (vertexPair = vertices(g); vertexPair.first != vertexPair.second; ++vertexPair.first)
        vertexVector.push_back(index[*vertexPair.first]);
    return vertexVector;
};


/* Return the a vector, containing the edges as integer pairs.
 * @return vector<pair<int, int>> edgeVector -- a vector containing the edges */
std::vector<std::pair<int, int>> GraphService::getEdges() {
    
    // create the vector and the pair
    std::vector<std::pair<int, int>> edgeVector;
    edgeVector.reserve(num_edges(g));
    std::pair<int, int> edgePair;
    
    // create a map, so the graph properties can be read
    typedef property_map<Graph, vertex_index_t>::type IndexMap;
    IndexMap index = get(vertex_index, g);
    
    // iterate over the edges in the map and save them to the vector
    EdgeIterator ei, eiEnd;
    for (tie(ei, eiEnd) = edges(g); ei != eiEnd; ++ ei) {
        edgePair.first = index[source(*ei, g)];
        edgePair.second = index[target(*ei, g)];
        edgeVector.push_back(edgePair);
    };
    return edgeVector;
    
};







/* Compute the node degree distribution
 * @return vector<int> degDist -- the node degree distribution */
std::vector<int> GraphService::computeDegreeDist() {
    
    // create the vector and a variable for the degree
    int degree;
    std::vector<int> degDist (num_vertices(g));
    
    
    
    // iterate over the vertices
    VertexIterator vi, vi_end, next;
    tie(vi, vi_end) = vertices(g);
    int i = 0;
    
    
    for (next = vi; vi != vi_end; vi = next) {
        ++next;
        degree = out_degree(*vi, g); // store their degree
        degDist[i] = degree; // and put it into the vector
        i++;
        
    }
    return degDist;
  
};

/* Returns the graph's label
 * @return string name -- the graph's name */
std::string GraphService::get_label() {
    std::string name = g[graph_bundle].label;
    return name;
    
};

/* Returns a vector of vertices adjacent to the vertex with the given id i
 * @param <int> the vertex id
 * @return <vector<int>> all vertices adjacent to i */
std::vector<int> GraphService::get_adjacent(int i) {
    AdjacencyIterator first, last;
    std::vector<int> vertex_vector = std::vector<int>();
    vertex_vector.push_back(i);
    
    
    
    
    for (tie(first,last) = adjacent_vertices(i,g); first != last; ++first) {
                
        vertex_vector.push_back(g[*first].id);
    }
        
    return vertex_vector;
}


/* Returns a nested vector containing vector which store the id of each 
 * vertex adjacent to the vertex with the id that corresponds to the position
 * in the vector
 * @return <vector<vector<int>>  */
std::vector<std::vector<int>> GraphService::get_adjacent_all() {
    
    std::vector<std::vector<int>> adj_all_vector = std::vector<std::vector<int>>();
    
    int n = getNumVertices();
    
    for (int i = 0; i < n; i++) {
        std::vector<int> adj_vector = get_adjacent(i);
        
        adj_all_vector.push_back(adj_vector);
    }
    
    return adj_all_vector;
    
}

/* Returns a nested vector containing the absloute counts of all graphlets in
 * graph g
 * @return <vector<vector<int>> the counts of all graphlets */
std::vector<std::vector<int>> GraphService::get_abs_counts() {
    std::vector<std::vector<int>> out_vec = gc.get_unlabeled_abs_counts();
    
    return out_vec;
}


/* Returns a nested vector containing the normalized counts of all graphlets in
 * graph g
 * @return <vector<vector<float>> */
std::vector<std::vector<float>> GraphService::get_norm_counts() {
    std::vector<std::vector<float>> out_vec = gc.get_normalized_counts();
    
    return out_vec;
}


/* Get a unordered_map containing the patterns spcified in the argument patterns as keys
 * and their corresponding counts as values.
 * @param <string> id - the property of the vertex info for which grahlet patterns should be counted
 * @param <vector<string>> - the patterns which should be counted
 * @return <unorded_map<string, vector<int>>> - a map where the absolute labeled counts can be accessed by their pattern */
std::unordered_map<std::string, std::vector<int>> GraphService::get_labeled_abs_counts(std::string id, std::vector<std::string> patterns) {
    //TODO: finish support for returning graphlet counts map
    
    // define the map to be filled
    std::unordered_map<std::string,std::vector<int>> map = std::unordered_map<std::string, std::vector<int>>();
    
    // copying graphlet identifier to field
    graphlet_identifier = id;
    
    
    //setting total counts
    total_counts = 0;
    
    //initialising vectors to store data
    std::vector<int> abs_counts2;
    std::vector<std::vector<int>> abs_counts3;
    std::vector<std::vector<int>> abs_counts4;
    
    std::vector<std::string> pat_vec2 = std::vector<std::string>();
    std::vector<std::string> pat_vec3 = std::vector<std::string>();
    std::vector<std::string> pat_vec4 = std::vector<std::string>();
    
    std::vector<std::vector<std::string>> pat_vecvec3 = std::vector<std::vector<std::string>>(2);
    std::vector<std::vector<std::string>> pat_vecvec4 = std::vector<std::vector<std::string>>();
    
    
    // separating patterns by size and copying them into vectors
    for (int i = 0; i < patterns.size(); i++) {
        
        std::string str = patterns[i];
        
        if (str.length() == 2) {
            pat_vec2.push_back(str);
            
        } else if (str.length() == 3) {
            pat_vec3.push_back(str);
            
            
        } else if (str.length() == 4) {
            pat_vec4.push_back(str);
            
            
        } else {
            std::cerr << "WARNING: Trying to look for graphlet of unsupported size. " << std::endl
                    << "Only labeled graphlets of with 2,3 or 4 vertices are supported." <<std::endl;
        }
        
    }
    
    pat_vecvec3[0] = pat_vec3;
    pat_vecvec3[1] = pat_vec3;
    pat_vecvec4[5] = pat_vec4;
    
    // computing lableled graphlet counts
    abs_counts2 = gc.get_labeled_2_countsABS(graphlet_identifier, pat_vec2);
    abs_counts3 = gc.get_labeled_3_countsABS(graphlet_identifier, pat_vecvec3);
    abs_counts4 = gc.get_labeled_4_countsABS(graphlet_identifier, pat_vecvec4);
    

    // copying graphlet counts into map
    for (int i = 0; i< pat_vec2.size(); i++) {
        map[pat_vec2[i]] = abs_counts2;
        
    }
    for (int i = 0; i< pat_vec3.size(); i++) {

            
            std::vector<int> counts_for_1_label = std::vector<int>();
            counts_for_1_label.push_back(abs_counts3[0][i]);
            counts_for_1_label.push_back(abs_counts3[1][i]);
                    
            map[pat_vec3[i]] = counts_for_1_label;
        
    }
    for (int i = 0; i < pat_vecvec4[5].size(); i++) {
        
        std::vector<int> vec5graphlets = std::vector<int>();
        vec5graphlets.push_back(abs_counts4[5][i]);
        
        map[pat_vecvec4[5][i]] = vec5graphlets;
    }
    
    
    
    return map;
}


std::unordered_map<std::string, std::vector<float>> GraphService::get_labeled_norm_counts(std::string id, std::vector<std::string> patterns) {
    std::unordered_map<std::string,std::vector<float>> map = std::unordered_map<std::string, std::vector<float>>();
    
    std::unordered_map<std::string, std::vector<int>> imap = get_labeled_abs_counts(id, patterns);
    
    
    graphlet_identifier = id;
    
    
    return map;
}

std::string GraphService::get_graphlet_identifier() {
    return graphlet_identifier;
}

std::vector<std::string> GraphService::get_patterns() {
    
}

GraphService & GraphService::operator =(const GraphService & serv) {
    
    if (this == &serv) {
        return *this;
    }
    //copy_graph(serv.getGraph(),g);
    
    g = serv.getGraph();
    gc = GraphletCounts(g);
    graphlet_identifier = "";
    graphlet_patterns = std::vector<std::string>();
    return *this;
}