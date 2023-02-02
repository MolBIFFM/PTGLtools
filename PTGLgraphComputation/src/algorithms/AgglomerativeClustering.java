/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Jan Niclas Wolf 2021. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author jnw
 */
package algorithms;

import datastructures.ClusteringResult;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.IntStream;
import java.util.stream.Collectors;
import java.util.Scanner;
import java.util.Random; 
import static java.util.stream.Collectors.toList;
import proteingraphs.ComplexGraphEdgeWeightTypes.EdgeWeightType;
import proteingraphs.ComplexGraphEdgeWeightTypes;
import proteinstructure.Atom;
import proteinstructure.Chain;
import static proteinstructure.Chain.AtomWeights;
import proteinstructure.Molecule;
import tools.DP;
import settings.Settings;

/**
 *
 * @author jnw
 * modified by Christian Ickes
 * modified by Florian Gisdon
 */
public class AgglomerativeClustering {
    private final static String CLASS_TAG = "CG";
    private static final Random RANDOM = new Random(); 
    
    private ArrayList<Edge> edges;
    private final Map<Integer, Double> chainProperties;
    private final EdgeWeightType weightType;
    private final Map<Integer, String> labelMap;
    private final List<Chain> chains;
    private Map<String, List<String>> mergedChains;
    private int curStepNum;
    
    class Edge implements Comparable<Edge> {
        int v1;  // always lower ID
        int v2;  // always greater ID
        private int absoluteWeight;
        private BigDecimal normalizedWeight;
        private final EdgeWeightType normalizationType;
        
        Edge (int v1, int v2, int absoluteWeight, EdgeWeightType normalizationType, BigDecimal normalizedWeight) {
            this.v1 = Math.min(v1, v2);
            this.v2 = Math.max(v1, v2);
            this.absoluteWeight = absoluteWeight;
            this.normalizationType = normalizationType;
            this.normalizedWeight = normalizedWeight;            
        }
        
        public Edge mergeEdgeIntoThis(Edge otherEdge, double propertyChainA, double propertyChainB) {
            if (this.normalizationType != otherEdge.normalizationType) {
                DP.getInstance().e(CLASS_TAG, "Cannot merge edges of different normalization. Returning first edge and trying to go on. "
                        + "Please report this bug to the developers.");
            }
            absoluteWeight += otherEdge.absoluteWeight;
            normalizedWeight = ComplexGraphEdgeWeightTypes.computeLengthNormalization(absoluteWeight, propertyChainA, propertyChainB, normalizationType);

            return this;
        }
        
        @Override
        public int compareTo(Edge otherEdge) {
            return normalizedWeight.compareTo(otherEdge.normalizedWeight);
        }
        
        @Override
        public int hashCode() {
            return Objects.hash(v1, v2);
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj) {
                return true;
            }
            if (obj == null) {
                return false;
            }
            if (getClass() != obj.getClass()) {
                return false;
            }
            final Edge other = (Edge) obj;
            if (this.v1 != other.v1) {
                return false;
            }
            return this.v2 == other.v2;
        }
        
        @Override
        public String toString() {
            return "Edge {" + v1 + "," + v2 + "}, absolute weight: " + absoluteWeight + " normalization: " + normalizationType.name + ", normalized weight: " + normalizedWeight.toString();
        }
        
        public int[] getVertices() {
            int[] vertices = {v1, v2};
            return vertices;
        }
        
        public int getAbsoluteWeight() {
            return absoluteWeight;
        }
    }
    
    class EdgeComparatorNormalizedWeight implements Comparator<Edge> {
        @Override
        public int compare(Edge e1, Edge e2) {
            return e1.compareTo(e2);
        }
    }
    
    /**
     * Creates a new instance of the clustering to copy it completely. 
     * @param aggloClust
     * @return complete new instance 
     */
    public static AgglomerativeClustering newInstance(AgglomerativeClustering aggloClust){
        return new AgglomerativeClustering(aggloClust);
    }
    
    /**
     * Constructor of an agglomerative clustering.Class cannot be static so that it can hold instances of nested class Edge.
     * @param edgeList Array of [vertex ID 1, vertex ID 2, edge weight]
     * @param chainProperties
     * @param weightType
     * @param labelMap
     * @param chains
     */   
    public AgglomerativeClustering(int[][] edgeList, Map<Integer, Double> chainProperties, EdgeWeightType weightType, Map<Integer, String> labelMap, List<Chain> chains) {
        this.curStepNum = 0;
        
        this.edges = new ArrayList<>();
        for (int[] edge : edgeList) {
            if (edge.length < 3) {
                DP.getInstance().e(CLASS_TAG, "Got an edge without two particpating vertices. Skipping edge and trying to go on. "
                        + "Please report this bug to the developers.");
                continue;
            }
            
            BigDecimal normalizedWeight;
            
            normalizedWeight = ComplexGraphEdgeWeightTypes.computeLengthNormalization(edge[2], chainProperties.get(edge[0]), chainProperties.get(edge[1]), 
                    weightType);

            edges.add(new Edge(edge[0], edge[1], edge[2], weightType, normalizedWeight));
        }
        
        this.chainProperties = new HashMap<>(chainProperties);  // Shallow copy to not change chain lengths outside of this class
        this.weightType = weightType;
        this.labelMap = labelMap;
        this.chains = chains;
        this.mergedChains = new HashMap<>();
        for (HashMap.Entry<Integer,String> entry : labelMap.entrySet()) {
            mergedChains.put(entry.getValue(), Arrays.asList(entry.getValue()));
        }
    }
    
    /**
     * Used to create copy by newInstance
     * @param aggloClust 
     */
    public AgglomerativeClustering(AgglomerativeClustering aggloClust){
        this.edges = new ArrayList<>(aggloClust.edges);
        this.chainProperties = new HashMap<>(aggloClust.chainProperties);
        this.weightType = aggloClust.weightType;
        this.labelMap = new HashMap<>(aggloClust.labelMap);
        this.curStepNum = 0;
        this.chains = aggloClust.chains;
        this.mergedChains = new HashMap<>(aggloClust.mergedChains);
    }
            
    /**
     * Performs a hierarchical clustering algorithm on the given complex graph. 
     * @return Clustered Tree
     */
    public ClusteringResult chainInformationClustering() {
        ClusteringResult clusteringResult = new ClusteringResult(true);
    
        while (edges.size() > 0) {
            chainInformationClusteringStep(clusteringResult);
            
            if (Settings.getInteger("PTGLgraphComputation_I_debug_level") > 1) {
                System.out.println("[DEBUG LV 2] Merges after step " + curStepNum);
                System.out.println(clusteringResult.toString());
            }
        }
   
        return clusteringResult;
    }
    
    private void chainInformationClusteringStep(ClusteringResult clusteringResult) {
        // 1) sort by normalized weight
        Collections.sort(edges, Collections.reverseOrder(new EdgeComparatorNormalizedWeight()));
        
        if (Settings.getInteger("PTGLgraphComputation_I_debug_level") > 2) {
            System.out.println("[DEBUG LV 3] Sorted edge list from step " + curStepNum);
            System.out.println("[DEBUG LV 3] Assembly type: " + 
                                Settings.getInteger("PTGLgraphComputation_I_type_assembly_prediction"));
            edges.forEach(edge -> {
                System.out.println("  " + edge.toString());
            });
            
        }
        
        if (Settings.getInteger("PTGLgraphComputation_I_debug_level") > 1) {
            System.out.println("[DEBUG LV 2] Edge with highest normalized weight: " + edges.get(0).toString());
        }
        
        // choose edge to merge: greedy or interactive
        int mergeEdgeIndex;
        switch(Settings.getInteger("PTGLgraphComputation_I_type_assembly_prediction")){
            case(1): mergeEdgeIndex = getEdgeIndexFromUser(6, clusteringResult); break;
            case(2): mergeEdgeIndex = getEdgeIndexStochastic(); break;
            case(3): mergeEdgeIndex = getEdgeIndexRandom(); break;
            default: mergeEdgeIndex = 0; break;
        }
        
        // add vertices of edges to merges of cluster result
        Edge curEdge = edges.get(mergeEdgeIndex);
        clusteringResult.addMerge(
                IntStream.of(curEdge.getVertices()).boxed().toArray(Integer[]::new), 
                curEdge.absoluteWeight,
                chainProperties.get(curEdge.v1),
                chainProperties.get(curEdge.v2));
        
        int v1 = curEdge.v1;
        int v2 = curEdge.v2;
        
        // 2.1) remove edge between merged vertices
        edges.remove(mergeEdgeIndex);
        
        // 2.2) update chain lengths or chain gyradii
        if(weightType.shortName.equals("addGyrInd") || weightType.shortName.equals("addGyr")){
            // update the mergedChains List: mergedChains.get(labelMap.get(v2)) should be added to the list of mergedChains.get(labelMap.get(v1))
            // and the entry mergedChains.get(labelMap.get(v2)) should be deleted
            ArrayList<String> merge = new ArrayList<>();
            merge.addAll(mergedChains.get(labelMap.get(v1)));
            merge.addAll(mergedChains.get(labelMap.get(v2)));
            mergedChains.put(labelMap.get(v1), merge);
            mergedChains.remove(labelMap.get(v2));
            // do the recalculation with the new cluster
            //System.out.println("Merged chains list: " + mergedChains.get(labelMap.get(v1)));
            double recalculatedChainGyradius = recomputeChainRadiusOfGyration("mass", mergedChains.get(labelMap.get(v1))); 
            chainProperties.replace(v1, recalculatedChainGyradius);  // v1 is representative
            chainProperties.remove(v2);
        } else {
            double combinedChainLength = chainProperties.get(v1) + chainProperties.get(v2);
            chainProperties.replace(v1, combinedChainLength);  // v1 is representative
            chainProperties.remove(v2);
        }
     
        // 2.3 update edges
        HashMap<Edge, Edge> edgeMap = new HashMap<>();
        edges.forEach(edge -> {
            Edge updatedEdge = edge;  // Set it to edge and update below if necessary
            if (IntStream.of(edge.getVertices()).anyMatch(v -> v == v2 || v == v1)) {
                // possibly update indices
                int updatedV1 = (edge.getVertices()[0] == v2 ? v1 : edge.getVertices()[0]);
                int updatedV2 = (edge.getVertices()[1] == v2 ? v1 : edge.getVertices()[1]);
                int absoluteWeight = edge.getAbsoluteWeight();

                BigDecimal updatedNomalizedWeight = ComplexGraphEdgeWeightTypes.computeLengthNormalization(absoluteWeight, chainProperties.get(updatedV1), chainProperties.get(updatedV2), 
                    weightType);
                updatedEdge = new Edge(updatedV1, updatedV2, absoluteWeight, weightType, updatedNomalizedWeight);
            } 
            
            // try to add updated edge
            if (edgeMap.containsKey(updatedEdge)) {
                // duplicate edge -> merge the two edges
                    edgeMap.get(updatedEdge).mergeEdgeIntoThis(updatedEdge, chainProperties.get(updatedEdge.v1), chainProperties.get(updatedEdge.v2));

            } else {
                edgeMap.put(updatedEdge, updatedEdge);
            }
        });
        
        edges = new ArrayList<>(edgeMap.values());    
        curStepNum++;
    }
    
    private double recomputeChainRadiusOfGyration(String method, List<String> mergedChainsList) {
        // The radius of gyration can be calculated with mass weighting or according to the geometrical center

        Double chainRadiusOfGyration = null;
        if ("mass".equals(method)) {
            double[] tmpCenter = new double[3];
            double[] centerOfMass = new double[3];
            tmpCenter[0] = tmpCenter[1] = tmpCenter[2] = 0.0;
            double tmpSumAtomWeights = 0.0;
            for (Chain currChainL1 : chains){
                //System.out.println("Merged chains list: " + mergedChainsList);
                //System.out.println("Current Chain ID: " + currChainL1.getPdbChainID());
                // look for chain matches with the mergedChainsList
                if (mergedChainsList.stream().anyMatch(s -> s.contains(currChainL1.getPdbChainID()))) {
                    //System.out.println("Match! - Chain ID: " + currChainL1.getPdbChainID());
                    // compute center of mass and the total mass for all chains
                    for (Molecule mol : currChainL1.getMolecules()) {
                        if(Settings.getBoolean("PTGLgraphComputation_B_CG_ignore_ligands")){
                            if(!mol.isLigand()){
                                for (Atom a : mol.getAtoms()) {
                                    //System.out.println("Begin current Atom: " + a);
                                    tmpCenter[0] += a.getCoordX()*AtomWeights.get(a.getChemSym());
                                    tmpCenter[1] += a.getCoordY()*AtomWeights.get(a.getChemSym());
                                    tmpCenter[2] += a.getCoordZ()*AtomWeights.get(a.getChemSym());
                                    tmpSumAtomWeights += AtomWeights.get(a.getChemSym());
                                }
                            }
                        } else {
                            for (Atom a : mol.getAtoms()) {
                                //System.out.println("Begin current Atom: " + a);
                                tmpCenter[0] += a.getCoordX()*AtomWeights.get(a.getChemSym());
                                tmpCenter[1] += a.getCoordY()*AtomWeights.get(a.getChemSym());
                                tmpCenter[2] += a.getCoordZ()*AtomWeights.get(a.getChemSym());
                                tmpSumAtomWeights += AtomWeights.get(a.getChemSym());
                            }
                        }
                    }
                }
            }
            
            // there could be the case that only RNA/DNA atoms are in a chain
            // in this case radius is set to -1
            if (tmpSumAtomWeights > 0) {
                centerOfMass[0] =  tmpCenter[0] / tmpSumAtomWeights;
                centerOfMass[1] =  tmpCenter[1] / tmpSumAtomWeights;
                centerOfMass[2] =  tmpCenter[2] / tmpSumAtomWeights;

                if (Settings.getInteger("PTGLgraphComputation_I_debug_level") > 0) {
                    System.out.println("[DEBUG] Center of mass of chains " + mergedChainsList + " is at " + Arrays.toString(centerOfMass));
                }
                // compute the radius of gyration
                double tmpGyradius;
                tmpGyradius = 0.0;
                
                // And here comes the same loop again since now we have the center of mass
                for (Chain currChainL2 : chains){
                    if (mergedChainsList.stream().anyMatch(s -> s.contains(currChainL2.getPdbChainID()))) {
                        for (Molecule mol : currChainL2.getMolecules()) {
                            for (Atom a : mol.getAtoms()) {
                                tmpGyradius +=  AtomWeights.get(a.getChemSym())*Math.pow(Math.sqrt(Math.pow((a.getCoordX()-centerOfMass[0]),2) + Math.pow((a.getCoordY()-centerOfMass[1]),2) + Math.pow((a.getCoordZ()-centerOfMass[2]),2)), 2);
                            }
                        }
                    }
                }
                // calculate chainRadiusOfGyration
                chainRadiusOfGyration = Math.sqrt(tmpGyradius / tmpSumAtomWeights);
            } else {
                DP.getInstance().w("The chains " + mergedChainsList + " seem not to hold protein atoms. No center can be detected.", 2);
                chainRadiusOfGyration =  -1.0;
            }  
        }
        //System.out.println("Gyradius: " + chainRadiusOfGyration);
        return chainRadiusOfGyration;
    }
    // TODO: implement Gyradius recalculation with geometrical center and implement this option in the settings!
 
    /**
     * Print a table of all edges and ask user for an input until it is a valid edge index.
     * @return valid edge index as int
     */
    private int getEdgeIndexFromUser(int numberIndentSpaces, ClusteringResult clusteringResult) {
        Integer edgeIndex = null;
        Scanner lineReader = new Scanner(System.in);
        System.out.println(" ".repeat(numberIndentSpaces) + "Step " + curStepNum + ": Please choose an edge for merging");
        for (String line : prettyFormatEdgeList(numberIndentSpaces + 2, clusteringResult).split("\n")) {
            System.out.println(line);
        }
        while (edgeIndex == null) {
            String input = lineReader.nextLine();
            try {
                edgeIndex = Integer.parseInt(input.trim());
                if (edgeIndex < 0 || edgeIndex > edges.size() - 1) {
                    System.out.println("Please choose a number between zero and " + (edges.size() - 1));
                    edgeIndex = null;
                }
            } catch (NumberFormatException e) {
                System.out.println("Please provide a valid number");
            }
        }
        return edgeIndex;
    }
    
    /**
     * Chooses an Edge by performing a roulette algorithm.
     * Probabilities are equal to edge weights. 
     * @return index of chosen edge in array. 
     */
    private int getEdgeIndexStochastic(){
        BigDecimal summed = new BigDecimal(0);
        for (Edge edge: edges){
            summed = summed.add(edge.normalizedWeight);
        }
        
        int edgeIndex = 0; 
        // generate a random BigDecimal in a range: max * random
        BigDecimal choice = summed.multiply(new BigDecimal(RANDOM.nextDouble()));
        
        summed = edges.get(0).normalizedWeight;
        // do binary search
        while(summed.compareTo(choice) < 0){
            summed = summed.add(edges.get(edgeIndex).normalizedWeight);
            edgeIndex++;
        }
        
        if (Settings.getInteger("PTGLgraphComputation_I_debug_level") > 2) {
            System.out.println("[DEBUG LV 3] Summed Weight: " + summed.toString());
            System.out.println("[DEBUG LV 3] Chosen Weight: " + choice.toString());
            System.out.println("[DEBUG LV 3] Chosen Index: " + edgeIndex);
        }
        return edgeIndex; 
    }
    
    /**
     * Chooses an edge by equal probability. 
     * @return index in edge list
     */
    private int getEdgeIndexRandom(){
        return RANDOM.nextInt(edges.size());
    }
    
    /**
     * Return a formatted string representation of the current edge list.
     * @return multi-line string
     */
    private String prettyFormatEdgeList(int numberIndentSpaces, ClusteringResult clusteringResult) {
        String formattedString = " ".repeat(numberIndentSpaces + 4) + "rounded normalized weight | #res-res contacts | sub assembly 1 - sub assembly 2\n";
        Integer lineNum = 0;
        for (Edge edge: edges) {
            formattedString += " ".repeat(numberIndentSpaces + 2) + lineNum.toString() + ": " + 
                    edge.normalizedWeight.setScale(4, RoundingMode.HALF_UP) + " | " + 
                    edge.absoluteWeight + " | " +
                    clusteringResult.getVerticesFromRepresentative(edge.v1).stream().map(ele -> labelMap.get(ele)).collect(toList()) + " - " +
                    clusteringResult.getVerticesFromRepresentative(edge.v2).stream().map(ele -> labelMap.get(ele)).collect(toList()) + "\n";
            lineNum++;
        }
        return formattedString;
    }
    
    /**
     * Get the labelMap. 
     * @return labelMap
     */
    public Map<Integer,String> getlabelMap(){
        return labelMap;
    }
    
    /**
     * Tests agglomerative clustering implementation.
     */
    public static void main(String[] args) {
        System.out.println("CLUSTERING TEST");
        Settings.set("PTGLgraphComputation_I_debug_level", "3");
        Settings.set("PTGLgraphComputation_I_type_assembly_prediction", "2");
        
        // Example 1
        //   Additive length normalization merges
        //     2,3
        //     0,1
        //     0,2
        /*
        int[][] edgeList = {
            {0, 1, 10},
            {0, 2, 10},
            {1, 2, 5},
            {2, 3, 60}
        };
        
        Map<Integer, Integer> chainLengths = new HashMap<>();
        chainLengths.put(0, 5);
        chainLengths.put(1, 5);
        chainLengths.put(2, 15);
        chainLengths.put(3, 10);
        */
        
        // Example 2
        //   Additive length normalization merges
        //     5,8
        //     1,2
        //     5,7
        //     1,3
        //     5,6
        //     0,4
        //     0,5
        //     0,1
        int[][] edgeList = {
            {0, 1, 5},
            {0, 2, 4},
            {0, 3, 3},
            {0, 4, 10},
            {0, 5, 8},
            {0, 6, 7},
            {1, 2, 15},
            {2, 3, 14},
            {5, 7, 18},
            {5, 8, 21},
            {6, 7, 12}
        };
        
        Map<Integer, Double> chainProperties = new HashMap<>();
        chainProperties.put(0, 40.0);
        chainProperties.put(1, 8.0);
        chainProperties.put(2, 10.0);
        chainProperties.put(3, 20.0);
        chainProperties.put(4, 5.0);
        chainProperties.put(5, 12.0);
        chainProperties.put(6, 15.0);
        chainProperties.put(7, 11.0);
        chainProperties.put(8, 7.0);
        
        List<Chain> dummyChainList = null; // dummy list to satisfy the constructor; this list is only needed for gyradius normalization, so these normalizations are not possible to test here unless we create some suitable chain instances for this example
        Map<Integer, String> labelMap = new HashMap<>();
        for (Integer vertexId : chainProperties.keySet()) {
            labelMap.put(vertexId, String.valueOf(vertexId));
        }

        AgglomerativeClustering clustering = new AgglomerativeClustering(edgeList, chainProperties, EdgeWeightType.ADDITIVE_LENGTH_NORMALIZATION, labelMap, dummyChainList);
        
        ClusteringResult cr = clustering.chainInformationClustering();        
        
        System.out.println(cr.toNewickString());
        
        clustering = new AgglomerativeClustering(edgeList, chainProperties, EdgeWeightType.MULTIPLICATIVE_LENGTH_NORMALIZATION, labelMap, dummyChainList);
        clustering.chainInformationClustering();
        
        clustering = new AgglomerativeClustering(edgeList, chainProperties, EdgeWeightType.ABSOLUTE_WEIGHT, labelMap, dummyChainList);
        clustering.chainInformationClustering();
        
        System.out.println("FINISHED CLUSTERING TEST AND EXITING");
        System.exit(0);
    }
}
