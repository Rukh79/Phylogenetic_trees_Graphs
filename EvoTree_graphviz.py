import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

# Step 1: Define DNA sequences
sequences = {
    "Human": "ACTGAC",
    "Chimp": "ACTGAT",
    "Mouse": "GCTGAC",
    "Dog": "ACTAAC"
}

# Step 2: Function to calculate Hamming distance
def hamming(seq1, seq2):
    return sum(a != b for a, b in zip(seq1, seq2))

# Step 3: Build distance matrix
species = list(sequences.keys())
n = len(species)
distance_matrix = np.zeros((n, n))

for i in range(n):
    for j in range(i+1, n):
        dist = hamming(sequences[species[i]], sequences[species[j]])
        distance_matrix[i, j] = dist
        distance_matrix[j, i] = dist

# Step 4: UPGMA algorithm implementation
def upgma(distance_matrix, species):
    n = len(species)
    clusters = [[i] for i in range(n)]
    tree = nx.DiGraph()
    
    # Add initial nodes
    for i in range(n):
        tree.add_node(i, label=species[i], height=0)
    
    while len(clusters) > 1:
        # Find minimum distance
        min_dist = float('inf')
        min_i, min_j = -1, -1
        
        for i in range(len(clusters)):
            for j in range(i+1, len(clusters)):
                # Calculate average distance between clusters
                dist = 0
                count = 0
                for x in clusters[i]:
                    for y in clusters[j]:
                        dist += distance_matrix[x, y]
                        count += 1
                avg_dist = dist / count
                
                if avg_dist < min_dist:
                    min_dist = avg_dist
                    min_i, min_j = i, j
        
        # Merge clusters
        new_cluster = clusters[min_i] + clusters[min_j]
        new_node = len(tree.nodes)
        
        # Add new node
        tree.add_node(new_node, height=min_dist/2)
        
        # Add edges to children
        tree.add_edge(new_node, clusters[min_i][0], length=min_dist/2 - tree.nodes[clusters[min_i][0]]['height'])
        tree.add_edge(new_node, clusters[min_j][0], length=min_dist/2 - tree.nodes[clusters[min_j][0]]['height'])
        
        # Update clusters
        clusters = [c for idx, c in enumerate(clusters) if idx not in [min_i, min_j]] + [new_cluster]
    
    return tree

# Step 5: Build the tree
tree = upgma(distance_matrix, species)

# Step 6: Draw the tree
plt.figure(figsize=(10, 8))
# Use spring layout instead of graphviz
pos = nx.spring_layout(tree, k=2, iterations=50)
nx.draw(tree, pos, with_labels=True, 
        labels={node: tree.nodes[node].get('label', '') for node in tree.nodes},
        node_size=3000, node_color="skyblue", font_size=10, font_weight='bold',
        arrows=True)
plt.title("UPGMA Phylogenetic Tree")
plt.show()
