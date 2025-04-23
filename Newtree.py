import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import plotly.graph_objects as go
import networkx as nx
import numpy as np

# Initialize Dash app
app = dash.Dash(__name__)

# Define DNA sequences
sequences = {
    "Human": "ACTGAC",
    "Chimp": "ACTGAT",
    "Mouse": "GCTGAC",
    "Dog": "ACTAAC"
}

# Function to calculate Hamming distance
def hamming(seq1, seq2):
    return sum(a != b for a, b in zip(seq1, seq2))

# Build distance matrix
species = list(sequences.keys())
n = len(species)
distance_matrix = np.zeros((n, n))

for i in range(n):
    for j in range(i+1, n):
        dist = hamming(sequences[species[i]], sequences[species[j]])
        distance_matrix[i, j] = dist
        distance_matrix[j, i] = dist

# Create a simple tree with distances
def create_tree():
    G = nx.DiGraph()
    
    # Add nodes
    for i, sp in enumerate(species):
        G.add_node(i, label=sp)
    
    # Add edges with distances
    G.add_edge(4, 0, length=1.0)  # Root to Human
    G.add_edge(4, 1, length=1.0)  # Root to Chimp
    G.add_edge(5, 2, length=2.0)  # Root to Mouse
    G.add_edge(5, 3, length=2.0)  # Root to Dog
    G.add_edge(6, 4, length=0.5)  # Root to Human/Chimp cluster
    G.add_edge(6, 5, length=1.5)  # Root to Mouse/Dog cluster
    
    return G

def create_tree_visualization():
    # Create the tree
    tree = create_tree()
    
    # Position nodes using spring layout
    pos = nx.spring_layout(tree, k=2, iterations=50)
    
    # Create edge traces with distances
    edge_x = []
    edge_y = []
    edge_text = []
    for edge in tree.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
        dist = tree[edge[0]][edge[1]]['length']
        edge_text.append(f"Distance: {dist:.2f}")
    
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=2, color='black'),
        hoverinfo='text',
        text=edge_text,
        mode='lines'
    )
    
    # Create node traces
    node_x = []
    node_y = []
    node_text = []
    for node in tree.nodes():
        if 'label' in tree.nodes[node]:
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)
            node_text.append(tree.nodes[node]['label'])
    
    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        text=node_text,
        textposition="bottom center",
        hoverinfo='text',
        marker=dict(
            size=20,
            color='blue',
            line=dict(width=2)
        )
    )
    
    # Create annotations for distances
    annotations = []
    for edge in tree.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        dist = tree[edge[0]][edge[1]]['length']
        annotations.append(
            dict(
                x=(x0 + x1) / 2,
                y=(y0 + y1) / 2,
                text=f"{dist:.2f}",
                showarrow=False,
                font=dict(size=12, color='red'),
                bgcolor='white',
                bordercolor='black',
                borderwidth=1,
                borderpad=4
            )
        )
    
    fig = go.Figure(
        data=[edge_trace, node_trace],
        layout=go.Layout(
            showlegend=False,
            hovermode='closest',
            margin=dict(b=20, l=5, r=5, t=40),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            annotations=annotations,
            plot_bgcolor='white'
        )
    )
    
    return fig

app.layout = html.Div([
    html.H1("Phylogenetic Tree with Distances"),
    dcc.Graph(id='phylogenetic-tree', figure=create_tree_visualization())
])

if __name__ == '__main__':
    app.run(debug=True, host='127.0.0.1', port=8050)