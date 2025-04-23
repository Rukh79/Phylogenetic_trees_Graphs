import networkx as nx
import numpy as np
import plotly.graph_objects as go
from Bio import Entrez, SeqIO, Phylo
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random
import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import plotly.express as px
import tempfile
import io

# Set your email for NCBI
Entrez.email = "f20220849@pilani.bits-pilani.ac.in"

# Define mitochondrial genes to analyze
MITO_GENES = {
    "COX1": {"start": 5904, "end": 7445},  # Human mitochondrial coordinates
    "COX2": {"start": 7586, "end": 8269},
    "COX3": {"start": 9207, "end": 9990},
    "CYTB": {"start": 14747, "end": 15887},
    "ND2": {"start": 4470, "end": 5511}
}

# Define species and their NCBI accession numbers
species_data = {
    "Human": "NC_012920.1",
    "Chimpanzee": "NC_001643.1",
    "Mouse": "NC_005089.1",
    "Dog": "NC_002008.4"
}

class PhylogeneticAnalyzer:
    def __init__(self):
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'global'
        self.aligner.match_score = 1
        self.aligner.mismatch_score = -1
        self.aligner.open_gap_score = -1
        self.aligner.extend_gap_score = -0.5
        
    def extract_gene_sequence(self, accession, gene_name):
        """Extract gene sequence from mitochondrial genome"""
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            
            for feature in record.features:
                if feature.type == "CDS" and gene_name in feature.qualifiers.get("gene", [""])[0]:
                    return str(feature.extract(record.seq))
            return None
        except Exception as e:
            print(f"Error extracting {gene_name} from {accession}: {str(e)}")
            return None

    def calculate_similarity(self, seq1, seq2):
        """Calculate sequence similarity with alignment"""
        alignment = self.aligner.align(seq1, seq2)[0]
        return 1 - (alignment.score / max(len(seq1), len(seq2)))

    def bootstrap_analysis(self, sequences, n_bootstraps=100):
        """Perform bootstrap analysis"""
        bootstrap_trees = []
        seq_length = len(next(iter(sequences.values())))
        
        for _ in range(n_bootstraps):
            # Resample alignment columns
            resampled_seqs = {}
            for species, seq in sequences.items():
                resampled_indices = random.choices(range(seq_length), k=seq_length)
                resampled_seqs[species] = ''.join(seq[i] for i in resampled_indices)
            
            # Build tree from resampled data
            tree = self.build_tree(resampled_seqs)
            bootstrap_trees.append(tree)
        
        return bootstrap_trees

    def build_tree(self, sequences):
        """Build phylogenetic tree using neighbor-joining"""
        # Create a distance matrix
        species = list(sequences.keys())
        n = len(species)
        distance_matrix = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i+1, n):
                dist = self.calculate_similarity(sequences[species[i]], sequences[species[j]])
                distance_matrix[i, j] = dist
                distance_matrix[j, i] = dist
        
        # Create a Phylo tree using the distance matrix
        tree = self._create_phylo_tree(distance_matrix, species)
        return tree

    def _create_phylo_tree(self, distance_matrix, species):
        """Create a Phylo tree from distance matrix"""
        # Create a Newick string
        newick = "("
        for i in range(len(species)):
            newick += f"{species[i]}:{distance_matrix[i,0]}"
            if i < len(species)-1:
                newick += ","
        newick += ");"
        
        # Parse the Newick string into a Phylo tree
        tree = Phylo.read(io.StringIO(newick), "newick")
        return tree

    def analyze_multiple_genes(self, species_data):
        """Analyze multiple mitochondrial genes"""
        all_sequences = {}
        bootstrap_support = {}
        
        for gene_name in MITO_GENES:
            print(f"Analyzing {gene_name}...")
            sequences = {}
            for species, accession in species_data.items():
                seq = self.extract_gene_sequence(accession, gene_name)
                if seq:
                    sequences[species] = seq
            
            if sequences:
                # Perform bootstrap analysis
                bootstrap_trees = self.bootstrap_analysis(sequences)
                bootstrap_support[gene_name] = bootstrap_trees
                all_sequences[gene_name] = sequences
        
        return all_sequences, bootstrap_support

def create_interactive_tree(tree):
    """Create interactive tree visualization using Plotly"""
    # Create traces for edges and nodes
    edge_trace = go.Scatter(
        x=[],
        y=[],
        line=dict(width=2, color='black'),
        hoverinfo='none',
        mode='lines'
    )
    
    node_trace = go.Scatter(
        x=[],
        y=[],
        text=[],
        mode='markers+text',
        hoverinfo='text',
        marker=dict(
            size=10,
            color='blue',
            line=dict(width=2)
        )
    )
    
    # Add edges
    for clade in tree.get_nonterminals():
        for child in clade.clades:
            edge_trace['x'] += [clade.branch_length, child.branch_length, None]
            edge_trace['y'] += [0, 0, None]
    
    # Add nodes
    for clade in tree.get_terminals():
        node_trace['x'] += [clade.branch_length]
        node_trace['y'] += [0]
        node_trace['text'] += [clade.name]
    
    # Create figure
    fig = go.Figure(
        data=[edge_trace, node_trace],
        layout=go.Layout(
            showlegend=False,
            hovermode='closest',
            margin=dict(b=20, l=5, r=5, t=40),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)
        )
    )
    
    return fig

# Initialize Dash app
app = dash.Dash(__name__)

app.layout = html.Div([
    html.H1("Interactive Phylogenetic Analysis"),
    dcc.Graph(id='phylogenetic-tree'),
    html.Div([
        html.Label("Select Gene:"),
        dcc.Dropdown(
            id='gene-selector',
            options=[{'label': gene, 'value': gene} for gene in MITO_GENES.keys()],
            value='COX1'
        ),
    ]),
    html.Div(id='bootstrap-support')
])

@app.callback(
    Output('phylogenetic-tree', 'figure'),
    [Input('gene-selector', 'value')]
)
def update_tree(selected_gene):
    analyzer = PhylogeneticAnalyzer()
    sequences, bootstrap_support = analyzer.analyze_multiple_genes(species_data)
    
    if selected_gene in sequences:
        tree = analyzer.build_tree(sequences[selected_gene])
        fig = create_interactive_tree(tree)
        return fig
    return {}

if __name__ == '__main__':
    # Run the Dash app with proper configuration
    app.run(debug=True, host='127.0.0.1', port=8050)