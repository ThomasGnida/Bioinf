import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import networkx as nx
import pandas as pd
import heapq

def dijkstra(G, source):
    """
    Dijkstra's algorithm using combined_score as edge weight.
    Returns distances and predecessors from source to all other nodes.
    :param G: NetworkX graph
    :param source: Source node
    """
    distances = {node: float('inf') for n in G.nodes()}
    distances[source] = 0
    predecessors = {node: None for n in G.nodes()}

    pq = [(0, source)]
    visited = set()

    while pq:
        current_dist, current_node = heapq.heappop(pq)

        if current_node in visited:
            continue
        else:
            visited.add(current_node)

        for neighbor, edge_data in G[current_node].items():
            # Use 1 - combined_score as weight (higher score = shorter distance)
            weight = 1 - edge_data.get('combined_score', 0)
            new_dist = current_dist + weight

            if new_dist < distances[neighbor]:
                distances[neighbor] = new_dist
                predecessors[neighbor] = current_node
                heapq.heappush(pq, (new_dist, neighbor))

    return distances, predecessors


def reconstruct_path(predecessors, source, target):
    """Reconstruct path from source to target using predecessor map."""
    path = []
    current = target
    while current is not None:
        path.append(current)
        current = predecessors[current]
    path.reverse()
    return path if path[0] == source else []


'''
A + B : Downloading and loading PPI network.
The dataset was downloaded from string-db.org, Homo Sapiens Hemoglobin was chosen and extended until at least
50 nodes were present.
'''

df = pd.read_csv('Assignment_5_Graph.tsv', sep='\t', comment=None)
# Strip the # from the first column name
df.columns = df.columns.str.lstrip('#')
G = nx.from_pandas_edgelist(df, source='node1', target='node2', edge_attr='combined_score')

'''
C: Basic graph statistics
Number of nodes and edges:

Node and edge connectivity:

Completeness:

Number and size of components:
'''
print(f"Nodes: {G.number_of_nodes()}")
print(f"Edges: {G.number_of_edges()}")
print(f"Node connectivity: {nx.node_connectivity(G)}")
print(f"Edge connectivity: {nx.edge_connectivity(G)}")
n = G.number_of_nodes()
max_edges = n * (n - 1) / 2
is_complete = G.number_of_edges() == max_edges
print(f"Is complete: {is_complete}")
components = list(nx.connected_components(G))
print(f"Number of connected components: {len(components)}")
print(f"Component sizes: {sorted([len(c) for c in components], reverse=True)}")


'''
D: Centrality measures:

Degree centrality:

Betweenness Centrality:
'''
degree_centrality = nx.degree_centrality(G)
print("Degree Centrality (Top 5):")
for node, score in sorted(degree_centrality.items(), key=lambda x: x[1], reverse=True)[:5]:
    print(f" {node}: {score:.3f}")

betweenness_centrality = nx.betweenness_centrality(G)
print("\nBetweenness Centrality (Top 5):")
for node, score in sorted(betweenness_centrality.items(), key=lambda x: x[1], reverse=True)[:5]:
    print(f" {node}: {score:.3f}")

'''
E: Louvain communities

Number of communities detected:

Modularity score:
'''
largest_cc = components[0]
G_lcc = G.subgraph(largest_cc).copy()

communities = nx.algorithms.community.louvain_communities(G_lcc, seed=42)

print(f"Number of communities: {len(communities)}")
for i, community in enumerate(communities):
    print(f"  Community {i+1} (size {len(community)}): {sorted(community)}")

mod_score = nx.algorithms.community.modularity(G_lcc, communities)
print(f"\nModularity score: {mod_score:.4f}")

'''
F: Shortest path between two proteins
Selected Proteins:
Cytochrome c (CYCS)
Albumin (ALB)
'''
source_protein = 'CYCS'
target_protein = 'ALB'

distances, predecessors = dijkstra(G, source_protein)
path = reconstruct_path(predecessors, source_protein, target_protein)

print(f"Shortest path from {source_protein} to {target_protein}:")
print(f"  Path: {' -> '.join(path)}")
print(f"  Weighted distance: {distances[target_protein]:.4f}")
print(f"  Number of hops: {len(path) - 1}")

'''
G: Average shortest path

Average shortest path length (aspl):

Diameter:
'''
avg_shortest_path = nx.average_shortest_path_length(G_lcc)
diameter = nx.diameter(G_lcc)

print(f"\nAverage shortest path length: {avg_shortest_path:.4f}")
print(f"Diameter: {diameter}")

'''
H: Visualisation

Node colour = community
Node size = degree centrality
Edge highlighting
'''
community_colors = {}
palette = list(mcolors.TABLEAU_COLORS.values())
for i, community in enumerate(communities):
    for node in community:
        community_colors[node] = palette[i % len(palette)]

node_sizes = [degree_centrality[node] * 3000 for node in G_lcc.nodes()]

node_colors = [community_colors[node] for node in G_lcc.nodes()]

path_edges = list(zip(path[:-1], path[1:]))
edge_colors = []
edge_widths = []
for u, v in G_lcc.edges():
    if (u, v) in path_edges or (v, u) in path_edges:
        edge_colors.append('#FF3333')
        edge_widths.append(4.0)
    else:
        edge_colors.append('#cccccc')
        edge_widths.append(0.8)

pos = nx.kamada_kawai_layout(G_lcc)

fig, ax = plt.subplots(figsize=(18, 14))
ax.set_facecolor('#f9f9f9')
fig.patch.set_facecolor('#f9f9f9')

nx.draw_networkx_edges(
    G_lcc, pos,
    edge_color=edge_colors,
    width=edge_widths,
    alpha=0.7,
    ax=ax
)

nx.draw_networkx_nodes(
    G_lcc, pos,
    node_color=node_colors,
    node_size=node_sizes,
    alpha=0.95,
    ax=ax
)

nx.draw_networkx_labels(
    G_lcc, pos,
    font_size=7,
    font_weight='bold',
    font_color='black',
    ax=ax
)

legend_patches = [
    mpatches.Patch(color=palette[i], label=f'Community {i+1} (n={len(c)})')
    for i, c in enumerate(communities)
]
legend_patches.append(mpatches.Patch(color='#FF3333', label=f'Shortest path: {source_protein} → {target_protein}'))

ax.legend(handles=legend_patches, loc='upper left', fontsize=9, framealpha=0.9)

ax.set_title(
    f'Myoglobin PPI Network (n={G_lcc.number_of_nodes()} proteins, {G_lcc.number_of_edges()} interactions)\n'
    f'Node size = degree centrality | Node color = Louvain community | Red edges = shortest path',
    fontsize=13, fontweight='bold', pad=20
)
ax.axis('off')
plt.tight_layout()
plt.savefig('ppi_network.png', dpi=150, bbox_inches='tight')
plt.show()
