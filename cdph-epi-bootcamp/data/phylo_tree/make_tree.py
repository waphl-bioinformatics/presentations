"""Simulate a 100-tip phylogenetic tree (coalescent) and render an interactive plotly HTML."""
import random, math
import plotly.graph_objects as go

random.seed(42)
N = 100

class Node:
    _i = 0
    def __init__(self, name=None, height=0.0):
        self.name = name
        self.height = height        # time before present
        self.children = []
        self.parent = None
        self.x = self.y = None
        self.id = Node._i; Node._i += 1

# --- Kingman coalescent: 100 tips at height 0, coalesce upward -------------
tips = [Node(name=f"taxon_{i+1:03d}") for i in range(N)]
active = tips[:]
t = 0.0
while len(active) > 1:
    k = len(active)
    t += random.expovariate(k * (k - 1) / 2.0)   # waiting time to next coalescence
    a, b = random.sample(active, 2)
    p = Node(height=t)
    p.children = [a, b]
    a.parent = b.parent = p
    active.remove(a); active.remove(b); active.append(p)
root = active[0]

# --- ladderize + assign coordinates ---------------------------------------
def n_leaves(n):
    if not n.children: return 1
    return sum(n_leaves(c) for c in n.children)

def ladderize(n):
    if n.children:
        n.children.sort(key=n_leaves)
        for c in n.children: ladderize(c)
ladderize(root)

leaf_counter = [0]
def layout(n):
    n.x = root.height - n.height          # distance from root
    if not n.children:
        n.y = leaf_counter[0]; leaf_counter[0] += 1
    else:
        for c in n.children: layout(c)
        n.y = sum(c.y for c in n.children) / len(n.children)
layout(root)

# --- color tips by clade defined by cutting the tree at 40% of root height --
cut = 0.40 * root.height
palette = ["#4C78A8","#F58518","#54A24B","#E45756","#72B7B2","#B279A2",
           "#EECA3B","#FF9DA6","#9D755D","#BAB0AC","#3E6D9C","#D67C3E"]
clade_of = {}
def assign(n, cid):
    if cid is None and n.height <= cut and n.parent is not None:
        cid = assign.next; assign.next += 1
    for c in n.children: assign(c, cid)
    if not n.children: clade_of[n.id] = cid if cid is not None else -1
assign.next = 0
assign(root, None)
color = lambda n: palette[clade_of[n.id] % len(palette)] if clade_of[n.id] >= 0 else "#888"

# --- build rectangular-phylogram line segments -----------------------------
fig = go.Figure()
segs_x, segs_y = [], []          # grey structural lines (drawn as one trace)
def walk(n):
    for c in n.children:
        segs_x.extend([n.x, n.x, None])       # vertical connector
        segs_y.extend([n.y, c.y, None])
        segs_x.extend([n.x, c.x, None])       # horizontal branch
        segs_y.extend([c.y, c.y, None])
        walk(c)
walk(root)

fig.add_trace(go.Scatter(x=segs_x, y=segs_y, mode="lines",
                         line=dict(color="#5a5a5a", width=1),
                         hoverinfo="skip", showlegend=False))

# dotted alignment lines out to the tip labels
tip_nodes = [n for n in tips]
maxx = max(n.x for n in tip_nodes)
dash_x, dash_y = [], []
for n in tip_nodes:
    dash_x.extend([n.x, maxx * 1.02, None]); dash_y.extend([n.y, n.y, None])
fig.add_trace(go.Scatter(x=dash_x, y=dash_y, mode="lines",
                         line=dict(color="#cccccc", width=0.6, dash="dot"),
                         hoverinfo="skip", showlegend=False))

# tips
fig.add_trace(go.Scatter(
    x=[n.x for n in tip_nodes], y=[n.y for n in tip_nodes],
    mode="markers+text",
    marker=dict(size=7, color=[color(n) for n in tip_nodes],
                line=dict(width=0.5, color="white")),
    text=[" " + n.name for n in tip_nodes], textposition="middle right",
    textfont=dict(size=8, color="#333"),
    customdata=[[n.name, round(n.x, 4)] for n in tip_nodes],
    hovertemplate="<b>%{customdata[0]}</b><br>root distance: %{customdata[1]}<extra></extra>",
    showlegend=False))

# internal nodes (hoverable)
internal = []
def collect(n):
    if n.children:
        internal.append(n)
        for c in n.children: collect(c)
collect(root)
fig.add_trace(go.Scatter(
    x=[n.x for n in internal], y=[n.y for n in internal],
    mode="markers", marker=dict(size=4, color="#333", opacity=0.55),
    customdata=[[round(n.height, 4), n_leaves(n)] for n in internal],
    hovertemplate="internal node<br>age: %{customdata[0]}<br>descendant tips: %{customdata[1]}<extra></extra>",
    showlegend=False))

fig.update_layout(
    title=dict(text="Simulated phylogeny — 100 tips (Kingman coalescent)", x=0.02, font=dict(size=18)),
    xaxis=dict(title="branch length (substitutions / time from root)",
               showgrid=True, gridcolor="#eee", zeroline=False, range=[-0.02 * maxx, maxx * 1.22]),
    yaxis=dict(showticklabels=False, showgrid=False, zeroline=False,
               range=[-1.5, N + 0.5]),
    plot_bgcolor="white", paper_bgcolor="white",
    width=1000, height=1500, margin=dict(l=40, r=20, t=60, b=50),
    hovermode="closest")

fig.write_html("/mnt/user-data/outputs/phylo_tree_100tips.html", include_plotlyjs="cdn")

# --- also write the Newick string ------------------------------------------
def newick(n):
    if not n.children:
        return f"{n.name}:{n.parent.height - n.height:.6f}"
    inner = ",".join(newick(c) for c in n.children)
    bl = "" if n.parent is None else f":{n.parent.height - n.height:.6f}"
    return f"({inner}){bl}"
with open("/mnt/user-data/outputs/tree_100tips.nwk", "w") as f:
    f.write(newick(root) + ";\n")

print("root height:", round(root.height, 4), "| tips:", len(tips), "| internal nodes:", len(internal))
