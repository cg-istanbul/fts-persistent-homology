# scripts/report/viz.sage

from sage.all import Graphics, point, Graph
load("scripts/report/outputs.sage")

def save_scatter_plot(points, name, outfile=None):
    if outfile is None:
        outfile = scatter_png_path(name)   
    G = Graphics()
    G += point(points, size=18)
    G.save(outfile)
    print("[plot] saved to:", outfile)
    return outfile


def save_mapper_graph_png(mapper_out, name, outfile=None, vertex_labels=True, layout="spring"):
    if outfile is None:
        outfile = mapper_png_path(name)
    G = mapper_out.get("G", None) if isinstance(mapper_out, dict) else None
    if G is None:
        nodes = mapper_out.get("nodes", []) if isinstance(mapper_out, dict) else []
        edges = mapper_out.get("edges", []) if isinstance(mapper_out, dict) else []
        G = Graph()
        G.add_vertices([int(nd.get("id", 0)) for nd in nodes])
        G.add_edges([(int(a), int(b)) for (a, b) in edges])
    P = G.plot(vertex_labels=bool(vertex_labels), layout=layout)
    P.save(outfile)
    return outfile
