def parse_networkx(instance_graph, start_node):
    nodes_map = {node: idx for idx, node in enumerate(instance_graph.nodes)}
    edges = [(nodes_map[u], nodes_map[v]) for u, v in instance_graph.edges]

    n = len(instance_graph.nodes)
    r = nodes_map[start_node]

    graph = [[0] * n for _ in range(n)]
    for u, v in edges:
        graph[u][v] = 1
        graph[v][u] = 1

    return n, r, graph
