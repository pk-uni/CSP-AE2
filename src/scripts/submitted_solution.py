import networkx as nx
from ortools.linear_solver import pywraplp

from minizinc import Instance, Model, Solver


def build_dzn_graph(n, edges):
    G = [[0] * n for _ in range(n)]
    for u, v in edges:
        G[u][v] = 1
        G[v][u] = 1
    return G


def standardise_graph(nodes, edges):
    nodes_map = {node: idx for idx, node in enumerate(nodes)}
    new_edges = [(nodes_map[u], nodes_map[v]) for u, v in edges]
    return new_edges, nodes_map


def parse_networkx(instance_graph, start_node):
    edges, nodes_map = standardise_graph(instance_graph.nodes, instance_graph.edges)

    n = len(instance_graph.nodes)
    r = nodes_map[start_node]

    graph = build_dzn_graph(n, edges)

    return n, r, graph


#
# This function should run your ILP implementation
# for infectious vaccine
# - instance_graph will be a networkx graph object
# - start_node is the node where the fire starts
# - timeout is the maximum time in ms you should let the search run for
# The function should return a dictionary that has 'num_saved' mapped to
# the number saved in the solution
# or None if the model does not halt in the time allowed
# (The dictionary structure is so you can return other things if it's
# useful for your pipeline)
def run_ilp(instance_graph, start_node=1, timeout=1000):
    n, r, graph = parse_networkx(instance_graph, start_node)

    (solver, burned, defended) = solve_firefighter(n, r, graph)

    print_solution(solver, burned, defended, n, n * 3)

    return {"num_saved": -1}


#
# This function should run your ILP implementation
# for infectious vaccine
# - instance_graph will be a networkx graph object
# - start_node is the node where the fire starts
# - timeout is the maximum time in ms you should let the search run for
# The function should return a dictionary that has 'num_saved' mapped to
# the number saved in the solution
# or None if the model does not halt in the time allowed
# (The dictionary structure is so you can return other things if it's
# useful for your pipeline)
# There are many ways to deal with the instance_graph, and
# this is intentionally left up to you.
# For example, you could create a .dzn file in whatever encoding you want
# and add it using the https://python.minizinc.dev/en/latest/api.html#minizinc.model.Model.add_file capability
def run_cp(instance_graph, start_node=1, timeout=1000):
    gecode = Solver.lookup("gecode")

    base = Model()
    base.add_file("src/cp/infectious-defence.mzn")

    instance = Instance(gecode, base)

    n, r, graph = parse_networkx(instance_graph, start_node)

    T = 3 * n  # T should be equal to the diameter of the graph
    instance["n"] = n
    instance["T"] = T
    instance["r"] = r
    instance["graph"] = graph

    result = instance.solve()
    num_saved = getattr(result.solution, "objective", -1)

    return {"num_saved": num_saved}


# ----------


def solve_firefighter(n, r, graph):

    solver = pywraplp.Solver.CreateSolver("SCIP")

    # TODO: use graph diameter (3 phases per round, diameter(graph) rounds)
    T = n * 3
    d = 1
    burned = {v: {} for v in range(1, n + 1)}
    defended = {v: {} for v in range(1, n + 1)}

    for v in range(1, n + 1):
        for t in range(T + 1):
            burned[v][t] = solver.IntVar(0, 1, f"burned_{v}_{t}")
            defended[v][t] = solver.IntVar(0, 1, f"defended_{v}_{t}")

    neighbours = lambda x, y: graph[x - 1][y - 1] == 1

    # Initial conditions
    # root starts burning
    solver.Add(burned[r][0] == 1)

    # all other vertices start unburned
    for v in range(1, n + 1):
        if v != r:
            solver.Add(burned[v][0] == 0)

    # no vertices start defended
    for v in range(1, n + 1):
        solver.Add(defended[v][0] == 0)

    # State exclusivity
    for t in range(T + 1):
        for v in range(1, n + 1):
            solver.Add(burned[v][t] + defended[v][t] <= 1)

    # Permanence constraints
    for t in range(1, T + 1):
        for v in range(1, n + 1):
            solver.Add(burned[v][t] >= burned[v][t - 1])
            solver.Add(defended[v][t] >= defended[v][t - 1])

    # phase 1: defence placement
    for t in range(1, T + 1, 3):
        solver.Add(
            sum(defended[v][t] - defended[v][t - 1] for v in range(1, n + 1)) <= d
        )
        # only defence placement in phase 1
        for v in range(1, n + 1):
            solver.Add(burned[v][t] == burned[v][t - 1])

    # phase 2: fire spread
    for t in range(2, T + 1, 3):
        for x in range(1, n + 1):
            for y in range(1, n + 1):
                if neighbours(x, y):
                    solver.Add(burned[x][t] + defended[x][t] >= burned[y][t - 1])

        # only fire spread in phase 2
        for x in range(1, n + 1):
            solver.Add(defended[x][t] == defended[x][t - 1])

    # phase 3: defence spread
    for t in range(3, T + 1, 3):
        for x in range(1, n + 1):
            for y in range(1, n + 1):
                if neighbours(x, y):
                    solver.Add(
                        defended[x][t]
                        >= defended[y][t - 1] + (1 - burned[x][t - 1]) - 1
                    )

        # only defence spread in phase 3
        for x in range(1, n + 1):
            solver.Add(burned[x][t] == burned[x][t - 1])

    # minimise total burned vertices at time T
    objective = solver.Objective()

    for v in range(1, n + 1):
        objective.SetCoefficient(burned[v][T], 1)

    objective.SetMinimization()

    return solver, burned, defended


# nuke if not used
def print_solution(solver, burned, defended, n, T):
    """Helper function to print the solution in a readable format"""
    if solver.Solve() == pywraplp.Solver.OPTIMAL:
        print(f"Optimal solution found!")
        for t in range(T + 1):
            print(f"\nt={t} (Round {t//3 + 1}, Phase {t%3 + 1})")
            print(
                "Burned vertices:",
                [v for v in range(1, n + 1) if burned[v][t].solution_value() > 0.5],
            )
            print(
                "Defended vertices:",
                [v for v in range(1, n + 1) if defended[v][t].solution_value() > 0.5],
            )
    else:
        print("No solution found.")
