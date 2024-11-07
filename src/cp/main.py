from minizinc import Instance, Model, Solver


def solve(instance, n):
    result = instance.solve()
    burned = getattr(result.solution, "objective", -1)
    if burned != -1:
        return n - burned

    return -1


def generate_solver(n, r, graph):
    gecode = Solver.lookup("gecode")

    base = Model()
    base.add_file("src/cp/infectious-defence.mzn")

    instance = Instance(gecode, base)

    T = 3 * n  # T should be equal to the diameter of the graph
    instance["n"] = n
    instance["T"] = T
    instance["r"] = r
    instance["graph"] = graph

    return instance
