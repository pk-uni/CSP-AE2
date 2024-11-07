from minizinc import Instance, Model, Solver


def generate_solver(n, r, D, graph):
    DEFENCE_BUDGET = 1
    PHASES = 3
    T = D * PHASES

    gecode = Solver.lookup("gecode")

    model = Model()
    model.add_file("src/cp/infectious-defence.mzn")

    instance = Instance(gecode, model)

    instance["n"] = n
    instance["r"] = r
    instance["T"] = T
    instance["d"] = DEFENCE_BUDGET
    instance["graph"] = graph

    return instance


def solve(instance, n):
    result = instance.solve()
    burned = getattr(result.solution, "objective", -1)
    if burned != -1:
        return n - burned

    return -1
