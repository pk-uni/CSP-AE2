from ortools.linear_solver import pywraplp


def generate_solver(n, r, D, graph):

    solver = pywraplp.Solver.CreateSolver("SCIP")

    DEFENCE_BUDGET = 1
    PHASES = 3
    T = D * PHASES

    burned = {}
    defended = {}

    for v in range(n):
        for t in range(T + 1):
            burned[(v, t)] = solver.IntVar(0, 1, f"burned_{v}_{t}")
            defended[(v, t)] = solver.IntVar(0, 1, f"defended_{v}_{t}")

    neighbours = lambda x, y: graph[x - 1][y - 1] == 1

    # root starts burning
    solver.Add(burned[(r, 0)] == 1)

    # all other vertices start unburned
    for v in range(n):
        if v != r:
            solver.Add(burned[(v, 0)] == 0)

    # no vertices start defended
    for v in range(n):
        solver.Add(defended[(v, 0)] == 0)

    # state exclusivity
    for t in range(T + 1):
        for v in range(n):
            solver.Add(burned[(v, t)] + defended[(v, t)] <= 1)
            # vertices can either be: burning, defended or unaffected
            solver.Add(burned[(v, t)] <= 1)
            solver.Add(burned[(v, t)] >= 0)
            solver.Add(defended[(v, t)] <= 1)
            solver.Add(defended[(v, t)] >= 0)

    # permanence
    for t in range(1, T + 1):
        for v in range(n):
            solver.Add(burned[(v, t)] >= burned[(v, t - 1)])
            solver.Add(defended[(v, t)] >= defended[(v, t - 1)])

    # phase 1: defense placement
    for t in range(1, T + 1, PHASES):
        solver.Add(
            sum(defended[(v, t)] - defended[(v, t - 1)] for v in range(n))
            <= DEFENCE_BUDGET
        )
        # only defence placement in phase 1
        for v in range(n):
            solver.Add(burned[(v, t)] == burned[(v, t - 1)])

    # phase 2: fire spread
    for t in range(2, T + 1, PHASES):
        for x in range(n):
            for y in range(n):
                if neighbours(x, y):
                    solver.Add(burned[(x, t)] + defended[(x, t)] >= burned[(y, t - 1)])

        # only fire spread in phase 2
        for x in range(n):
            solver.Add(defended[(x, t)] == defended[(x, t - 1)])

    # phase 3: defence spread
    for t in range(3, T + 1, PHASES):
        for x in range(n):
            for y in range(n):
                if neighbours(x, y):
                    solver.Add(
                        defended[(x, t)]
                        >= defended[(y, t - 1)] + (1 - burned[(x, t - 1)]) - 1
                    )

        # only defence spread in phase 3
        for x in range(n):
            solver.Add(burned[(x, t)] == burned[(x, t - 1)])

    # minimise total burned vertices at time T
    objective = solver.Objective()

    for v in range(n):
        objective.SetCoefficient(burned[(v, T)], 1)

    objective.SetMinimization()

    return solver, burned, defended


def solve(solver, burned, defended, n, T):
    status = solver.Solve()

    for ti in range(T + 1):
        defended_a = [f"{defended[(v, ti)].solution_value():.3f}" for v in range(n)]
        burned_a = [f"{burned[(v, ti)].solution_value():.3f}" for v in range(n)]
        print(f"-- {defended_a}")
        print(f"-- {burned_a}")
        print()

    if status == pywraplp.Solver.OPTIMAL:
        final_burned = [v for v in range(n) if burned[(v, T)].solution_value() > 0.5]
        return n - len(final_burned)

    elif status == pywraplp.Solver.INFEASIBLE:
        print("Problem is infeasible")
        return -1

    elif status == pywraplp.Solver.UNBOUNDED:
        print("Problem is unbounded")
        return -1
