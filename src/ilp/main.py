from ortools.linear_solver import pywraplp


def solve_firefighter(n, r, graph):

    solver = pywraplp.Solver.CreateSolver("SCIP")

    # TODO: use graph diameter (3 phases per round, diameter(graph) rounds)
    T = n * 3
    d = 1
    burned = {}
    defended = {}

    for v in range(1, n + 1):
        for t in range(T + 1):
            burned[v, t] = solver.IntVar(0, 1, f"burned_{v}_{t}")
            defended[v, t] = solver.IntVar(0, 1, f"defended_{v}_{t}")

    neighbours = lambda x, y: graph[x - 1][y - 1] == 1

    # Initial conditions
    # root starts burning
    solver.Add(burned[r, 0] == 1)

    # all other vertices start unburned
    for v in range(1, n + 1):
        if v != r:
            solver.Add(burned[v, 0] == 0)

    # no vertices start defended
    for v in range(1, n + 1):
        solver.Add(defended[v, 0] == 0)

    # State exclusivity
    for t in range(T + 1):
        for v in range(1, n + 1):
            solver.Add(burned[v, t] + defended[v, t] <= 1)

    # Permanence constraints
    for t in range(1, T + 1):
        for v in range(1, n + 1):
            solver.Add(burned[v, t] >= burned[v, t - 1])
            solver.Add(defended[v, t] >= defended[v, t - 1])

    # phase 1: defense placement
    for t in range(1, T + 1, 3):
        solver.Add(
            sum(defended[v, t] - defended[v, t - 1] for v in range(1, n + 1)) <= d
        )
        # only defence placement in phase 1
        for v in range(1, n + 1):
            solver.Add(burned[v, t] == burned[v, t - 1])

    # phase 2: fire spread
    for t in range(2, T + 1, 3):
        for x in range(1, n + 1):
            for y in range(1, n + 1):
                if neighbours(x, y):
                    solver.Add(burned[x, t] + defended[x, t] >= burned[y, t - 1])

        # only fire spread in phase 2
        for x in range(1, n + 1):
            solver.Add(defended[x, t] == defended[x, t - 1])

    # phase 3:
    for t in range(3, T + 1, 3):
        for x in range(1, n + 1):
            for y in range(1, n + 1):
                if graph[x - 1][y - 1] == 1:
                    solver.Add(
                        defended[x, t]
                        >= defended[y, t - 1] + (1 - burned[x, t - 1]) - 1
                    )

        # Fire can't spread during defense spread
        for x in range(1, n + 1):
            solver.Add(burned[x, t] == burned[x, t - 1])

    # minimise total burned vertices at time T
    objective = solver.Objective()

    for v in range(1, n + 1):
        objective.SetCoefficient(burned[v, T], 1)

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
                [v for v in range(1, n + 1) if burned[v, t].solution_value() > 0.5],
            )
            print(
                "Defended vertices:",
                [v for v in range(1, n + 1) if defended[v, t].solution_value() > 0.5],
            )
    else:
        print("No solution found.")
