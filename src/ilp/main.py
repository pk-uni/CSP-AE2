from ortools.linear_solver import pywraplp
import networkx as nx
from typing import Dict, Tuple, List


class InfectiousFirefighter:
    def __init__(self, graph: nx.Graph, root: int, defenders: int):
        """
        Initialize the Infectious Firefighter solver.

        Args:
            graph (nx.Graph): The input graph
            root (int): The starting vertex for the fire
            defenders (int): Number of defenders that can be placed per round
        """
        # Store input parameters
        self.graph = graph
        self.root = root
        self.defenders = defenders

        # Derived parameters
        self.n_vertices = graph.number_of_nodes()
        self.phases = 3
        self.T = self.n_vertices * self.phases  # Maximum possible timesteps

        # Initialize solver
        self.solver = pywraplp.Solver.CreateSolver("SCIP")

        # Initialize variables dictionary
        self.burned: Dict[Tuple[int, int], pywraplp.Variable] = {}
        self.defended: Dict[Tuple[int, int], pywraplp.Variable] = {}

        # Create variables
        self._create_variables()

    def _create_variables(self):
        """Create binary variables for burned and defended states."""
        # Create variables for each vertex and timestep
        for v in self.graph.nodes():
            for t in range(self.T + 1):  # Include time 0
                # Binary variable for burned state
                self.burned[v, t] = self.solver.BoolVar(f"b_{v}_{t}")
                # Binary variable for defended state
                self.defended[v, t] = self.solver.BoolVar(f"d_{v}_{t}")

    def get_phase(self, t: int) -> int:
        """Get the phase (1-3) for a given timestep."""
        if t == 0:
            return 0
        return ((t - 1) % self.phases) + 1

    def get_round(self, t: int) -> int:
        """Get the round number for a given timestep."""
        if t == 0:
            return 0
        return ((t - 1) // self.phases) + 1

    def neighbors(self, v: int) -> List[int]:
        """Get neighbors of vertex v."""
        return list(self.graph.neighbors(v))

    def add_initial_conditions(self):
        """Add constraints for initial conditions."""
        # Root vertex starts burning
        self.solver.Add(self.burned[self.root, 0] == 1)

        # All other vertices start unburned
        for v in self.graph.nodes():
            if v != self.root:
                self.solver.Add(self.burned[v, 0] == 0)

        # No initial defenses
        for v in self.graph.nodes():
            self.solver.Add(self.defended[v, 0] == 0)

    def solve(self) -> bool:
        """
        Solve the ILP.

        Returns:
            bool: True if a solution was found, False otherwise
        """
        status = self.solver.Solve()
        return status == pywraplp.Solver.OPTIMAL

    def get_solution(self) -> Dict[str, List[List[int]]]:
        """
        Get the solution as a dictionary of burned and defended vertices at each timestep.

        Returns:
            Dict containing lists of burned and defended vertices at each timestep
        """
        if not self.solve():
            return None

        solution = {
            "burned": [[] for _ in range(self.T + 1)],
            "defended": [[] for _ in range(self.T + 1)],
        }

        for t in range(self.T + 1):
            for v in self.graph.nodes():
                if self.burned[v, t].solution_value() > 0.5:
                    solution["burned"][t].append(v)
                if self.defended[v, t].solution_value() > 0.5:
                    solution["defended"][t].append(v)

        return solution
