from utils import parse_networkx
import cp.main as cp
import ilp.main as ilp


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

    (solver, ctx) = ilp.generate_solver(n, r, graph)
    num_saved = ilp.solve(solver, ctx, n, n * 3)

    return {"num_saved": num_saved}


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
    n, r, graph = parse_networkx(instance_graph, start_node)

    instance = cp.generate_solver(n, r, graph)
    num_saved = cp.solve(instance, n)

    return {"num_saved": num_saved}
