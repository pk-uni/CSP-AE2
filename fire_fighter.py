from pulp import *

lp = LpProblem("Firefighter", LpMinimize)

x = LpVariable.dicts("x", range(10), cat="Binary")

print(x)
