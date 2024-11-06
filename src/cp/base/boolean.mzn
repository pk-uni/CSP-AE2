% setup
int: n;                                     % no vertices
int: T=2*n;                                 % no phases
int: r;                                     % root
int: d;                                     % no defenders per timestep
array[1..n, 1..n] of bool: graph;           % using bool instead of 0..1 for graph

% variables
array[1..n, 0..T] of var bool: burned;      % burned[x,t] = x == burned at t
array[1..n, 0..T] of var bool: defended;    % defended[x,t] = x == defended at t

% helper functions for boolean comparisons
function var bool: bool_geq(var bool: a, var bool: b) = 
    a \/ not b;  % equivalent to a >= b

function var bool: bool_leq(var bool: a, var bool: b) = 
    not a \/ b;  % equivalent to a <= b

% helper function to count true values in an array
function var int: count_true(array[int] of var bool: arr) =
    sum([bool2int(x) | x in arr]);

% helper function to get burning neighbors
function var int: burning_neighbors(int: vertex, var int: time) =
    count_true([burned[y,time] | y in 1..n where graph[vertex,y]]);

% constraints

% initial conditions
constraint forall(x in 1..n)(not defended[x,0]);
constraint burned[r,0];
constraint forall(x in 1..n where x != r)(not burned[x,0]);

% permanence
constraint forall(t in 1..T)(forall(x in 1..n)(bool_geq(defended[x,t], defended[x,t-1])));
constraint forall(t in 1..T)(forall(x in 1..n)(bool_geq(burned[x,t], burned[x,t-1])));

% fire spreads
constraint forall(t in 1..T)(forall(x in 1..n)(forall(y in 1..n where graph[x,y])(bool_geq(burned[x,t] \/ defended[x,t], burned[y,t-1]))));

% single state per vertex
constraint forall(t in 1..T)(forall(x in 1..n)(not (burned[x,t] /\ defended[x,t])));

% no spontaneous combustion
constraint forall(t in 1..T)(forall(x in 1..n where not burned[x,t-1])(bool_leq(burned[x,t], burning_neighbors(x,t-1) > 0)));

% at most d new defenders per timestep
constraint forall(t in 1..T)(count_true([defended[x,t] /\ not defended[x,t-1] | x in 1..n]) <= d);

% objective
solve minimize count_true([burned[x,T] | x in 1..n]);

output [
  "Time " ++ show(t) ++ ":\n" ++
  "Burned: " ++ show([if fix(burned[i,t]) then i else 0 endif | i in 1..n]) ++ "\n" ++
  "Defended: " ++ show([if fix(defended[i,t]) then i else 0 endif | i in 1..n]) ++ "\n"
  | t in 0..T
];