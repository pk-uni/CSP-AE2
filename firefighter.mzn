% setup
int: n;                                    % number of vertices
int: T;                                    % number of timesteps
int: r;                                    % root vertex where fire starts
int: d;                                    % number of defenders per timestep
array[1..n, 1..n] of 0..1: graph;         % adjacency matrix

% variables
array[1..n, 0..T] of var 0..1: burned;    % burned[x,t] = if vertex x is burned at time t then 1 else 0
array[1..n, 0..T] of var 0..1: defended;  % defended[x,t] = if vertex x is defended at time t then 1 else 0

% constraints

% initial conditions
constraint forall(x in 1..n)(defended[x,0]==0);
constraint burned[r,0]==1;
constraint forall(x in 1..n where x != r)(burned[x,0]==0);

% permanence
constraint forall(t in 1..T)(forall(x in 1..n)(defended[x,t] >= defended[x,t-1]));
constraint forall(t in 1..T)(forall(x in 1..n)(burned[x,t] >= burned[x,t-1]));

% fire spreads
constraint forall(t in 1..T)(forall(x in 1..n)(forall(y in 1..n where graph[x,y]==1)(burned[x,t] + defended[x,t] >= burned[y,t-1])));

% single state per vertex
constraint forall(t in 1..T)(forall(x in 1..n)(burned[x,t] + defended[x,t] <= 1));
constraint forall(t in 1..T)(forall(x in 1..n)(burned[x,t-1] + defended[x,t] <= 1));

% no spontaneous combustion
constraint forall(t in 1..T)(forall(x in 1..n)(burned[x,t] >= sum(y in 1..n where graph[x,y]==1)(burned[y,t-1])));

% at most d defenders
constraint forall(t in 1..T)(sum(x in 1..n)(defended[x,t] - defended[x,t-1]) <= d);

% objective
solve minimize sum(x in 1..n)(burned[x,T]);

% output
output [
  "Time " ++ show(t) ++ ":\n" ++
  "Burned: " ++ show([if fix(burned[i,t]) == 1 then i else 0 endif | i in 1..n]) ++ "\n" ++
  "Defended: " ++ show([if fix(defended[i,t]) == 1 then i else 0 endif | i in 1..n]) ++ "\n"
  | t in 0..T
];