% setup
int: n;                                   % no vertices
int: T;                                   % no timesteps
int: r;                                   % root
int: d;                                   % no defenders per timestep
array[1..n, 1..n] of 0..1: graph;         

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


function var bool: defended_before(var int: x, var int: t) = 
    defended[x,t-1]==1;

function var bool: defended_just_now(var int: x, var int: t) = 
    defended[x,t-1]==0 /\ defended[x,t]==1;


% defence spreads
constraint forall(t in 1..T)(
  forall(x in 1..n where defended[x,t-1]==1)(
    forall(y in 1..n where graph[x,y]==1 /\ (burned[y,t-1]==0))(defended[y,t]==1)));

% for every vertex that was defended in the previous round
% check all its neighbours
% for every one of its neighbours not currently burning -> defend them

% single state per vertex
constraint forall(t in 1..T)(forall(x in 1..n)(burned[x,t] + defended[x,t] <= 1));
% constraint forall(t in 1..T)(forall(x in 1..n)(burned[x,t-1] + defended[x,t] <= 1));

% no spontaneous combustion
constraint forall(t in 1..T)(forall(x in 1..n where burned[x,t-1] == 0)(burned[x,t] <= sum(y in 1..n where graph[x,y]==1)(burned[y,t-1])));

% at most d defenders per timestep
constraint forall(t in 1..T)(sum(x in 1..n)(defended[x,t] - defended[x,t-1]) <= d);

% objective
solve minimize sum(x in 1..n)(burned[x,T]);

output [
  "Time " ++ show(t) ++ ":\n" ++
  "Burned: " ++ show([if fix(burned[i,t]) == 1 then i else 0 endif | i in 1..n]) ++ "\n" ++
  "Defended: " ++ show([if fix(defended[i,t]) == 1 then i else 0 endif | i in 1..n]) ++ "\n"
  | t in 0..T
];