% setup
int: n;                                   % no vertices
int: r;                                   % root
int: d;     
int:steps=3;                              % no defenders per timestep
int: T = n*steps;
array[1..n, 1..n] of 0..1: graph;         

% variables
array[1..n, 0..T] of var 0..1: burned;    % burned[x,t] = if vertex x is burned at time t then 1 else 0
array[1..n, 0..T] of var 0..1: defended;  % defended[x,t] = if vertex x is defended at time t then 1 else 0

function var bool: time_step(var int: t, var int: s) = t mod steps == s;
function var int: real_time(var int: t) = t div steps;
function var int: prev_time(var int: t) = real_time(t-1);

% constraints

% setup
constraint forall(t in 0..steps)(burned[r,t]==1);
constraint forall(t in 0..steps)(forall(x in 1..n where x != r)(burned[x,t]==0));
constraint forall(t in 0..steps)(forall(x in 1..n)(defended[x,t]==0));


% peristence
constraint forall(t in 1..T)(forall(x in 1..n)(defended[x,t] >= defended[x,t-1]));
constraint forall(t in 1..T)(forall(x in 1..n)(burned[x,t] >= burned[x,t-1]));


% 4. single state per vertex
constraint forall(t in 1..T where time_step(t,3))(forall(x in 1..n)(burned[x,t] + defended[x,t] <= 1));
% constraint forall(t in 1..T)(forall(x in 1..n)(burned[x,t-1] + defended[x,t] <= 1));

% 3. force decision
constraint forall(t in 1..T where time_step(t,3))(forall(x in 1..n)(forall(y in 1..n where graph[x,y]==1)(burned[x,t] + defended[x,t] >= burned[y,prev_time(t)])));

% 6. at every timestep.1 make sure the total number of new defended nodes is at most equal to d
constraint forall(t in 1..T where time_step(t,3))(sum(x in 1..n)(defended[x,t] - defended[x,prev_time(t)]) <= d);

% 5. a vertex cannot spontaneously combust
% constraint forall(t in 1..T where time_step(t,3))(forall(x in 1..n)(sum(y in 1..n)(burned[y,real_time(t-1)]) >= burned[x,t]));
% no spontaneous combustion
constraint forall(t in 1..T where time_step(t,3))(forall(x in 1..n where burned[x,t-1] == 0)(burned[x,t] <= sum(y in 1..n where graph[x,y]==1)(burned[y,t-1])));


% constraint forall(t in 1..T where time_step(t,1))(
%   forall(x in 1..n)(
%     forall(y in 1..n where graph[x,y]==1)(
%       burned[x,real_time(t)] + defended[x,real_time(t)] >= burned[y,real_time(t-1)]
%     )
%   )
% );

% constraint forall(t in 1..T where t mod 3 == 1)(
%   forall(x in 1..n)(
%     forall(y in 1..n where graph[x,y]==1)(
%       burned[x,real_time(t-1)] + defended[x,real_time(t-1)] >= burned[y,t]
%     )));


% t=k
% direct defence
% at most d defenders per timestep
% constraint forall(t in 1..T where t mod 3 == 1)(
  % sum(x in 1..n)(defended[x,t] - defended[x,t-1]) <= d);


% t=2k
% fire spreads
% constraint forall(t in 1..T where t mod 3 == 2)(
  % forall(x in 1..n)(
  %   forall(y in 1..n where graph[x,y]==1)(burned[x,t] + defended[x,t] >= burned[y,t-1])));


% t=3k
% constraint forall(t in 1..T where t mod 3 == 0)(
  % forall(x in 1..n)(
    % forall(y in 1..n where graph[x,y]==1)(burned[x,t] + defended[x,t] >= defended[y,t-1])));
% defence spreads
% for every vertex that was defended in the previous round
% check all its neighbours
% for every one of its neighbours not currently burning -> defend them


% objective
solve minimize sum(x in 1..n)(burned[x,T]);

output [
  "Time: " ++ show(t div 3) ++ "." ++ show((t mod 3) + 1) ++ "\n" ++
  "Burned:\t\t" ++ show([if fix(burned[i,t]) == 1 then i else 0 endif | i in 1..n]) ++ "\n" ++
  "Defended:\t" ++ show([if fix(defended[i,t]) == 1 then i else 0 endif | i in 1..n]) ++ "\n\n"
  | t in 0..T where time_step(t,3)
];

