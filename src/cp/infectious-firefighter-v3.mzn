% setup
int: n;                                   % no vertices
int: r;                                   % root
int: d;     
int: phases=3;                              % no defenders per timestep
int: T = n*phases;
array[1..n, 1..n] of 0..1: graph;         

% variables
array[1..n, 0..T] of var 0..1: burned;    % burned[x,t] = if vertex x is burned at time t then 1 else 0
array[1..n, 0..T] of var 0..1: defended;  % defended[x,t] = if vertex x is defended at time t then 1 else 0

function var bool: phase(var int: t, var int: p) = (t-1) mod phases == (p-1) mod phases;
function var int: real_time(var int: t) = ((t-1) div phases)+1;
function var int: prev_time(var int: t) = real_time(t)-1;


% constraints

% setup
constraint (burned[r,0]==1);
constraint (forall(x in 1..n where x != r)(burned[x,0]==0));
constraint (forall(x in 1..n)(defended[x,0]==0));

% state exclusivity
constraint forall(t in 0..T)(
  forall(x in 1..n)(
    burned[x,t] + defended[x,t] <= 1
  )
);

% permanence
constraint forall(t in 1..T)(
  forall(x in 1..n)(
    burned[x,t-1] <= burned[x,t] /\
    defended[x,t-1] <= defended[x,t]
  )
);


% phase 1: defense budget
constraint forall(t in 1..T where phase(t, 1))(
    sum(x in 1..n)(defended[x,t] - defended[x,t-1]) <= d
);

% phase 1: only defense placement can change
constraint forall(t in 1..T where phase(t, 1))(
    forall(x in 1..n)(
        burned[x,t] = burned[x,t-1]
    )
);


% phase 2: fire spread
constraint forall(t in 1..T where phase(t, 2))(
    forall(x in 1..n)(
        forall(y in 1..n where graph[x,y]==1)(
            burned[x,t] + defended[x,t] >= burned[y,t-1]
        )
    )
);

% phase 2: only fire can spread
constraint forall(t in 1..T where phase(t, 2))(
    forall(x in 1..n)(
        defended[x,t] = defended[x,t-1]
    )
);


% phase 3: defence spread
constraint forall(t in 1..T where phase(t, 3))(
    forall(x in 1..n)(
        forall(y in 1..n where graph[x,y]==1)(
          defended[y,t-1]==1 /\ burned[x,t-1]==0 -> defended[x,t]==1
        )
    )
);

% phase 3: only defense can spread
constraint forall(t in 1..T where phase(t, 3))(
    forall(x in 1..n)(
        burned[x,t] = burned[x,t-1]
    )
);

% objective
solve minimize sum(x in 1..n)(burned[x,T]);

output [
  "t: " ++ show(t) ++ "\n" ++
  "Round: " ++ show(real_time(t)) ++ "." ++ show(((t-1) mod 3)+1) ++ "\n" ++
  "Burned:\t\t" ++ show([if fix(burned[i,t]) == 1 then i else 0 endif | i in 1..n]) ++ "\n" ++
  "Defended:\t" ++ show([if fix(defended[i,t]) == 1 then i else 0 endif | i in 1..n]) ++ "\n\n"
  | t in 0..T
];

