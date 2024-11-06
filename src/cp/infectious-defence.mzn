% setup
int: n;                                   % no vertices
int: r;                                   % root
int: d;
int: phases=3;
int: T = n*phases;
array[1..n, 1..n] of 0..1: graph;         

% variables
array[1..n, 0..T] of var 0..1: burned;    % burned[x,t] = if x == burned at t then 1 else 0
array[1..n, 0..T] of var 0..1: defended;  % defended[x,t] = if x == defended at t then 1 else 0

function var int: phase(var int: t) = ((t-1) mod phases) + 1;
function var int: round(var int: t) = ((t-1) div phases) + 1;

% constraints

% initial conditions
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
constraint forall(t in 1..T where phase(t) == 1)(
    sum(x in 1..n)(defended[x,t] - defended[x,t-1]) <= d
);

% phase 1: only defense placement can change
constraint forall(t in 1..T where phase(t) == 1)(
    forall(x in 1..n)(
        burned[x,t] = burned[x,t-1]
    )
);


% phase 2: fire spread
constraint forall(t in 1..T where phase(t) == 2)(
    forall(x in 1..n)(
        forall(y in 1..n where graph[x,y]==1)(
            burned[x,t] + defended[x,t] >= burned[y,t-1]
        )
    )
);

% phase 2: only fire can spread
constraint forall(t in 1..T where phase(t) == 2)(
    forall(x in 1..n)(
        defended[x,t] = defended[x,t-1]
    )
);


% phase 3: defence spread
constraint forall(t in 1..T where phase(t) == 3)(
    forall(x in 1..n)(
        forall(y in 1..n where graph[x,y]==1)(
          defended[y,t-1]==1 /\ burned[x,t-1]==0 -> defended[x,t]==1
        )
    )
);

% phase 3: only defense can spread
constraint forall(t in 1..T where phase(t) == 3)(
    forall(x in 1..n)(
        burned[x,t] = burned[x,t-1]
    )
);

% objective
solve minimize sum(x in 1..n)(burned[x,T]);

output [
  "t: " ++ show(t) ++ "\n" ++
  "Round: " ++ show(round(t)) ++ "." ++ show(phase(t)) ++ "\n" ++
  "Burned:\t\t" ++ show([if fix(burned[i,t]) == 1 then i else 0 endif | i in 1..n]) ++ "\n" ++
  "Defended:\t" ++ show([if fix(defended[i,t]) == 1 then i else 0 endif | i in 1..n]) ++ "\n\n"
  | t in 0..T
];

