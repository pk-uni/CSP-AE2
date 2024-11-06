#  Infectious Firefighter Problem

We are given a graph $G = (V,E)$ with a root $r$ denoting the starting point of the fire. Then, on at each round the following happens in this order: 

1. An unburning vertex is chosen to be defended
2. All undefended vertices that are adjacent to a burning vertex catch fire
3. All unburning and undefended vertices that are adjacent to a defended vertex become defended.

The process ends when there are no burning vertices adjacent to undefended unburning vertices - that is, the fire can spread no further. The number of vertices saved is the number of unburning vertices at the end of the process. 

The main difference between this process and the standard firefighter problem is that here the defence spreads.


## Solution

Let $t$ represent an absolute timestep. Given that a round requires 3 timesteps, we can define a function for determining the phase we are in:
```
phase t = ((t-1) mod 3) + 1
```

We have a way to determine what phase we're in. We now need to determine what action takes place at each phase in a round:

- phase 1: Direct defense placement
- phase 2: Fire spreading
- phase 3: Defence spreading

### Variables
- $b_{x,t} =$ `if` vertex $x$ is burning at $t=0$ `then` $1$ `else` $0$
- $d_{x,t} =$ `if` vertex $x$ is defended at $t=0$ `then` $1$ `else` $0$

### Objective
Minimise the number of burned vertices at the final time-step:

$$\text{Min } \sum_{x \in V} b_{x,T}$$

### Constraints

1. **Initial Conditions**
   $$b_{x,0} = \begin{cases} 1 & \text{if } x = r \\ 0 & \text{otherwise} \end{cases}$$
   $$d_{x,0} = 0 \quad \forall x \in V$$

2. **Permanence**
   $$b_{x,t} \geq b_{x,t-1} \quad \forall x \in V, t \geq 1$$
   $$d_{x,t} \geq d_{x,t-1} \quad \forall x \in V, t \geq 1$$

3. **State Exclusivity**
   $$b_{x,t} + d_{x,t} \leq 1 \quad \forall x \in V, t \geq 0$$

4. **Defense Budget** (only during phase 1)
   $$\sum_{x \in V} (d_{x,3t} - d_{x,3t-1}) \leq d \quad \forall t \geq 0$$

5. **Fire Spread** (during phase 2)
   $$b_{x,3t+1} + d_{x,3t+1} \geq b_{y,3t} \quad \forall x \in V, y \in N(x), t \geq 0$$

6. **No Spontaneous Combustion** (during phase 2)
   $$b_{x,3t+1} \leq \sum_{y \in N(x)} b_{y,3t} \quad \forall x \in V \text{ where } b_{x,3t} = 0, t \geq 0$$

7. **Defence Spread** (during phase 3)
   $$d_{y,3t+2} \wedge \neg b_{x,3t+2} \rightarrow d_{x,3t+2} = 1 \quad \forall x \in V, y \in N(x), t \geq 0$$

8. **Phase Consistency** (states only change during their respective phases)
   - During phase 1: Only defence can change
   $$b_{x,3t} = b_{x,3t-1} \quad \forall x \in V, t \geq 1$$
   
   - During phase 2: Only fire can spread
   $$d_{x,3t+1} = d_{x,3t} \quad \forall x \in V, t \geq 0$$
   
   - During phase 3: Only defence can spread
   $$b_{x,3t+2} = b_{x,3t+1} \quad \forall x \in V, t \geq 0$$
   peco