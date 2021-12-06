# QuadrupedBalance.jl


## How to run
At the Package directory, start jupyter notebook by running `jupyter notebook`. All experiments are under the `notebooks` directory, 

### Finding Equilibrium pose
 Run `Equilibrium finder (ipopt).ipynb`. This should save a file called `ipopt_eq_poimt.toml` that contains `x_eq`, `u_eq`, and `\lambda_eq` for eq state, control, and ground reaction forces. 

### Finding LQR Gains 
 Run `MaximalCoordinateLqr.ipynb`. This should output a delimited file called `maximal_lqr_gain.txt` that stores the `36x12` error feedback gain matrix.

### Simulation 
 Run `BalanceSim.ipynb` for a constrained semi-implicit euler integration simulation. It will load in `ipopt_eq_poimt.toml` and `maximal_lqr_gain.txt`.  

### Function docs 
For some documentation and examples on some of the core functions, run `UsefulFunctions.ipynb`.