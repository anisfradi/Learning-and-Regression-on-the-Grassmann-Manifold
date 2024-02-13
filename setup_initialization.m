


% Setup 
params.Ys=Q;
params.ts=time;
params.alpha = 0;               % alpha for balancing the prior knowledge of the slope of the geodesic
params.sigmaSqr = 1;            % sigma square for balancing the mismatching term
params.useODE45 = true;         % chose ODE45 for solving the PDEs, more stable
params.useRansac=0;
params.wi=ones(N,1);
