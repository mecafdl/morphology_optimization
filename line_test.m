clear
clc
close all

% Add path
addpath('/home/dennis/casadi_matlab')
import casadi.*


%% Define parameters


%% Optimization
disp('Creating training data...');

in = -1:0.01:1;
out = 5.*in.^2;

% Create basis function parameters
n_weights = 10;
a = SX.sym('a',n_weights,1);
b = SX.sym('b',n_weights,1);
c = SX.sym('c',n_weights,1);

x = SX.sym('x',1,1);

y = rbf_norm(x, a,b,c);
f = Function('f', {x,a,b,c}, {y});

% Create OCP
opti = Opti();

% Create decision variables for initial step
p = opti.variable(3*n_weights,1);
opti.set_initial(p, ones(3*n_weights,1));

cost = 0.;
for i=1:max(size(in))
    cost = cost + (out(i) - f(in(i), p(1:n_weights),p(n_weights+1:2*n_weights),p(2*n_weights+1:3*n_weights))).^2;
end


opti.minimize(cost );

casadi_opts = struct();
solver_opts = struct('max_iter', 10000, 'linear_solver', 'ma86');
opti.solver('ipopt',casadi_opts, solver_opts); % set numerical backend

opti.solve();

%% Evaluation
a_sol = opti.debug.value(A);
b_sol = opti.debug.value(B);
c_sol = opti.debug.value(C);


plot(in, out);
hold on
plot(in, rbf_norm(in,a_sol,b_sol,c_sol))
