clear
clc
close all

% Add path
addpath('/home/dennis/casadi_matlab')
import casadi.*


%% Define parameters


%% Optimization
disp('Creating training data...');

N_sample = 10;
N_grid   = N_sample*N_sample;
samples = linspace(-pi/6, pi/6, N_sample);
[q1,q2] = meshgrid(samples);

q1 = reshape(q1, N_grid,1);
q2 = reshape(q2, N_grid,1);

q = [q1,q2];


tau_cc = zeros(N_grid,2);
for i = 1:N_grid
    tau_cc(i,:) = -75.*get_M(q(i,:)') *q(i,:)';
end

%%

% Create basis function parameters
n_w = 5;
a = SX.sym('a',n_w,1);
b = SX.sym('b',n_w,1);
c = SX.sym('c',n_w,1);

x = SX.sym('x',1,1);

y = rbf(x, a,b,c);
f = Function('f', {x,a,b,c}, {y});

% Create OCPf
opti = Opti();

% Create decision variables for initial step
p1 = opti.variable(n_w*3,1);
p2 = opti.variable(n_w*3,1);
l0 =get_l(zeros(2,1));

%p3 = opti.variable(n_w*3,1);
%p4 = opti.variable(n_w*3,1);

% opti.set_initial(p1, 0.01*ones(n_w*3,1));
 %opti.set_initial(p2, 0.01*ones(n_w*3,1));
% opti.set_initial(p3, 0.01*ones(n_w*3,1));
% opti.set_initial(p4, 0.01*ones(n_w*3,1));

delta_l_debug = {};
F_debug = {};
cost = 0.;
for i=1:N_grid
    l = get_l(q(i,:)');
    Delta_l = l0 - l;
    
    F= [f(Delta_l(1), p1(1:n_w),p1(n_w+1:2*n_w),p1(2*n_w+1:3*n_w));
        f(Delta_l(2), p1(1:n_w),p1(n_w+1:2*n_w),p1(2*n_w+1:3*n_w));
        f(Delta_l(3), p2(1:n_w),p2(n_w+1:2*n_w),p2(2*n_w+1:3*n_w));
        f(Delta_l(4), p2(1:n_w),p2(n_w+1:2*n_w),p2(2*n_w+1:3*n_w))];
    
F_cc = -get_J_lever_pinv(q(i,:)') * tau_cc(i,:)';

    
  %  cost = cost + (tau_cc(i,1)' - tau_ant(1)).^2  + (tau_cc(i,2)' - tau_ant(2)).^2;

    cost = cost + (F_cc(1) - F(1)).^2  + ...
        (F_cc(2) - F(2)).^2  + ...
        (F_cc(3) - F(3)).^2  + ...
        (F_cc(4) - F(4)).^2  ;


delta_l_debug{end+1} = Delta_l;
F_debug{end+1} = F;

end
% 
% F0= [f(0, p1(1:n_w),p1(n_w+1:2*n_w),p1(2*n_w+1:3*n_w));
%         f(0, p1(1:n_w),p1(n_w+1:2*n_w),p1(2*n_w+1:3*n_w));
%         f(0, p2(1:n_w),p2(n_w+1:2*n_w),p2(2*n_w+1:3*n_w));
%         f(0, p2(1:n_w),p2(n_w+1:2*n_w),p2(2*n_w+1:3*n_w))];
% 
% opti.subject_to(F0 <= 0.01*ones(4,1));
% opti.subject_to(F0 >= -0.01*ones(4,1));
% 

  delta_l_debug =  [delta_l_debug{:}];
  F_debug =  [F_debug{:}];

opti.minimize(cost );

casadi_opts = struct();
solver_opts = struct('max_iter', 10000, 'linear_solver', 'ma97');
opti.solver('ipopt',casadi_opts, solver_opts); % set numerical backend

opti.solve();

%% Evaluation

%tau_ant = opti.debug.value(tau_ant_debug)';
delta_l = opti.debug.value(delta_l_debug)';
F = opti.debug.value(F_debug)';
l0 = opti.debug.value(l0);


for i = 1:N_grid
        tau_ant(i,:) = -get_tau_F(q(i,:)',F(i,:)');

end


%plot(delta_l(:,3), F(:,3))

x = -0.025:0.001:0.025
a_sol = opti.debug.value(p1(1:n_w));
 b_sol = opti.debug.value(p1(n_w+1:2*n_w));
 c_sol = opti.debug.value(p1(2*n_w+1:3*n_w));
 
 plot(x, rbf(x,a_sol,b_sol,c_sol))
 