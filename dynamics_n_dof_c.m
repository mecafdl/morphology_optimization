clear all; clc;
n_dof = 2;
d = zeros(n_dof,1);
a = 0.5*ones(n_dof,1);
alpha = zeros(n_dof,1);


m{1} = 1.;
I{1} = diag([0.1,0.1,0.1]);
m{2} = 0.5;
I{2} = 0.5*diag([0.1,0.1,0.1]);

r = [0.25, 0, 0];


for i= 1:n_dof

        L(i) = Link('revolute', 'd', d(i), 'a', a(i), 'alpha', alpha(i), 'm', m{i}, 'I', I{i}, 'r', r);
end
robot = SerialLink(L, 'name', 'robot arm', 'gravity', [0 0 -9.81]);

%%

rng(0)
RR = randrot(3,2);
period = 5;


%% Symbolic vectors
q = sym('q', [1,n_dof], 'real');
dq = sym('dq', [1,n_dof], 'real');
tau = sym('tau', [1,n_dof], 'real');


%% Extract kinematics terms
J = robot.jacob0(q); %Jacobian in world coordinate frame
Jt = J(1:3,:);
%% Extract dynamics terms

% 
M = simplify(robot.inertia(q));
 C = robot.coriolis(q,dq);
%   G = robot.gravload(q);

a = sym('a', [2,1], 'real');
b = sym('b', [2,1], 'real');
c = sym('c', [2,1], 'real');

a = [0.075; 0.075; 0.075];
b = [0.2; 0.2; 0.2];
c = [0.075; 0.075; 0.075];

r1_ori = [0.;
    a(1);
    0];

r2_ori = [0.;
    -a(1);
    0];

r1_ins = [b(1)* cos(q(1)) +  c(1) *cos(q(1) + pi/2) ;
             b(1)* sin(q(1)) +  c(1) *sin(q(1)+ pi/2);
    0];

r2_ins = [b(1)* cos(q(1)) + c(1) *cos(q(1) - pi/2) ;
             b(1)* sin(q(1)) + c(1) *sin(q(1) - pi/2);
    0];



r3_ori = [ a(2) *cos(q(1)+ pi/2)  +  0.5 * cos(q(1));
             a(2) *sin(q(1)+ pi/2)   + 0.5*  sin(q(1));
             0];
r4_ori = [ a(2) *cos(q(1) - pi/2)+ 0.5 * cos(q(1));
              a(2) *sin(q(1)- pi/2)+  0.5 * sin(q(1));
              0];       

r3_ins = [0.5 * cos(q(1)) + b(2) *  cos(q(1)+ q(2)) + c(2) *cos(pi/2+q(1) + q(2)) ;
            0.5 * sin(q(1)) + b(2) *  sin(q(1)+ q(2)) + c(2) *sin(pi/2+q(1) + q(2));
            0];
r4_ins = [0.5 * cos(q(1)) + b(2) *  cos(q(1)+ q(2)) + c(2) *cos(-pi/2+q(1) + q(2)) ;
             0.5 * sin(q(1)) + b(2) *  sin(q(1)+ q(2)) + c(2) *sin(-pi/2+q(1) + q(2));
             0];

r5_ori = [0.;
    a(3);
    0];

r6_ori = [0.;
    -a(3);
    0];
         
r5_ins = [0.5 * cos(q(1)) + b(3) *  cos(q(1)+ q(2)) + c(3) *cos(pi/2+q(1) + q(2)) ;
            0.5 * sin(q(1)) + b(3) *  sin(q(1)+ q(2)) + c(3) *sin(pi/2+q(1) + q(2));
            0];
r6_ins = [0.5 * cos(q(1)) + b(3) *  cos(q(1)+ q(2)) + c(3) *cos(-pi/2+q(1) + q(2)) ;
             0.5 * sin(q(1)) + b(3) *  sin(q(1)+ q(2)) + c(3) *sin(-pi/2+q(1) + q(2));
             0];         
         
     
r1  = r1_ins-r1_ori;
r2  = r2_ins-r2_ori;
r3  = r3_ins-r3_ori;
r4  = r4_ins-r4_ori;
r5  = r5_ins-r5_ori;
r6  = r6_ins-r6_ori;


J1 = jacobian(r1, q');
J2 = jacobian(r2, q');
J3 = jacobian(r3, q');
J4 = jacobian(r4, q');
J5 = jacobian(r5, q');
J6 = jacobian(r6, q');


l1 = simplify(norm(r1));
l2 = simplify(norm(r2));
l3 = simplify(norm(r3));
l4 = simplify(norm(r4));
l5 = simplify(norm(r5));
l6 = simplify(norm(r6));

unilat  = @(disp)(tanh(1000*disp) - 0.5) +0.5;

% disp1 = k(1)*((l0(1) - l1)) ;%*unilat((l0(1) - l1));
% disp2 = k(2)*((l0(2) - l2)) ;%* unilat((l0(2) - l2));
% disp3 = k(3)*((l0(3) - l3)) ;%* unilat((l0(3) - l3));
% disp4 = k(4)*((l0(4) - l4)) ;%* unilat((l0(4) - l4));

F = sym('F', [6,1], 'real');


% F1 = -(r1_plus /l1)  *F(1)* unilat((l0(1) - l1));
% F2 = -(r1_minus/l2)  *F(2)* unilat((l0(2) - l2));
% F3 = -(r2_plus/l3)   *F(3)* unilat((l0(3) - l3));
% F4 = -(r2_minus /l4) *F(4)* unilat((l0(4) - l4));

F1 = -(r1/l1)  *F(1);%* unilat((l0(1) - l1));
F2 = -(r2/l2)  *F(2);%* unilat((l0(2) - l2));
F3 = -(r3/l3)  *F(3);%* unilat((l0(3) - l3));
F4 = -(r4/l4)  *F(4);%* unilat((l0(4) - l4));

  J_lever = [ J1' * (r1 /l1),...
       J2' * (r2 /l2),...
      J3'  * (r3 /l3) , ...
      J4'  *  (r4 /l4),...
      J5'  * (r5 /l5) , ...
      J6'  *  (r6 /l6) ];
  
tau = -J_lever* F;
J_lever = simplify(J_lever);
tau = simplify(tau);  

l_vect = [l1;l2;l3;l4;l5;l6];
%
%% Create matlab functions
matlabFunction(tau,'File','get_tau_F', 'Vars', {q', F});
matlabFunction(l_vect,'File','get_l', 'Vars', {q'});


matlabFunction((pinv(J_lever)),'File','get_J_lever_pinv', 'Vars', {q'});

% matlabFunction(l_vect,'File','get_l', 'Vars', {q' ,a, b,c});
 
% matlabFunction(l1,'File','get_l1', 'Vars', {q' ,a, b});
%  matlabFunction(l2,'File','get_l2', 'Vars', {q', a, b});
%    matlabFunction(l3,'File','get_l3', 'Vars', {q', a, b});
%   matlabFunction(l4,'File','get_l4', 'Vars', {q', a, b});

matlabFunction(M,'File','get_M', 'Vars', {q'});
matlabFunction(C,'File','get_C', 'Vars', {q',dq'});
% matlabFunction(G','File','get_G', 'Vars', {q'});
%matlabFunction(ddq','File','inv_dyn', 'Vars', {q',dq', tau'});
  %matlabFunction(Jt,'File','get_EE_velocity_jacobian', 'Vars', {q'});


%% Show Robot
%q = zeros(1,n_dof);
%q(1) = 0.2;

%  W = [xmn, xmx ymn ymx zmn zmx]
 %W = [-1 1 -1 1 -1 1];
%obot.plot(q, 'workspace', W);

defl = -1:0.01:1;
plot(defl, exp_spring(defl))

function F = exp_spring(defl)
F = exp(defl) - exp(-defl);
end
