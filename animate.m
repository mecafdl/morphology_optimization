q = out.q.Data;
dq = out.dq.Data;

%tau_c = out.tau_c.Data;
%disp = out.displacement.Data;

% q(1:10000,:) = [];
%  tau_c(1:10000,:) = [];
 %disp(1:10000,:) = [];
% 
 %plot(disp(:,2), tau_c(:,2))

%plot(q(:,1), tau_c(:,1))
%plot(q(:,2), tau_c(:,2))

%plot(dq(:,2), q(:,2))

%%


q(1:2:end,:) = [];
q(1:2:end,:) = [];
q(1:2:end,:) = [];
q(1:2:end,:) = [];


%q(2:2:end,:) = []; 
%q(2:2:end,:) = []; 
%q(2:2:end,:) = []; 
%q(2:2:end,:) = []; 
%q(2:2:end,:) = []; 

%q(1:30000,:) = [];


full_time = max(size(q)) / 1e3;

n_dof = 2;
d = zeros(n_dof,1);
a = 0.4*ones(n_dof,1);
alpha = zeros(n_dof,1);


m = 1.;
I = diag([0.1,0.1,0.2]);


r = [0.05, 0, 0.];


for i= 1:n_dof

        L(i) = Link('revolute', 'd', d(i), 'a', a(i), 'alpha', alpha(i), 'm', m, 'I', I, 'r', r);
end
robot = SerialLink(L, 'name', 'ant. system', 'gravity', [0 0 -9.81]);


%% Animate the robot motion
% Set frame rate
clc
close all
frame_rate = max(size(q))/10; % duration = 10 s

 W = [-0.25 1.75 -1 1 -0.2 1];
  W = [-0.1 0.9 -0.4 0.6 -0.05 0.4];

% Plot robot
figure();
%robot.plot(q)
robot.plot(q,'loop',  'notiles', 'workspace', W, 'fps', frame_rate,  'trail', 'b--');

%robot.plot(q, 'loop', 'notiles', 'workspace',  W, 'fps', frame_rate, 'trail', 'b--', 'movie', '/home/dennis/matlab-tail/R2_Controlled/movie');
%robot.plot([q; zeros(1, size(q,2))]', 'loop', 'view', 'top', 'notiles', 'workspace',  [-1.4*h_joint1 1.4*h_joint1 -1.4*h_joint1 1.4*h_joint1 -1.4*h_joint1 1.4*h_joint1], 'fps', frame_rate, 'trail', 'b--', 'movie', 'animation.gif');

%ffmpeg -f image2 -i %4d.png test.avi
%ffmpeg -i test.avi -filter:v "setpts=0.4*PTS" -c:v mpeg4 -q:v 2  -an       out.avi
