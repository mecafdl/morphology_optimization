function [Gamma, dGamma] = compute_geodesics(q0, N_geo, N_int)

tspan = 0:0.01:50;
%Create unit vectors
R = @(angle)[cos(angle) -sin(angle);
    sin(angle)  cos(angle)];
dQ0 = {};
for angle = 0.: (2*pi)/N_geo : 2*pi
    dq0 = R(angle)*[1;0];
    dQ0{end+1} = dq0;
end

%%

Xinterp = {};
Gamma = {};
dGamma= {};

colors = hsv(N_geo);

s = 0:1/N_int:1;
for i=1:length(dQ0)-1
    dq0 = dQ0{i};
x0 = [q0;dq0];

% Optionen für ode45
options = odeset('Events', @event_function);

% Integrieren mit ode45
[t, x, te, ~, ~] = ode45(@ode_function, tspan, x0, options);

t_normed = t/te; % Normed time vector t_norm: 0->1
x_interp = interp1(t_normed, x, s);


gamma = x_interp(:,1:2);
dgamma = x_interp(:,3:4);
%dgamma_orth = [null(dgamma'), null(dgamma')];


Gamma{end+1} = gamma;
dGamma{end+1} = dgamma;

%Gamma_orth{end+1} = [null(dgamma'), null(dgamma')];

% for j = 1:length(gamma)
%     gamma_curr = gamma(j,:)';
%     dgamma_orth = [null(gamma_curr'), null(gamma_curr')];
%     rhs = get_M(gamma_curr) \ dgamma_orth;
% end

Xinterp{end+1} = x_interp;
plot(x_interp(:,1), x_interp(:,2), 'color', colors(i,:));
hold on
end
% Plot Ergebnisse
end

function dxdt = ode_function(t, x)

q = x(1:2);
dq = x(3:4);

ddq = get_M(q) \ (-get_C(q,dq) *dq);
dxdt  = [dq; ddq];
end

function [value, isterminal, direction] = event_function(t, x)
% Ereignisfunktion: Wenn
q = x(1:2);
q_edge = [ 1.5;1.5];

value = [q - q_edge; q + q_edge] ;
isterminal = [1;1;1;1]; % Stoppt die Integration
direction = [0;0;0;0];
end