
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%                                                         %%%%%%%%%
%%%%%%%%%       ANALISI E SIMULAZIONE DI SISTEMI AEROSPAZIALI     %%%%%%%%%
%%%%%%%%%                       A.A. 2024-2025                    %%%%%%%%%
%%%%%%%%%                                                         %%%%%%%%%
%%%%%%%%%                          GRUPPO 26                      %%%%%%%%%
%%%%%%%%%                        LABORATORIO 1                    %%%%%%%%%
%%%%%%%%%                                                         %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
clear
clc
close all

%% CODE

% Task 1
MPend.n = 10;            % number of modes
MPend.xi = 0.003;        % modal damping ratio
MPend.d = 1;             % tank height [m]
MPend.h = 1;             % liquid height [m]
MPend.rho = 1000;        % water density [kg/m^3]
MPend.g = 9.81;          % gravity [m/s^2]

MPend = parameters(MPend);  % parameters calculation

% Task 3

MPend.u_type = 'impulse';    % input type: step
MPend.t0 = 0;                % initial simulation time
MPend.tf = 10;               % final simulation time
MPend.n = 1;                 % number of modes

% built-in impulse response
[y_imp_built_in,t_imp, MPend.sys] = slosh_response(MPend);

figure;
subplot(2,2,1);
plot(t_imp, y_imp_built_in,'LineWidth',2);
xlabel('t [s]');
ylabel('F_x(t) [N]');
grid;
title('Built-in Impulse Response (n = 1)');

% analytical impulse response

A = MPend.sys.A;
B = MPend.sys.B;
C = MPend.sys.C;
D = MPend.sys.D;

t_imp_analytical = linspace(0,10,300);
N = length(t_imp_analytical);
y_imp_analytical = zeros(N,1);

y_imp_analytical(1) = 0;
for i = 2:N
    y_imp_analytical(i) = C*expm(A.*t_imp_analytical(i))*B + D*MPend_input(t_imp_analytical(i),MPend);  
end

subplot(2,2,2);
plot(t_imp_analytical,y_imp_analytical,'LineWidth',2);
xlabel('t [s]');
ylabel('F_x(t) [N]');
grid;
title('Analytical Impulse Response (n = 1)');

% analytical impulse response in modal form

[V,lambda] = eig(A);
inv_V = V \ eye(2*MPend.n); 
y_imp_modal_analytical = zeros(N,1);

y_imp_modal_analytical(1) = 0;
for i = 2:N
    e_lambda = diag(exp(diag(lambda).*t_imp_analytical(i)));
    y_imp_modal_analytical(i) = C*V*e_lambda*inv_V*B + D.*MPend_input(t_imp_analytical(i),MPend);   
end

subplot(2,2,3);
plot(t_imp_analytical,y_imp_modal_analytical,'LineWidth',2);
xlabel('t [s]');
ylabel('F_x(t) [N]');
grid;
title('Modal Analytical Impulse Response (n = 1)');

% numerical integration technique

[t_imp_ode,y_imp_ode] = MPend_Ode_Calc(MPend);

subplot(2,2,4);
plot(t_imp_ode,y_imp_ode,'LineWidth',2);
xlabel('t [s]');
ylabel('F_x(t) [N]');
grid;
title('Numerical Impulse Response (n = 1)');


% task 4

tol = 10;  % tolerance for opt_modes evaluation
[y_imp_opt,t_imp_opt,sys_opt,n_opt,err,MPend_opt] = opt_modes(MPend,tol);

figure
plot(t_imp_opt,y_imp_opt,'LineWidth',2);
xlabel('t [s]');
ylabel('F_x(t) [N]');
grid
title('Impulse Response for Optimized Modes Model');

clear MPend;

MPend_opt.n = n_opt;


% task 5
MPend_opt.x0 = zeros(2*MPend_opt.n,1);
MPend_opt.u_type = 'step';   % input
MPend_opt.step_t0 = 0;       % initial value
MPend_opt.tf = 100;          % final value
MPend_opt.step_amp = 1;      % step amplitude
 

% steady_state
s = tf('s');
Vf = dcgain(MPend_opt.sys);  % steady state response value

% built-in step response

[y_step_built_in,t_step] = slosh_response(MPend_opt);
Vf = Vf*ones(1,length(t_step));

figure;
subplot(1,3,1);
plot(t_step,y_step_built_in);
hold on
plot(t_step,Vf,LineWidth=2,LineStyle="--");
xlabel('t [s]');
ylabel('F_x(t) [N]');
grid;
title('Built-in Step Response');
legend('','steady-state')

% analytical step response

A = MPend_opt.sys.A;
B = MPend_opt.sys.B;
C = MPend_opt.sys.C;
D = MPend_opt.sys.D;
I = eye(2*MPend_opt.n);

inv_A = A \ I;
N = length(t_step);
y_step_analytical = zeros(N,1);

for i = 1:N
    y_step_analytical(i) = C*inv_A*(expm(A.*t_step(i)) - I)*B + D;
end

subplot(1,3,2);
plot(t_step,y_step_analytical);
hold on
plot(t_step,Vf,LineWidth=2,LineStyle="--");
xlabel('t [s]');
ylabel('F_x(t) [N]');
grid;
title('Analytical Step Response');
legend('','steady-state')

% numerical integration technique

[t_step_ode,y_step_ode] = MPend_Ode_Calc(MPend_opt);

subplot(1,3,3);
plot(t_step_ode,y_step_ode);
hold on
plot(t_step,Vf,LineWidth=2,LineStyle="--");
xlabel('t [s]');
ylabel('F_x(t) [N]');
grid;
title('Numerical Step Response');
legend('','steady-state')



% task 6

H = tf(MPend_opt.sys); % transfer function;
H_disp = H*s^2;       % transfer funtion with u as displacement

m_liq = MPend_opt.m_liq;  % liquid mass
H_frozen = m_liq*s^2;     % frozen liquid transfer function

% FRF evaluation
w = linspace(1,100,1e3);
H_disp_mag = freqresp(H_disp,w);
H_disp_mag = 20*log10(abs(H_disp_mag(:)));
H_frozen_mag = freqresp(H_frozen,w);
H_frozen_mag = 20*log10(abs(H_frozen_mag(:)));

figure
semilogx(w, H_disp_mag, 'b', 'LineWidth', 2);
hold on;
semilogx(w, H_frozen_mag, 'r--', 'LineWidth', 2);
xlabel('Frequency [rad/s]');
ylabel('Amplitude [dB]');
legend('Dynamic model', 'Frozen liquid model');
grid on;
title('Frequency response comparison');

% task 7

sim = sim('Model_26_01_1.slx');

% extracting array outputs
t_step_simulink = sim.tout;    
y_step_simulink = sim.y;   

figure
subplot(1,2,1)
plot(t_step_simulink,y_step_simulink);
title('Simulink Step Response');
xlabel('t [s]');
ylabel('F_x(t) [N]');
grid on
subplot(1,2,2)
plot(t_step,y_step_built_in);
title('Built-in Step Response');
xlabel('t [s]');
ylabel('F_x(t) [N]');
grid on


%% FUNCTIONS 
function [MPend] = parameters(MPend)
% function [MPend] = parameters(MPend)
%
% This function returns the parameters of an equivalent mechanical model of 
% the sloshing phenomenon in a cylindrical tank represented through pendulums.
%
% INPUT
% Name            Type                Size [when applicable]
% MPend           struct              1x1    - Input structure with tank and liquid parameters
%
% OUTPUT
% Name            Type                Size [when applicable]
% MPend           struct              1x1    - Updated structure with computed dynamic parameters
%
% Extract input parameters
d = MPend.d;             % Tank diameter [m]
r = d/2;
h = MPend.h;             % Liquid height [m]
rho = MPend.rho;         % Density of liquid [kg/m^3]
g = MPend.g;             % Gravitational acceleration [m/s^2]
n_opt = MPend.n;         % Number of modes to compute
m_liq = rho * pi * d^2 * h / 4;  % Total mass of the liquid [kg]



% Initialization of vectors for zeros of derivative of J1, natural frequencies, 
% associated masses, pendulum lengths, and hinge positions.
lambda = zeros(n_opt,1);   % Zeros of the derivative of Bessel function J1 (lambda_n)
omega_n = zeros(n_opt,1);  % Natural frequencies (rad/s) for each mode
m_n = zeros(n_opt,1);      % Masses associated with each mode [kg]
L_n = zeros(n_opt,1);      % Equivalent pendulum lengths [m]
H_n = zeros(n_opt,1);      % Hinge positions [m]

% Calculation of the zeros (lambda_n) using the derivative of Bessel function
% J'_1(x) = 0.5*(J0(x) - J2(x))
d1_fJ = @(x) 0.5*(besselj(0,x) - besselj(2,x));

guess = 1.8;  % Initial guess for the first mode
for n = 1:n_opt
    lambda(n) = fzero(d1_fJ, guess);
    guess = lambda(n) + pi;  % Successive guess increases roughly by pi
end

% Calculation of dynamic parameters for each mode
for n = 1:n_opt
    m_n(n) = m_liq * (d * tanh(2 * lambda(n) * h / d) / (lambda(n) * (lambda(n)^2 - 1) * h));
    L_n(n) = d / (2 * lambda(n) * tanh(2 * lambda(n) * h / d));
    H_n(n) = h / 2 - d / (2 * lambda(n)) * (tanh(lambda(n) * h / d) - 1 / sinh(2 * lambda(n) * h / d));
    omega_n(n) = sqrt(g * lambda(n) / r * tanh(lambda(n) * h / r));
end

% Calculation of the mass that is rigidly anchored and its corresponding hinge position
m_0 = m_liq - sum(m_n);  % Rigidly anchored mass [kg]
H_0 = sum(m_n .* (H_n - L_n) / m_0);  % Hinge position of the anchored mass [m]

% Update structure MPend with the computed parameters
MPend.omega_n = omega_n;
MPend.L_n = L_n;
MPend.m_n = m_n;
MPend.H_n = H_n;
MPend.m_0 = m_0;
MPend.H_0 = H_0;
MPend.x0 = [0; 0];                    % Initial conditions [displacement; velocity]
MPend.f_n = 2*pi * omega_n;           % Natural frequencies in Hz
MPend.m_liq = m_liq;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y,t,sys] = slosh_response(MPend)
% function [y,t,sys] = slosh_response(MPend)
%
% slosh_response
% This function computes the response of the sloshing system using a 
% state-space model which represents the equivalent pendulum system.
%
% INPUT
% Name            Type                Size [when applicable]
% MPend           struct              1x1    - Structure containing the system parameters
%
% OUTPUT
% Name            Type                Size [when applicable]
% y               vector              nx1    - System output (force response)
% t               vector              nx1    - Time vector
% sys             ss                  - State-space object representing the system dynamics
%
n = MPend.n;
omega_n = MPend.omega_n(1:n);
L_n = MPend.L_n(1:n);
xi = MPend.xi;
m_n = MPend.m_n(1:n);
m_0 = MPend.m_0;
t_0 = MPend.t0;
t_final = MPend.tf;

% Initialize state-space matrices
A = zeros(n, n);
B = zeros(n, 1);
C = zeros(1, n);
D = -m_0;

% Loop to construct block-diagonal system matrices for each mode
for i = 1 : n
    A_tmp = [0, 1; -omega_n(i)^2, -2 * xi * omega_n(i)];
    B_tmp = [0; -1 / L_n(i)];
    C_tmp = [m_n(i) * L_n(i) * omega_n(i)^2, 2 * m_n(i) * L_n(i) * xi * omega_n(i)];
    
    if i == 1
        A = A_tmp;
        B = B_tmp;
        C = C_tmp;
    else
        A = blkdiag(A, A_tmp);
        B = [B; B_tmp];
        C = [C, C_tmp];
    end
end

% Create state-space system object
sys = ss(A, B, C, D);

% Compute the system response based on the specified input type
switch MPend.u_type
    case {'impulse'}
        [y, t] = impulse(sys, [t_0, t_final]);
    case {'step'}
        [y, t] = step(sys, [t_0, t_final]);
    otherwise
        error('Input type not defined/found');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function force = MPend_input(t, MPend)
% function force = MPend_input(t, MPend)
%
% MPend_input
% This function defines the input force acting on the sloshing system based on 
% the specified input type (impulse or step).
%
% INPUT
% Name            Type                Size [when applicable]
% t               scalar              1x1    - Time instant at which the force is computed
% MPend           struct              1x1    - Structure containing input parameters and type
%
% OUTPUT
% Name            Type                Size [when applicable]
% force           scalar              1x1    - Computed input force at time t
%
switch MPend.u_type
    case {'impulse'}
        if t <= 1e-5
            force = 1e5;
        else
            force = 0;
        end
    case {'step'}
        if t >= MPend.step_t0
            force = MPend.step_amp;
        else
            force = 0;
        end    
    otherwise
        error('Input type not defined/found');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t, y] = MPend_Ode_Calc(MPend)
% function [t, y] = MPend_Ode_Calc(MPend)
%
% MPend_Ode_Calc
% This function computes the system response by solving the system of ODEs 
% associated with the sloshing model.
%
% INPUT
% Name            Type                Size [when applicable]
% MPend           struct              1x1    - Structure containing system parameters and initial conditions
%
% OUTPUT
% Name            Type                Size [when applicable]
% t               vector              nx1    - Time vector of the integration interval
% y               vector              nx1    - System output computed using the state-space equations
%
    ODE_obj = ode;  % Create ODE object (assumed custom class for ODE solving)
    ODE_obj.ODEFcn = @(t, x) MPend_f(t, x, @MPend_input, MPend);
    ODE_obj.InitialValue = MPend.x0;
    ODE_obj.Solver = 'ode45';
    ODEResults_obj = solve(ODE_obj, MPend.t0, MPend.tf);

   
    t = ODEResults_obj.Time;
    x = ODEResults_obj.Solution;
    
    % Extract the output using the state-space matrices stored in MPend.sys
    C = MPend.sys.C;
    D = MPend.sys.D;
    
    % Compute the output y as a function of the state x and the input force
    y = C * x + D * MPend_input(t, MPend);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xdot = MPend_f(t, x, MPend_input, MPend)
% function xdot = MPend_f(t, x, MPend_input, MPend)
%
% MPend_f
% This function defines the set of ordinary differential equations (ODEs)
% for the sloshing system modeled with equivalent pendulums.
%
% INPUT
% Name            Type                Size [when applicable]
% t               scalar              1x1    - Current time instant
% x               vector              2n x 1 - State vector [displacements; velocities]
% MPend_input     function_handle     - Handle to the function computing the input force
% MPend           struct              1x1    - Structure containing the system parameters
%
% OUTPUT
% Name            Type                Size [when applicable]
% xdot            vector              2n x 1 - Time derivative of the state vector
%
    n = MPend.n;
    u = MPend_input(t, MPend);
    xdot = zeros(2 * n, 1);

    for i = 1:n
        xdot(i) = x(n + i);  % Derivative of displacement is velocity
        xdot(n + i) = -MPend.omega_n(i)^2 * x(i) - 2 * MPend.omega_n(i) * MPend.xi * x(n + i) - u / MPend.L_n(i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y, t, sys, n_opt, err, MPend_opt] = opt_modes(MPend, tol)
% function [y, t, sys, n_opt, err, MPend_opt] = opt_modes(MPend, tol)
%
% opt_modes
% This function performs an optimization by increasing the number of modes 
% until the maximum error between the response of the system of order n and
% that of order n+1 is below a specified tolerance.
%
% INPUT
% Name            Type                Size [when applicable]
% MPend           struct              1x1    - Structure containing the initial system parameters
% tol             scalar              1x1    - Tolerance for the maximum allowable error between successive models
%
% OUTPUT
% Name            Type                Size [when applicable]
% y               vector              nx1    - Response output corresponding to the optimized number of modes
% t               vector              nx1    - Time vector for the simulation
% sys             ss                  - State-space model of the optimized system
% n_opt           scalar              1x1    - Optimized number of modes
% err             scalar              1x1    - Final maximum error between successive model responses
% MPend_opt       struct              1x1    - Updated structure with optimized parameters
%
err = 1e3;
y_mat = [];
MPend.n = 0;
MPend_opt = MPend;

while err > tol
    MPend_opt.n = MPend_opt.n + 1;
    MPend_opt = parameters(MPend_opt);
    [y, t, sys] = slosh_response(MPend_opt);
    
    if MPend_opt.n > 1
        err = max(abs(y - y_mat(:, end)));
    end
    
    y_mat = [y_mat, y];  % Store successive outputs for error evaluation
    MPend_opt.sys = sys;
end

n_opt = MPend_opt.n;

% Print summary of optimized model configuration
fprintf('The optimized modes number for type ''%s'' input and %f accuracy is %d.\n',...
    MPend.u_type, tol, n_opt);

end
