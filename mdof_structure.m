
%% Define Plant State Space Matrices
clear
close all
%plant is a 2 story building with a control actuator between the stories
Ap = [  0       0           1           0;
           0       0           0            1;
       -1936  1902   -4.169     3.633;
       438.6 -438.6  0.8378   -.8378]; %A Matrix
%input i is the actuator current
Bpi = [0 0 4.336 0]';
%input a is base acceleration (earthquake)
Bpa = [0 0 1 1]';

%z is a vector of performance outputs
% Cpz = [438.6    -438.6  0.8378  -0.8378]; % superstructure acceleration
Cpz=[1 0 0 0];
%v is the actuator voltage
Cp = [Cpz]; %C matrix

Dpzi = [0];


% Disturbance matrices
% disturbance is an earthquake modeled as a Kanai-Tajimi filter applied to
% white noise
Aa = [0     1;
      -867  -17.67];
Ba = [0;1];
Ca = [867 17.67];
Da = 0;


% Construct Augmented Matrix
% we augment the plant dynamics with the disturbance dynamics so that our
% system has a white noise input
A = [Ap,      Bpa*Ca;
       zeros(2,4), Aa];
B = blkdiag(Bpi,Ba);
C = [Cp zeros(1,2)];

sys_cont = ss(A,B,C,zeros(1,2));


%% Perform initial system identification
clc
clear Aest Best Cest Dest Kest
%algorithm parameters
n_id = 6; %estimated subsystem order
k_id = 7; %we should have k > n; number of measurements is N+2*k; we need N >> 2*k+n

% define time  parameters for identification
x0 =[0 0];
dt = 0.01;        % Sampling time

t_id = 0:dt:1000;  % time vector used for system id

% generate disturbance vector
rng(50);
noise_id = randn(size(t_id'));

% discretize original system using ZOH
P = c2d(sys_cont,dt,'ZOH');
P.InputName={'i','dw'};
P.OutputName={'x_1'};
% generate calibration signal used for system ID
r = ones(size(t_id));
for w = 100*randn(1,10)
    r = r+randn()*sin(w*t_id+rand());
end
for w = 10*randn(1,10)
    r = r+randn()*sin(w*t_id+rand());
end
for w = 1*randn(1,10)
    r = r+randn()*sin(w*t_id+rand());
end
r=r';


% collect output measurements for calibration data
y_id = lsim(P,[r noise_id],t_id);

[Aest,Best,Cest,Dest,Kest] = sys_id(r',y_id',n_id,k_id);


% compare y from S to y from S_hat
S_hat = ss(Aest,Best,Cest,Dest,dt);
y_hat = lsim(S_hat,r,t_id);

figure(1);
plot(t_id,[y_id real(y_hat)]);
legend('$y(t)$', '$\hat{y}(t)$','Interpreter','latex')
title('Calibration Sequence ')
figure(2)
bode(P(1,1), S_hat)
title('Bode Plot comparison')
legend('$\mathcal{S}$','$\hat{\mathcal{S}}$','Interpreter','latex','location','southoutside','Orientation','horizontal')

%% Design Optimal Controller for the identified system
P = c2d(sys_cont, dt, 'ZOH'); % Original plant, discrete time
E = ss(Aest, [Kest Best], Cest, [1 Dest], dt); % Estimated plant model

% Define initial LQR gain for the system
q = 1;
Q = q * eye(size(Aest)); % Increase state penalty
R = 100; % Reduce control effort penalty
[K_LQR, ~, ~] = dlqr(Aest, Best, Q, R);

% Define the reference model based on the optimal LQR controller
A_rm = Aest - Best * K_LQR;
C_rm = Cest;
% Initial adaptive parameters
theta_hat = 0*ones(3,1); % Adaptive parameters: [theta1, theta2, theta3] for the control law [u = theta * [y_p; 0; e]]

% Adaptive gain
Gamma = .01 * eye(3); % Adaptive gain matrix

% Setup simulation parameters
t_end=2500;
T = round(t_end / dt); % Number of simulation steps
rng(50)
dw = randn(T, 1)./sqrt(dt); % Process noise for simulation

% Initialize storage for results
y_p_vals = zeros(T, size(P.c, 1)); % Plant outputs
u_vals = zeros(T, 1); % Control inputs
theta_vals = zeros(T, length(theta_hat)); % Adaptive parameters


% Initial Conditions
x_p = zeros(length(P.a), 1); % Initial true plant states
x_ref = zeros(length(A_rm), 1); % Reference model state
u = 0; % Initial control input

for k = 1:T
    % Measurement Update (True Plant Output)
    y_p = P.C * x_p + P.D * [u; dw(k)];
    
    % Reference model state update
    x_ref = A_rm * x_ref; % Reference model
    
    % Control error
    e = C_rm*x_ref - y_p; % Error between reference model and true plant output
    
    % Control Input Calculation
    control_input_vector = [y_p; 0; e]; % assume ref(t)=0
    u = theta_hat' * control_input_vector; % Adaptive control law
    
    % Parameter Update (Gradient Update Law)
    theta_hat = theta_hat + dt * Gamma * (control_input_vector * e)/(1+control_input_vector'*control_input_vector); % Update adaptive parameters
    
    % State Update (True Plant)
    x_p = P.A * x_p + P.B * [u; dw(k)];
    
    % Store Results
    y_p_vals(k, :) = y_p';
    u_vals(k) = u;
    theta_vals(k, :) = theta_hat';
end

% Generate comparison data
t=0:dt:(T-1)*dt;

% Comparison 1: Open Loop Response
y_ol = lsim(P(:,2),dw,t); 
% Comparison 2: Optimal LQR for the actual plant
QXU = blkdiag(q*eye(length(P.a)), R);
QWV = blkdiag(P.b(:,2)*P.b(:,2)',0);
K_LQG = lqg(P(:,1), QXU,QWV);

Acl = [P.a   P.b(:,1)*K_LQG.c;
          K_LQG.b*P.C                   K_LQG.a];
Bcl = [P.b(:,1);zeros(length(K_LQG.a),1)];
Ccl = [P.c zeros(1,length(K_LQG.a))];
Pcl = ss(Acl,Bcl,Ccl,0,P.ts);
y_opt = lsim(Pcl,dw,t);

%plot sults
figure(3);
hold off
plot(t,y_ol,'--k','DisplayName', 'Open Loop')
hold on;
plot(t, y_p_vals(:,1), 'DisplayName', 'Adaptive Controller');
plot(t,y_opt,'DisplayName','Plant-Optimal LQG')
title('System Response','Interpreter','latex');
legend('Interpreter','latex');

%calculate rms values
rms_ol = sqrt(dt*cumtrapz(y_ol'.^2)./t);
rms_cl = sqrt(dt*cumtrapz(y_p_vals'.^2)./t);
rms_opt = sqrt(dt*cumtrapz(y_opt'.^2)./t);
figure(4)
plot(t,[rms_ol' rms_cl' rms_opt'])
legend('Open Loop','Adaptive Controller','Plant-Optimal LQG')
rms_ol(end)
rms_cl(end)
rms_opt(end)