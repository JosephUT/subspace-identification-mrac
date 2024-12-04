%% State Space Model WAAM Temperature Control

clear; close all


Cs = 526.3;             % substrate heat capacity
kc = 7.25;             % Thermal conductivity
As = 2.5e-4;             % contact area 
Kc = 10;              % heat transfer coefficient for bulk
Ac = 0.000250;        % contact area for bulk cooling
Beta = 1e-6;          % cooling rate for bulk
eta = .84;            % arc efficiency
Vdc = 15;         % DC component of voltage signal

% System in original coordinates
A = [-kc*As/Cs      Kc*As/Cs;
         Kc*Ac/Cs      -Beta];

Bi = [eta*Vdc/Cs;0];

K = [sqrt(2)*Bi(1)/Vdc;0];

C = [1 0];
sys_cont = ss(A,[Bi K],C,[0 sqrt(1e-5)]);
sys_cont_modal = modalreal(sys_cont,'Normalize',true);
%% Obtain Initial Estimate

%algorithm parameters
n = 1; %estimated subsystem order
k = 2; %we should have k > n where n is the system order N >> 2*k+n

% define time  parameters for identification
x0 =[300 300]; %we need zero ic's for the algorithm, so we will shift all temps by -300;
dt = 0.001;        % Sampling time
t_id = 0:dt:100;  % time vector used for system id

% generate disturbance vector
rng(50);
noise_id = randn(size(t_id'));

% discretize original system using ZOH
sys_disc = c2d(sys_cont,dt,'tustin');
sys_disc.InputName={'u','dV'};
sys_disc.OutputName={'T_s'};
% generate calibration signal used for system ID
r = 1575*ones(size(t_id));
% for w = 10*randn(1,500)
%     r = r+randn()*sin(w*t_id+rand());
% end
r=r';


% collect output measurements for calibration data
y_id = lsim(sys_disc,[r noise_id],t_id);

[Aest,Best,Cest,Dest,Kest] = sys_id(r',y_id',n,k);

% compare y from S to y from S_hat
S_hat = ss(Aest,Best,Cest,Dest,dt);
y_hat = lsim(S_hat,r,t_id);

figure(1);
plot(t_id,[y_id y_hat]);
legend('$y(t)$', '$\hat{y}(t)$','Interpreter','latex')
title('Calibration Sequence ')
figure(2)
bode(sys_disc(1,1), S_hat)
title('Bode Plot comparison')
legend('$\mathcal{S}$','$\hat{\mathcal{S}}$','Interpreter','latex','location','southoutside','Orientation','horizontal')

%% Design initial PI gain for the system
% Define initial LQR gain for the system
rng(50)
% Extract state-space matrices of the original system
[A,B,C,D] = ssdata(sys_disc);
K = B(:,2); % Process noise influence matrix
B = B(:,1); % Control input matrix

%Only x1, the deterministic subspace, is controllable from u. We want to
%design a control input u such that x1->ref, which we can do with PI
%control. First, we will augment our estimated plant with an integrator
%state, where dot{x3}= r-y = r-Cx

A_aug = [Aest zeros(2,1);
              -Cest 0];
B_aug = [Best zeros(2,1);-Dest 1]; %input 1 is u, input 2 is ref

C_aug = [Cest 0];
D_aug = [Dest 0];

%next we design a PI controller of the form u = -Lx_hat+I(r-y_hat)
L = [.9313 0];
I = 10;
%setup simulation
t_end = 1000; %sim length
T = t_end/sys_disc.Ts; % Number of simulation time steps

%reference signal (ramp up to constant)
ramp_time = floor(T/30);
Tref = 1575; % (offset by 300 degrees)
ref = Tref*ones(T,1); % Reference signal 
ref(1:ramp_time)=linspace(0,Tref,ramp_time);
% Results vectors
y_vals = zeros(T, 1); %plant outputs
y_hat_vals = zeros(T,1);   %estimated plant outputs
u_vals = zeros(T, 1); %plant inputs
x_hat_vals = zeros(T, 3); %states of estimated plant
x_vals = zeros(T,2); %states of plant
n_vals = zeros(T,1); %error signal

%initial conditions
x = [0;0]; % plant states (offset by 300 degrees)
x_hat = [0;0;0]; %
y=0;
y_hat=0;
for k = 1:T
    %controller update
    

    u = L*(x_hat(1:2))+I*x_hat(3);


    % Measurement updates
    y = C*x + D*[0;randn()]; 
    y_hat = C_aug*x_hat+D_aug*[u;ref(k)];
    
    
    %state updates
    x = A*x+B*u+K*randn();

    %estimate updates
    x_hat = A_aug*x_hat+B_aug*[u;ref(k)];
    
    
    % Store results
    y_vals(k, :) = y';
    u_vals(k, :) = u';
    x_vals(k, :) = x';
    y_hat_vals(k,:) = y_hat';
end

  

% Plot results
t = ((1:T) - 1) * sys_disc.Ts;
figure(3);
hold off
plot(t, 300+y_vals, 'DisplayName', '$y(t)$');
hold on;
plot(t, 300+y_hat_vals,'DisplayName', '$\hat{y}(t)$')
plot(t, 300+ref, '--', 'DisplayName', '$r(t)$');
title('Output y(t) and Reference r');
legend('Interpreter','latex');



