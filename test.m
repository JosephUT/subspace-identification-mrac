clear; close all

rng(5);

A = [0   1;
     -2 -1];

B = [1;0];

K = [1;1];

C = [1 0];

dt = .001;
 t = 0:dt:1;
dw = randn(size(t))';

sys = ss(A,[B K],C,[0 0]);
sys_d = c2d(sys,dt);


r = zeros(size(t));
for w = 100*rand(1,10)
    r = r+rand()*sin(w*t+rand());
end

r=r';

y = lsim(sys_d,[r, dw],t,[0 0]');
% r = (r-mean(r))./std(r);
% y = (y-mean(y))./std(y);

[Aest,Best,Cest,Dest, Kest] = sys_id(r',y',1,2);
sys_est = ss(Aest,[Best],Cest,[Dest],dt);
% sys_n4 = n4sid(r,y,2,'Ts',dt); %check against matlab sys id toolbox

%% Check ID'd system
yhat = lsim(sys_est,[r],t);
% y_chk = lsim(sys_n4,r,t);
% plot results
figure(1)

plot(t,[y yhat])
legend('Plant', 'ORT Estimate')

figure(2)
bode(sys_d, sys_est)

%% LQR design example
% Define LQR weighting matrices
Q = 2*eye(size(Aest)); % State weight
R = eye(size(Best, 2)); % Control weight

% Compute LQR gain
[K_lqr, ~, ~] = dlqr(Aest,Best, Q, R);

% Design Kalman filter
[L, P, E] = kalman(ss(Aest, [Best Kest], Cest, [Dest zeros(size(Cest, 1), size(Kest, 2))]),1,1,0);
t = 0:dt:100;
dw = randn(size(t))';
r = ones(length(t),1);

sys_cl = ss(sys_d.A-sys_d.B(:,1)*K_lqr,sys_d.B,sys_d.C,sys_d.D,dt);


y2 = lsim(sys_cl,[r dw],t);

figure(3)
plot(t,[y2 r])





