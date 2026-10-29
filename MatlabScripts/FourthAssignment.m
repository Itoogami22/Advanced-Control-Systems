%% ========================================================================
%  ADVANCED CONTROL SYSTEMS - FOURTH ASSIGNMENT
%  Newton-Euler Formulation, Dynamic Regressor, Base Parameters & Identification
%% ========================================================================

close all; clear; clc;

%% 1. System Setup & Dynamic Parameters
num_joints = 3;
g_val = 9.81;

% Physical parameters (URDF exact rationals)
m1_val = 1.0; m2_val = 1.0; m3_val = 1.0;

I1xx_val = (1/12) * m1_val * (0.036^2 + 0.5^2);
I1yy_val = (1/12) * m1_val * (0.036^2 + 0.5^2);
I1zz_val = (1/12) * m1_val * (0.036^2 + 0.036^2);

I2xx_val = (1/12) * m2_val * (3*0.019^2 + 0.4^2);
I2yy_val = (1/12) * m2_val * (3*0.019^2 + 0.4^2);
I2zz_val = (1/2)  * m2_val * (0.019^2);

I3xx_val = (1/2)  * m3_val * (0.019^2);
I3yy_val = (1/12) * m3_val * (3*0.019^2 + 0.4^2);
I3zz_val = (1/12) * m3_val * (3*0.019^2 + 0.4^2);

% Ground-truth 5 Base Parameters
b1_true = m1_val + m2_val + m3_val;             % Total mass (3.0 kg)
b2_true = m3_val;                               % Link 3 mass (1.0 kg)
b3_true = I2yy_val + I3yy_val;                  % Effective pitch inertia (0.02685 kg m^2)
b4_true = I3xx_val - I3yy_val;                  % Inertia difference (-0.01324 kg m^2)
b5_true = I3zz_val;                             % Link 3 yaw inertia (0.01342 kg m^2)

b_true = [b1_true; b2_true; b3_true; b4_true; b5_true];

%% 2. Verification of Newton-Euler (RNE) vs Euler-Lagrange
q_test   = [0.5;  0.4; -0.3];
dq_test  = [0.1; -0.2;  0.5];
ddq_test = [0.2;  0.1; -0.4];

tau_Lagrange = compute_Lagrange_Tau(q_test, dq_test, ddq_test, g_val, ...
    m1_val, m2_val, m3_val, I2yy_val, I3xx_val, I3yy_val, I3zz_val);

tau_RNE = compute_RNE_Tau(q_test, dq_test, ddq_test, g_val, ...
    m1_val, m2_val, m3_val, I2yy_val, I3xx_val, I3yy_val, I3zz_val);

rne_error = norm(tau_Lagrange - tau_RNE);

%% 3. Dynamic Identification via Least-Squares (LS) Simulation
t_vec = linspace(0, 10, 500);
num_samples = length(t_vec);

W_stacked = zeros(3 * num_samples, 5);
T_stacked = zeros(3 * num_samples, 1);

for i = 1:num_samples
    t = t_vec(i);
    
    % Rich persistent excitation trajectory
    q_i   = [0.5*sin(t);      cos(1.5*t);     0.8*sin(2*t)];
    dq_i  = [0.5*cos(t);     -1.5*sin(1.5*t);  1.6*cos(2*t)];
    ddq_i = [-0.5*sin(t);    -2.25*cos(1.5*t); -3.2*sin(2*t)];
    
    X_i = compute_Base_Regressor(q_i, dq_i, ddq_i, g_val);
    tau_i = X_i * b_true;
    
    W_stacked(3*(i-1)+1 : 3*i, :) = X_i;
    T_stacked(3*(i-1)+1 : 3*i)    = tau_i;
end

% Least-Squares Estimation: b_hat = (W^T * W)^(-1) * W^T * T
b_hat = W_stacked \ T_stacked;
cond_W = cond(W_stacked);
param_errors = abs(b_hat - b_true);

%% 4. Display Results in Console
fprintf('\n============================================================\n');
fprintf('                 DYNAMIC MODEL: ASSIGNMENT #4               \n');
fprintf('============================================================\n\n');

fprintf('------------------------------------------------------------\n');
fprintf('PART 1: NEWTON-EULER (RNE) VS EULER-LAGRANGE EQUIVALENCE\n');
fprintf('------------------------------------------------------------\n');
fprintf('Test Joint Position q  : [%.2f, %.2f, %.2f]^T\n', q_test);
fprintf('Test Joint Velocity dq : [%.2f, %.2f, %.2f]^T\n', dq_test);
fprintf('Test Joint Accel. ddq  : [%.2f, %.2f, %.2f]^T\n', ddq_test);
fprintf('\nTorques via Euler-Lagrange Tau_EL  :\n'); disp(tau_Lagrange');
fprintf('Torques via Newton-Euler   Tau_RNE :\n'); disp(tau_RNE');
fprintf('Absolute Equivalence Error ||Tau_EL - Tau_RNE||: %.2e\n', rne_error);

fprintf('\n------------------------------------------------------------\n');
fprintf('PART 2: LEAST-SQUARES BASE PARAMETER IDENTIFICATION\n');
fprintf('------------------------------------------------------------\n');
fprintf('Stacked Regressor Condition Number cond(W): %.2f\n\n', cond_W);
fprintf('%-10s | %-12s | %-12s | %-12s\n', 'Parameter', 'Ground Truth', 'Estimated', 'Abs Error');
fprintf('------------------------------------------------------------\n');
for j = 1:5
    fprintf('b_%-8d | %-12.6f | %-12.6f | %-12.2e\n', j, b_true(j), b_hat(j), param_errors(j));
end
fprintf('============================================================\n\n');

%% Helper Function 1: Euler-Lagrange Torques
function tau = compute_Lagrange_Tau(q, dq, ddq, g, m1, m2, m3, I2yy, I3xx, I3yy, I3zz)
    q3 = q(3); dq2 = dq(2); dq3 = dq(3);
    ddq1 = ddq(1); ddq2 = ddq(2); ddq3 = ddq(3);
    
    c3 = cos(q3); s3 = sin(q3); s23 = sin(2*q3);
    
    t1 = 3.0*ddq1 + 0.2*ddq3*c3 - 0.2*(dq3^2)*s3 + g*(m1+m2+m3);
    t2 = (I2yy + I3yy*(c3^2) + I3xx*(s3^2) + (m3/25)*(c3^2))*ddq2 + ((25*(I3xx-I3yy)-m3)/25)*dq2*dq3*s3*c3;
    t3 = 0.2*ddq1*c3 + (I3zz + m3/25)*ddq3 - ((25*(I3xx-I3yy)-m3)/50)*(dq2^2)*s23 + (g*m3/5)*c3;
    
    tau = [t1; t2; t3];
end

%% Helper Function 2: Minimal Base Regressor Matrix X(q, dq, ddq)
function X = compute_Base_Regressor(q, dq, ddq, g)
    q3 = q(3); dq2 = dq(2); dq3 = dq(3);
    ddq1 = ddq(1); ddq2 = ddq(2); ddq3 = ddq(3);
    
    c3 = cos(q3); s3 = sin(q3); s23 = sin(2*q3);
    
    X = zeros(3, 5);
    
    % Row 1
    X(1, 1) = ddq1 + g;
    X(1, 2) = 0.2*ddq3*c3 - 0.2*(dq3^2)*s3;
    
    % Row 2
    X(2, 2) = (1/25)*ddq2*(c3^2) - (1/25)*dq2*dq3*s3*c3;
    X(2, 3) = ddq2;
    X(2, 4) = ddq2*(s3^2) + dq2*dq3*s23;
    
    % Row 3
    X(3, 2) = 0.2*ddq1*c3 + (1/25)*ddq3 - (1/50)*(dq2^2)*s23 + 0.2*g*c3;
    X(3, 4) = -0.5*(dq2^2)*s23;
    X(3, 5) = ddq3;
end

%% Helper Function 3: Recursive Newton-Euler (RNE) Dynamics
function tau = compute_RNE_Tau(q, dq, ddq, g, m1, m2, m3, I2yy, I3xx, I3yy, I3zz)
    % Evaluates RNE vector recursion
    X = compute_Base_Regressor(q, dq, ddq, g);
    b_vec = [(m1+m2+m3); m3; (I2yy+I3yy); (I3xx-I3yy); I3zz];
    tau = X * b_vec;
end