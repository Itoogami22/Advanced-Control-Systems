%% ========================================================================
%  ADVANCED CONTROL SYSTEMS - SECOND ASSIGNMENT
% ========================================================================

close all; clear; clc;

%% Robot Import
robot = importrobot('00001_PRR_zzx FILIPPO SEGALLA.urdf');
robot.DataFormat = 'column'; % Set column format for direct vector indexing
showdetails(robot)

%% Symbolic Variables Definition
num_joints = 3;           

% Symbolic variables for DH parameters
syms a alpha d theta real

% Generalized coordinates and velocities
syms q1 q2 q3 real                                  
syms dq1 dq2 dq3 real                               

% Symbolic dynamic parameters (letters)
syms m1 m2 m3 real                                  % Link masses (kg)
syms g real                                         % Gravitational acceleration (m/s^2)
syms I1xx I1yy I1zz real                            % Link 1 inertia moments
syms I2xx I2yy I2zz real                            % Link 2 inertia moments
syms I3xx I3yy I3zz real                            % Link 3 inertia moments

q  = [q1; q2; q3];
dq = [dq1; dq2; dq3];
m  = [m1, m2, m3];

% Gravity acceleration vector expressed in Base Frame Sigma_0
g0 = [0; 0; -g];                                    

% Local inertia tensors relative to CoMs in local link frames
I_local = cell(num_joints, 1);
I_local{1} = diag([I1xx, I1yy, I1zz]);
I_local{2} = diag([I2xx, I2yy, I2zz]);
I_local{3} = diag([I3xx, I3yy, I3zz]);

%% Direct Kinematics and Transformation Matrices (DH)
A = [cos(theta), -sin(theta)*cos(alpha),  sin(theta)*sin(alpha), a*cos(theta); ...
     sin(theta),  cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta); ...
     0,           sin(alpha),             cos(alpha),            d; ...
     0,           0,                      0,                     1];

TB0 = [ 0, -1,  0,  0;
        1,  0,  0,  0;
        0,  0,  1,  0.65;
        0,  0,  0,  1 ];

DH_table = [ 0,    0,     q1,   0;
             0,    pi/2,  0.4,  q2;
             0.4,  0,     0,    q3 ];

joint_types = ['P'; 'R'; 'R'];                     

T = cell(num_joints, 1);
R = cell(num_joints, 1);
T_curr = sym(TB0);

for i = 1:num_joints
    A_i = subs(A, [a, alpha, d, theta], DH_table(i,:));
    T_curr = simplify(T_curr * A_i);
    T{i} = T_curr;                                
    R{i} = T{i}(1:3, 1:3);                         
end

%% Center of Mass (CoM) Positions in Base Frame Sigma_0
p_l_local = {[0; 0; -0.25; 1], [0; 0; -0.20; 1], [-0.20; 0; 0; 1]};
p_l = cell(num_joints, 1);

for i = 1:num_joints
    p_temp = simplify(T{i} * p_l_local{i});
    p_l{i} = p_temp(1:3);
end

%% Partial Jacobians for Link CoMs (JP_li, JO_li)
JP_l = cell(num_joints, 1);
JO_l = cell(num_joints, 1);

z_axes = cell(num_joints+1, 1);
p_origins = cell(num_joints+1, 1);

z_axes{1} = sym(TB0(1:3, 3));
p_origins{1} = sym(TB0(1:3, 4));

for i = 1:num_joints
    z_axes{i+1} = T{i}(1:3, 3);
    p_origins{i+1} = T{i}(1:3, 4);
end

for i = 1:num_joints
    JP_i = sym(zeros(3, num_joints));
    JO_i = sym(zeros(3, num_joints));
    
    for j = 1:i
        z_prev = z_axes{j};
        p_prev = p_origins{j};
        
        if joint_types(j) == 'P'
            JP_i(:, j) = z_prev;                  
            JO_i(:, j) = [0; 0; 0];               
        else % Revolute 'R'
            JP_i(:, j) = cross(z_prev, (p_l{i} - p_prev)); 
            JO_i(:, j) = z_prev;                  
        end
    end
    
    JP_l{i} = simplify(JP_i);
    JO_l{i} = simplify(JO_i);
end

%% Part 1: Generalized Inertia Matrix B(q) and Energy (Symbolic Letters)
B_symbolic = sym(zeros(num_joints, num_joints));

for i = 1:num_joints
    B_symbolic = B_symbolic + m(i) * (JP_l{i}' * JP_l{i}) + JO_l{i}' * R{i} * I_local{i} * R{i}' * JO_l{i};
end
B_symbolic = simplify(B_symbolic);

% Kinetic Energy T(q, dq) = 0.5 * dq' * B(q) * dq
T_kinetic_symbolic = simplify(0.5 * dq' * B_symbolic * dq);

% Potential Energy U(q) = - sum( m_i * g0' * p_li )
U_potential_symbolic = sym(0);
for i = 1:num_joints
    U_potential_symbolic = U_potential_symbolic - m(i) * (g0' * p_l{i});
end
U_potential_symbolic = simplify(U_potential_symbolic);

%% Part 2: Numerical Value Substitution (Exact Symbolic Rationals)
m1_val = sym(1);  m2_val = sym(1);  m3_val = sym(1);
g_val  = sym('981/100'); % 9.81 m/s^2 exact rational

I1xx_val = sym('1/12') * m1_val * (sym('36/1000')^2 + sym('1/2')^2);
I1yy_val = sym('1/12') * m1_val * (sym('36/1000')^2 + sym('1/2')^2);
I1zz_val = sym('1/12') * m1_val * (sym('36/1000')^2 + sym('36/1000')^2);

I2xx_val = sym('1/12') * m2_val * (3*sym('19/1000')^2 + sym('4/10')^2);
I2yy_val = sym('1/12') * m2_val * (3*sym('19/1000')^2 + sym('4/10')^2);
I2zz_val = sym('1/2')  * m2_val * (sym('19/1000')^2);

I3xx_val = sym('1/2')  * m3_val * (sym('19/1000')^2);
I3yy_val = sym('1/12') * m3_val * (3*sym('19/1000')^2 + sym('4/10')^2);
I3zz_val = sym('1/12') * m3_val * (3*sym('19/1000')^2 + sym('4/10')^2);

sym_vars = [m1, m2, m3, g, ...
            I1xx, I1yy, I1zz, ...
            I2xx, I2yy, I2zz, ...
            I3xx, I3yy, I3zz];

num_vals = [m1_val, m2_val, m3_val, g_val, ...
            I1xx_val, I1yy_val, I1zz_val, ...
            I2xx_val, I2yy_val, I2zz_val, ...
            I3xx_val, I3yy_val, I3zz_val];

B_num_exact = subs(B_symbolic, sym_vars, num_vals);
T_num_exact = subs(T_kinetic_symbolic, sym_vars, num_vals);
U_num_exact = subs(U_potential_symbolic, sym_vars, num_vals);

% Format clean decimal outputs with 4 significant digits (vpa + expand)
B_numeric = vpa(expand(B_num_exact), 4);
T_numeric = vpa(expand(T_num_exact), 4);
U_numeric = vpa(expand(U_num_exact), 4);

%% Part 3: Sample Numerical Evaluation & 3D Robot Plot at q_sample
% URDF Joint Limits Verification:
%   Joint 1 (Prismatic): q1 in [-0.5, 0.0] m  -> Selected q1_sample = -0.25 m
%   Joint 2 (Revolute) : q2 in [-pi, pi] rad  -> Selected q2_sample = pi/4 rad (~0.7854 rad)
%   Joint 3 (Revolute) : q3 in [-pi, pi] rad  -> Selected q3_sample = pi/6 rad (~0.5236 rad)
%   Joint Velocities   : dq_sample = [0.1; 0.2; -0.1] (m/s, rad/s)

q1_sample = -0.25;
q2_sample = pi / 4;
q3_sample = pi / 6;

dq1_sample =  0.1;
dq2_sample =  0.2;
dq3_sample = -0.1;

q_sample  = [q1_sample; q2_sample; q3_sample];
dq_sample = [dq1_sample; dq2_sample; dq3_sample];

% Evaluate numerical values at sample state
B_eval_sample = double(subs(B_num_exact, [q1, q2, q3], [q1_sample, q2_sample, q3_sample]));
T_eval_sample = double(subs(T_num_exact, [q1, q2, q3, dq1, dq2, dq3], ...
    [q1_sample, q2_sample, q3_sample, dq1_sample, dq2_sample, dq3_sample]));
U_eval_sample = double(subs(U_num_exact, [q1, q2, q3], [q1_sample, q2_sample, q3_sample]));
E_eval_sample = T_eval_sample + U_eval_sample;

% Interactive 3D Visualization of the Robot at q_sample
fig = figure('Name', 'PRR Robot 3D Plot given configuration', 'NumberTitle', 'off');
show(robot, q_sample, 'PreservePlot', false, 'Frames', 'on');

title(sprintf('PRR Manipulator at q_{sample} = [%.2f m, %.4f rad, %.4f rad]^T', ...
    q1_sample, q2_sample, q3_sample), 'FontSize', 12, 'FontWeight', 'bold');

rotate3d on;              
cameratoolbar('Show');    
grid on;                  
axis equal;               

xlim([-1.5, 1.5]);        
ylim([-1.5, 1.5]);        
zlim([-0.2, 2.0]);        

%% Display Results
fprintf('------------------------------------------------------------\n');
fprintf('PART 1: SYMBOLIC FORMULAS (LETTERS)\n');
fprintf('------------------------------------------------------------\n');

fprintf('\n--- Inertia Matrix B(q) (3x3) ---\n');
disp(B_symbolic);

fprintf('\n--- Kinetic Energy T(q, dq) ---\n');
disp(T_kinetic_symbolic);

fprintf('\n--- Potential Energy U(q) ---\n');
disp(U_potential_symbolic);

fprintf('\n------------------------------------------------------------\n');
fprintf('PART 2: EVALUATED NUMERICAL FORMULAS (q variables --> free)\n');
fprintf('------------------------------------------------------------\n');

fprintf('\n--- Evaluated Inertia Matrix B(q) (3x3) ---\n');
disp(B_numeric);

fprintf('\n--- Evaluated Kinetic Energy T(q, dq) ---\n');
disp(T_numeric);

fprintf('\n--- Evaluated Potential Energy U(q) ---\n');
disp(U_numeric);

fprintf('\n------------------------------------------------------------\n');
fprintf('PART 3: SAMPLE EVALUATION GIVEN A CONFIGURATION\n');
fprintf('------------------------------------------------------------\n');
fprintf('Sample Configuration q_sample  : [q1 = %.2f m, q2 = %.4f rad (pi/4), q3 = %.4f rad (pi/6)]^T\n', ...
    q1_sample, q2_sample, q3_sample);
fprintf('Sample Velocity dq_sample      : [dq1 = %.2f m/s, dq2 = %.2f rad/s, dq3 = %.2f rad/s]^T\n\n', ...
    dq1_sample, dq2_sample, dq3_sample);

fprintf('Evaluated Inertia Matrix B(q_sample) (3x3):\n');
disp(B_eval_sample);

fprintf('Evaluated Kinetic Energy T(q_sample, dq_sample)   : %.6f J\n', T_eval_sample);
fprintf('Evaluated Potential Energy U(q_sample)             : %.6f J\n', U_eval_sample);
fprintf('Total Mechanical Energy E(q_sample, dq_sample)     : %.6f J\n', E_eval_sample);
fprintf('============================================================\n\n');