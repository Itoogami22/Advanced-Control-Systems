%% ========================================================================
%  ADVANCED CONTROL SYSTEMS - SECOND ASSIGNMENT (TRAJECTORY VERSION)
% ========================================================================
% This script calculates the kinetic and potential energies, as well as the
% generalized inertia matrix B(q), for a 3-DOF PRR manipulator.
% It uses a Quintic Polynomial trajectory to derive realistic joint
% positions and velocities for numerical evaluation and animation.
% ========================================================================

close all; clear; clc;

%% 1. Robot Import and Setup
% Import the robot model from the URDF file (Robot DNA containing masses and dimensions)
robot = importrobot('00001_PRR_zzx FILIPPO SEGALLA.urdf');
% Set DataFormat to 'column' to ensure joint states are vectors (3x1) for matrix math
robot.DataFormat = 'column';
showdetails(robot)

%% 2. Symbolic Variables Definition
% We use symbolic variables to derive general mathematical formulas before plugging in numbers
num_joints = 3;

% Standard DH parameters for general matrix construction
syms a alpha d theta real

% Generalized coordinates (joint positions q) and their time derivatives (velocities dq)
syms q1 q2 q3 real
syms dq1 dq2 dq3 real

% Physical parameters of the robot links (Masses, Gravity, and Inertia Moments)
syms m1 m2 m3 real                                  % Link masses (kg)
syms g real                                         % Gravitational acceleration (m/s^2)
syms I1xx I1yy I1zz real                            % Link 1 inertia tensor components
syms I2xx I2yy I2zz real                            % Link 2 inertia tensor components
syms I3xx I3yy I3zz real                            % Link 3 inertia tensor components

% Group variables into vectors for streamlined matrix operations
q  = [q1; q2; q3];
dq = [dq1; dq2; dq3];
m  = [m1, m2, m3];

% Gravity vector defined in the base frame Sigma_0 (pointing downwards along Z)
g0 = [0; 0; -g];

% Define local inertia tensors relative to the Center of Mass (CoM) of each link
I_local = cell(num_joints, 1);
I_local{1} = diag([I1xx, I1yy, I1zz]);
I_local{2} = diag([I2xx, I2yy, I2zz]);
I_local{3} = diag([I3xx, I3yy, I3zz]);

%% 3. Direct Kinematics and DH Transformation Matrices
% Standard Denavit-Hartenberg (DH) transformation matrix template
A = [cos(theta), -sin(theta)*cos(alpha),  sin(theta)*sin(alpha), a*cos(theta); ...
     sin(theta),  cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta); ...
     0,           sin(alpha),             cos(alpha),            d; ...
     0,           0,                      0,                     1];

% Transformation from World Frame to Robot Base Frame (Sigma_0)
TB0 = [ 0, -1,  0,  0;
        1,  0,  0,  0;
        0,  0,  1,  0.65;
        0,  0,  0,  1 ];

% DH Parameter Table: [a, alpha, d, theta]
DH_table = [ 0,    0,     q1,   0;       % Link 1: Prismatic
             0,    pi/2,  0.4,  q2;     % Link 2: Revolute
             0.4,  0,     0,    q3 ];    % Link 3: Revolute

% Identify joint types to apply different Jacobian rules (Prismatic vs Revolute)
joint_types = ['P'; 'R'; 'R'];

T = cell(num_joints, 1); % Stores cumulative transformations T_i^0
R = cell(num_joints, 1); % Stores rotation matrices R_i^0
T_curr = sym(TB0);

% Compute cumulative transformation matrices iteratively: T_i = T_{i-1} * A_i
for i = 1:num_joints
    A_i = subs(A, [a, alpha, d, theta], DH_table(i,:));
    T_curr = simplify(T_curr * A_i);
    T{i} = T_curr;
    R{i} = T{i}(1:3, 1:3);
end

%% 4. Center of Mass (CoM) Positions in Base Frame
% CoM positions defined in the local link frame [x; y; z; 1]
p_l_local = {[0; 0; -0.25; 1], [0; 0; -0.20; 1], [-0.20; 0; 0; 1]};
p_l = cell(num_joints, 1);

% Transform local CoM coordinates to the base frame Sigma_0 using T matrices
for i = 1:num_joints
    p_temp = simplify(T{i} * p_l_local{i});
    p_l{i} = p_temp(1:3);
end

%% 5. Partial Jacobians for Link CoMs (Linear and Angular)
JP_l = cell(num_joints, 1); % Linear Jacobians (How CoM moves linearly)
JO_l = cell(num_joints, 1); % Angular Jacobians (How CoM rotates)

z_axes = cell(num_joints+1, 1);    % Z-axis of each coordinate frame
p_origins = cell(num_joints+1, 1);  % Origin of each coordinate frame

z_axes{1} = sym(TB0(1:3, 3));
p_origins{1} = sym(TB0(1:3, 4));

for i = 1:num_joints
    z_axes{i+1} = T{i}(1:3, 3);
    p_origins{i+1} = T{i}(1:3, 4);
end

% Compute Jacobian columns based on joint type
for i = 1:num_joints
    JP_i = sym(zeros(3, num_joints));
    JO_i = sym(zeros(3, num_joints));

    for j = 1:i
        z_prev = z_axes{j};
        p_prev = p_origins{j};

        if joint_types(j) == 'P'
            % Prismatic joint: purely translational velocity along Z
            JP_i(:, j) = z_prev;
            JO_i(:, j) = [0; 0; 0];
        else % Revolute joint: velocity is cross product of axis and distance to CoM
            JP_i(:, j) = cross(z_prev, (p_l{i} - p_prev));
            JO_i(:, j) = z_prev;
        end
    end

    JP_l{i} = simplify(JP_i);
    JO_l{i} = simplify(JO_i);
end

%% 6. Generalized Inertia Matrix B(q) and Energy Formulas
% Compute Total Inertia Matrix B(q) based on the sum of translational and rotational components
B_symbolic = sym(zeros(num_joints, num_joints));

for i = 1:num_joints
    B_symbolic = B_symbolic + m(i) * (JP_l{i}' * JP_l{i}) + JO_l{i}' * R{i} * I_local{i} * R{i}' * JO_l{i};
end
B_symbolic = simplify(B_symbolic);

% Total Kinetic Energy T = 1/2 * dq' * B(q) * dq
T_kinetic_symbolic = simplify(0.5 * dq' * B_symbolic * dq);

% Total Potential Energy U = sum of (m_i * g * height_i)
U_potential_symbolic = sym(0);
for i = 1:num_joints
    U_potential_symbolic = U_potential_symbolic - m(i) * (g0' * p_l{i});
end
U_potential_symbolic = simplify(U_potential_symbolic);

%% 7. Numerical Value Substitution (Exact Rational Constants)
% Physical constants (using fractions for maximum symbolic precision)
m1_val = sym(1); m2_val = sym(1); m3_val = sym(1);
g_val  = sym('981/100'); % 9.81 m/s^2

% Inertia tensors for cylinders/bars based on physical dimensions
I1xx_val = sym('1/12') * m1_val * (sym('36/1000')^2 + sym('1/2')^2);
I1yy_val = sym('1/12') * m1_val * (sym('36/1000')^2 + sym('1/2')^2);
I1zz_val = sym('1/12') * m1_val * (sym('36/1000')^2 + sym('36/1000')^2);

I2xx_val = sym('1/12') * m2_val * (3*sym('19/1000')^2 + sym('4/10')^2);
I2yy_val = sym('1/12') * m2_val * (3*sym('19/1000')^2 + sym('4/10')^2);
I2zz_val = sym('1/2')  * m2_val * (sym('19/1000')^2);

I3xx_val = sym('1/2')  * m3_val * (sym('19/1000')^2);
I3yy_val = sym('1/12') * m3_val * (3*sym('19/1000')^2 + sym('4/10')^2);
I3zz_val = sym('1/12') * m3_val * (3*sym('19/1000')^2 + sym('4/10')^2);

% Map symbolic variables to their numerical values
sym_vars = [m1, m2, m3, g, I1xx, I1yy, I1zz, I2xx, I2yy, I2zz, I3xx, I3yy, I3zz];
num_vals = [m1_val, m2_val, m3_val, g_val, I1xx_val, I1yy_val, I1zz_val, I2xx_val, I2yy_val, I2zz_val, I3xx_val, I3yy_val, I3zz_val];

% Final numerical symbolic expressions
B_num_exact = subs(B_symbolic, sym_vars, num_vals);
T_num_exact = subs(T_kinetic_symbolic, sym_vars, num_vals);
U_num_exact = subs(U_potential_symbolic, sym_vars, num_vals);

%% 8. Trajectory Generation (Realistic Values)
% Define movement time window
t_start = 0;
t_end = 2;
t_sample = 1.0; % Mid-point evaluation

% Start and end positions (matching URDF limits)
q1_0 = -0.5; q1_f = 0.0; % Prismatic
q2_0 = 0; q2_f = pi/4;    % Revolute
q3_0 = 0; q3_f = pi/6;     % Revolute

% Quintic Polynomial function: ensures smooth motion with zero velocity and acceleration at start and end
% Formula: q(t) = a0 + a1*t + a2*t^2 + a3*t^3 + a4*t^4 + a5*t^5
get_quintic = @(q0, qf, tf, t) [ ...
    q0 + 10*(qf-q0)/tf^3*t^3 - 15*(qf-q0)/tf^4*t^4 + 6*(qf-q0)/tf^5*t^5, ...
    15*(qf-q0)/tf^3*t^2 - 24*(qf-q0)/tf^4*t^3 + 15*(qf-q0)/tf^5*t^4 ...
];

% Calculate the specific state (pos and vel) at the sample time
res1 = get_quintic(q1_0, q1_f, t_end, t_sample);
res2 = get_quintic(q2_0, q2_f, t_end, t_sample);
res3 = get_quintic(q3_0, q3_f, t_end, t_sample);

q_sample = [res1(1); res2(1); res3(1)];
dq_sample = [res1(2); res2(2); res3(2)];

%% 9. Numerical Evaluation at Sample State
% Use 'double' to convert symbolic results to floating point numbers
B_eval_sample = double(subs(B_num_exact, [q1; q2; q3], q_sample));
T_eval_sample = double(subs(T_num_exact, [q1; q2; q3; dq1; dq2; dq3], [q_sample; dq_sample]));
U_eval_sample = double(subs(U_num_exact, [q1; q2; q3], q_sample));
E_eval_sample = T_eval_sample + U_eval_sample;

%% 10. Display Detailed Results
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
disp(vpa(expand(B_num_exact), 4));
fprintf('\n--- Evaluated Kinetic Energy T(q, dq) ---\n');
disp(vpa(expand(T_num_exact), 4));
fprintf('\n--- Evaluated Potential Energy U(q) ---\n');
disp(vpa(expand(U_num_exact), 4));

fprintf('\n------------------------------------------------------------\n');
fprintf('PART 3: TRAJECTORY-BASED REALISTIC EVALUATION\n');
fprintf('------------------------------------------------------------\n');
fprintf('Sample Time t = %.2f s (of %.2f s total)\n', t_sample, t_end);
fprintf('Configuration q_sample: [%.4f, %.4f, %.4f]^T\n', q_sample);
fprintf('Velocities dq_sample: [%.4f, %.4f, %.4f]^T\n\n', dq_sample);
fprintf('Evaluated Inertia Matrix B(q_sample):\n');
disp(B_eval_sample);
fprintf('Kinetic Energy T: %.6f J\n', T_eval_sample);
fprintf('Potential Energy U: %.6f J\n', U_eval_sample);
fprintf('Total Energy E: %.6f J\n', E_eval_sample);
fprintf('============================================================\n\n');

%% 11. Trajectory Animation
figure('Name', 'PRR Robot Trajectory Animation', 'NumberTitle', 'off');
t_vec = linspace(t_start, t_end, 60); % Define 60 frames for smooth motion

for t = t_vec
    % Calculate q for the current time t using the quintic polynomial
    res1 = get_quintic(q1_0, q1_f, t_end, t);
    res2 = get_quintic(q2_0, q2_f, t_end, t);
    res3 = get_quintic(q3_0, q3_f, t_end, t);
    q_current = [res1(1); res2(1); res3(1)];

    % Update the robot's pose in the plot
    show(robot, q_current, 'PreservePlot', false, 'Frames', 'off');

    title(sprintf('PRR Manipulator Motion - Time: %.2f s', t), 'FontSize', 12);
    axis equal; grid on;
    xlim([-1.5, 1.5]); ylim([-1.5, 1.5]); zlim([-0.2, 2.0]);

    drawnow;    % Force immediate update of the figure
    pause(0.05); % Control animation speed
end
