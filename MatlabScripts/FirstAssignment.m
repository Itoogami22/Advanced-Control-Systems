%% DH-Table, Direct Kinematics, Inverse Kinematics, Geometric and Analytical Jacobian
close all; clear; clc;

%% Display Robot
robot = importrobot('00001_PRR_zzx FILIPPO SEGALLA.urdf');
showdetails(robot)

fig = figure;
config = homeConfiguration(robot);
show(robot, config);

% Enable 3D plot controls (zoom, pan, rotate)
rotate3d on;              % Enable 3D dynamic rotation with left-click drag
cameratoolbar('Show');    % Display camera toolbar at the top
grid on;                  % Enable grid lines for spatial reference
axis equal;               % Keep proportional axis scaling

% Set larger axis limits to center the robot in an expanded grid
xlim([-1.0, 1.0]);        % Expand X-axis limits
ylim([-1.0, 1.0]);        % Expand Y-axis limits
zlim([ 0.0, 1.4]);        % Expand Z-axis limits

%% Symbolic variables for the computations
syms a alpha d theta real                           % Generic DH-Table parameters
syms q1 q2 q3 real                                  % Active joint variables (q_i)
syms px py pz real                                  % Target position for Inverse Kinematics

%% Direct Kinematics
% Generic homogeneous transformation matrix A_i based on standard DH parameters
A = [cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha) a*cos(theta); ...
     sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta); ...
     0           sin(alpha)             cos(alpha)            d; ...
     0           0                      0                     1];

% Base frame transformation matrix (World frame to Frame 0)
TB0 = [ 0 -1  0  0;
        1  0  0  0;
        0  0  1  0.65;
        0  0  0  1 ];

% End-Effector frame transformation matrix (Frame 3 to Tool/EE frame)
EE_baseTransformation = [ 0 1 0 0;
                          0 0 1 0;
                          1 0 0 0;
                          0 0 0 1 ];

% Vector of joint variables
q = [q1; q2; q3];

% Joint types: 'P' for Prismatic, 'R' for Revolute
joint_types = ['P'; 'R'; 'R']';

% Denavit-Hartenberg parameter table [a, alpha, d, theta]
DH_table = [ 0    0     q1   0;
             0    pi/2  0.4  q2;
             0.4  0     0    q3 ];

% Compute cumulative transformation matrices T_0i
num_joints = size(DH_table, 1);
T = cell(num_joints, 1);
T_current = sym(TB0); % Cast to symbolic matrix for cumulative multiplication

for i = 1:num_joints
    % Substitute link i DH parameters into generic matrix A
    A_i = simplify(subs(A, [a alpha d theta], DH_table(i,:)));
    
    % Compute cumulative homogeneous transformation: T_0i = T_0(i-1) * A_i
    T_current = simplify(T_current * A_i);
    T{i} = T_current; % Store transformation T_0^i in cell array
end

% Apply End-Effector transformation to obtain final matrix T_final
T_final = simplify(T{end} * EE_baseTransformation);

fprintf("\n--- Base Transformation Matrix ---\n");
fprintf("Matrix TB0 (World Frame -> Joint 0 Frame):\n");
disp(vpa(TB0, 4));
fprintf("------------------------------------------------------------\n");

% Display intermediate transformation matrices in simplified decimal format
fprintf("\n--- Intermediate Transformation Matrices ---\n");
for i = 1:num_joints
    fprintf("Matrix T0%d (Base -> Joint %d):\n", i, i);
    disp(vpa(T{i}, 4));
    fprintf("------------------------------------------------------------\n");
end

% Display final transformation matrix T0E in simplified decimal format
fprintf("\nMatrix TBE Final (Base -> End-Effector):\n");
disp(vpa(T_final, 4));
fprintf("============================================================\n\n");

%% Inverse Kinematics (Closed-Form Analytical Solution for PRR)
r_xy = simplify(sqrt(px^2 + py^2));
cos_q3_val = simplify(r_xy / 0.4);

% Solution Branch 1 (Elbow-Up)
sin_q3_b1 = simplify(sqrt(1 - cos_q3_val^2));
q3_b1 = simplify(atan2(sin_q3_b1, cos_q3_val));
q1_b1 = simplify(pz - 1.05 - 0.4 * sin_q3_b1);
q2_b1 = simplify(atan2(-px, py));

% Solution Branch 2 (Elbow-Down)
sin_q3_b2 = simplify(-sqrt(1 - cos_q3_val^2));
q3_b2 = simplify(atan2(sin_q3_b2, cos_q3_val));
q1_b2 = simplify(pz - 1.05 - 0.4 * sin_q3_b2);
q2_b2 = simplify(atan2(-px, py));

q1_ik_sol = [q1_b1; q1_b2];
q2_ik_sol = [q2_b1; q2_b2];
q3_ik_sol = [q3_b1; q3_b2];

% Display analytical inverse kinematics solutions with decimal constants
fprintf("\n--- Analytical Closed-Form Inverse Kinematics Solutions ---\n");
for k = 1:2
    fprintf("Solution Branch %d:\n", k);
    fprintf("q1 = %s\n", char(vpa(q1_ik_sol(k), 4)));
    fprintf("q2 = %s\n", char(vpa(q2_ik_sol(k), 4)));
    fprintf("q3 = %s\n", char(vpa(q3_ik_sol(k), 4)));
    fprintf("------------------------------------------------------------\n");
end
fprintf("============================================================\n\n");

%% Geometric Jacobian
pe = T_final(1:3, 4);           % End-Effector position vector p_e = [px; py; pz]
JG = sym(zeros(6, num_joints)); % Initialize Geometric Jacobian matrix (6 x N)

% Extract z_0 and p_0 from TB0
z_prev = sym(TB0(1:3, 3));
p_prev = sym(TB0(1:3, 4));

% Column-by-column computation of the Geometric Jacobian JG = [Jv; Jw]
for i = 1:num_joints
    if i == 1
        zi_prev = z_prev;
        pi_prev = p_prev;
    else
        zi_prev = T{i-1}(1:3, 3);
        pi_prev = T{i-1}(1:3, 4);
    end

    if joint_types(i) == 'R'
        Jv_i = simplify(cross(zi_prev, (pe - pi_prev)));
        Jw_i = zi_prev;
    elseif joint_types(i) == 'P'
        Jv_i = zi_prev;
        Jw_i = sym([0; 0; 0]);
    else
        error('Invalid joint type. Use ''R'' or ''P''.');
    end

    JG(:, i) = [Jv_i; Jw_i];
end

JG = simplify(JG);

% Display Geometric Jacobian matrix in simplified decimal format
fprintf("\n--- Geometric Jacobian JG (6x%d) ---\n", num_joints);
disp(vpa(JG, 4));
fprintf("============================================================\n\n");

%% Analytical Jacobian
% Compute linear position component: JaP = d(p_e)/dq
JaP = simplify(jacobian(pe, q));

% Extract End-Effector rotation matrix R (3x3)
R = T_final(1:3, 1:3);

% Extract orientation parameters using Euler angles representation (ZYZ sequence)
phi_euler   = simplify(atan2(R(2,3), R(1,3)));
theta_euler = simplify(atan2(sqrt(R(1,3)^2 + R(2,3)^2), R(3,3)));
psi_euler   = simplify(atan2(R(3,2), -R(3,1)));

Euler = [phi_euler; theta_euler; psi_euler];

% Compute orientation component: JaPhi = d(Euler)/dq
JaPhi = simplify(jacobian(Euler, q));

% Complete Analytical Jacobian: Ja = [JaP; JaPhi]
Ja = simplify([JaP; JaPhi]);

% Display Analytical Jacobian matrix in simplified decimal format
fprintf("\n--- Analytical Jacobian Ja (6x%d) ---\n", num_joints);
disp(vpa(Ja, 4));
fprintf("============================================================\n\n");