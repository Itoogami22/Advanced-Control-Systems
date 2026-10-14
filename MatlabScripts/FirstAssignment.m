%% DH-Table, Direct Kinematics, Inverse Kinematics, Geometric and Analytical Jacobian
close all; clc;

% Symbolic variables for the computations
syms a alpha d theta real                           % Generic DH-Table parameters
syms q1 q2 q3 real                                  % Active joint variables (q_i)
syms a1 a2 a3 alpha1 alpha2 alpha3 d1 d2 d3 real    % Kinematic robot parameters (link lengths and offsets)
syms px py pz phi real                              % Inverse Kinematics parameters

%% Direct Kinematics
% Generic homogeneous transformation matrix A_i based on standard DH parameters
A = [cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha) a*cos(theta); ...
     sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta); ...
     0           sin(alpha)             cos(alpha)            d; ...
     0           0                      0                     1];

% Base frame transformation matrix (World frame to Frame 0)
% Default is the identity matrix eye(4) when no base offset/rotation exists
TB0 = eye(4);

% End-Effector frame transformation matrix (Frame n to Tool frame)
% Default is the identity matrix eye(4) when no tool offset/rotation exists
EE_baseTransformation = eye(4);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%                                                                     %
%  DH_table = [ a | alpha | d | theta ] --> Structure of the DH-Table %
%                                                                     %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Vector of joint variables matching the DH table dimensions
q = [q1; q2; q3];

% Joint types: 'R' for Revolute (theta is variable), 'P' for Prismatic (d is variable)
joint_types = ['R'; 'R'; 'R']';

% Denavit-Hartenberg parameter table [a, alpha, d, theta]
DH_table = [ a1  alpha1  d1  q1;
             a2  alpha2  d2  q2;
             a3  alpha3  d3  q3 ];

% Compute cumulative transformation matrices T_0i
num_joints = size(DH_table, 1);
T = cell(num_joints, 1);
T_current = TB0; % Start cumulative multiplication from base frame

for i = 1:num_joints
    % Substitute link i DH parameters into generic matrix A
    A_i = subs(A, [a alpha d theta], DH_table(i,:));
    
    % Compute cumulative homogeneous transformation: T_0i = T_0(i-1) * A_i
    T_current = simplify(T_current * A_i);
    T{i} = T_current; % Store transformation T_0^i in cell array
end

% Apply End-Effector transformation to obtain final matrix T_0n
T_final = simplify(T{end} * EE_baseTransformation);

fprintf("\n--- Base Transformation Matrix ---\n");
fprintf("Matrix TB0 (World Frame -> Joint 0 Frame):\n");
disp(TB0);
fprintf("------------------------------------------------------------\n");

% Display intermediate transformation matrices T_0i
fprintf("\n--- Intermediate Transformation Matrices ---\n");
for i = 1:num_joints
    fprintf("Matrix TB%d (Base -> Joint %d):\n", i, i);
    disp(T{i}); % Display T_0i matrix
    fprintf("------------------------------------------------------------\n");
end

% Display final transformation matrix T_0^n
fprintf("\nMatrix T0E Final (Base -> End-Effector):\n");
disp(T_final);
fprintf("============================================================\n\n");

%% Inverse Kinematics
% Desired homogeneous transformation matrix T_des
T_des = [ cos(phi) -sin(phi)  0  px;
          sin(phi)  cos(phi)  0  py;
          0           0       1  pz;
          0           0       0  1 ];

% Symbolic extraction of individual transformation matrices A_i
A1 = subs(A, [a alpha d theta], DH_table(1,:));
A2 = subs(A, [a alpha d theta], DH_table(2,:));
A3 = subs(A, [a alpha d theta], DH_table(3,:));

% Analytical inverse of A1 using rigid body property: inv(A1) = [R1'    -R1'*p1;
%                                                                0 0 0   1]
R1 = A1(1:3, 1:3);
p1 = A1(1:3, 4);
A1_inv = [R1', -R1'*p1; 0 0 0 1];

% Matrix decomposition: LHS = A1_inv * T_des, RHS = A2 * A3
LHS = simplify(A1_inv * T_des);
RHS = simplify(A2 * A3);

% Kinematic decoupling: position of Joint 3 origin (subtracting link 3 contribution)
px2 = px - a3 * cos(phi); % Coordinate X of Frame 2 --> Joint 3
py2 = py - a3 * sin(phi); % Coordinate Y of Frame 2 --> Joint 3

% Closed-form analytical solution for q2 (Law of Cosines on reduced 2R arm)
cos_q2_val = simplify((px2^2 + py2^2 - a1^2 - a2^2) / (2 * a1 * a2));
sin_q2_pos = sqrt(1 - cos_q2_val^2);
sin_q2_neg = -sqrt(1 - cos_q2_val^2);

% Vector of solutions for q2 (Branch 1: Elbow-Up, Branch 2: Elbow-Down)
q2_ik_sol = [atan2(sin_q2_pos, cos_q2_val); 
    atan2(sin_q2_neg, cos_q2_val)];

% Closed-form analytical solution for q1
q1_ik_sol = [atan2(py2, px2) - atan2(a2 * sin_q2_pos, a1 + a2 * cos_q2_val);
    atan2(py2, px2) - atan2(a2 * sin_q2_neg, a1 + a2 * cos_q2_val)];

% Closed-form analytical solution for q3 (derived from total orientation phi)
q3_ik_sol = [phi - q1_ik_sol(1) - q2_ik_sol(1);
    phi - q1_ik_sol(2) - q2_ik_sol(2)];

% Display analytical inverse kinematics solutions
fprintf("\n--- Analytical Closed-Form Inverse Kinematics Solutions ---\n");
for k = 1:2
    fprintf("Solution Branch %d:\n", k);
    fprintf("q1 = %s\n", char(q1_ik_sol(k)));
    fprintf("q2 = %s\n", char(q2_ik_sol(k)));
    fprintf("q3 = %s\n", char(q3_ik_sol(k)));
    fprintf("------------------------------------------------------------\n");
end
fprintf("============================================================\n\n");

%% Geometric Jacobian
pe = T_final(1:3, 4);           % End-Effector position vector p_e = [px; py; pz]
JG = sym(zeros(6, num_joints)); % Initialize Geometric Jacobian matrix (6 rows x N joints)

z_prev = sym([0 0 1]');         % Initial z-axis direction vector (z_0) in base frame
p_prev = sym([0 0 0]');         % Initial frame origin position vector (p_0) in base frame

% Column-by-column computation of the Geometric Jacobian JG = [Jv; Jw]
for i = 1:num_joints
    if i == 1
        % For the first joint, use base frame z_0 axis and origin p_0
        zi_prev = z_prev;
        pi_prev = p_prev;
    else
        % For subsequent joints, extract z_{i-1} and p_{i-1} from matrix T_0{i-1}
        zi_prev = T{i-1}(1:3, 3); % Third column of T_0{i-1} (z-axis)
        pi_prev = T{i-1}(1:3, 4); % Fourth column of T_0{i-1} (origin p)
    end

    % Compute linear velocity (Jv) and angular velocity (Jw) contributions
    if joint_types(i) == 'R'
        % Revolute Joint ('R'):
        % Jv_i = z_{i-1} x (p_e - p_{i-1})
        % Jw_i = z_{i-1}
        Jv_i = simplify(cross(zi_prev, (pe - pi_prev)));
        Jw_i = zi_prev;
    elseif joint_types(i) == 'P'
        % Prismatic Joint ('P'):
        % Jv_i = z_{i-1}
        % Jw_i = [0; 0; 0]
        Jv_i = zi_prev;
        Jw_i = sym([0; 0; 0]);
    else
        error('Invalid joint type. Use ''R'' or ''P''.');
    end

    % Assign the i-th column to the full Geometric Jacobian matrix
    JG(:, i) = [Jv_i; Jw_i];
end

% Algebraic simplification of the Geometric Jacobian
JG = simplify(JG);

% Display Geometric Jacobian matrix
fprintf("\n--- Geometric Jacobian JG (6x%d) ---\n", num_joints);
disp(JG);
fprintf("============================================================\n\n");

%% Analytical Jacobian
% Compute linear position component of Analytical Jacobian: J_aP = d(p_e)/dq
JaP = simplify(jacobian(pe, q));

% Extract End-Effector rotation matrix R (3x3) from final transformation matrix
R = T_final(1:3, 1:3);

% Extract orientation parameters using Euler angles representation (ZYZ sequence)
phi_euler     = atan2(R(2,3), R(1,3));
theta_euler = atan2(sqrt(R(1,3)^2 + R(2,3)^2), R(3,3));
psi_euler     = atan2(R(3,2), -R(3,1));

% Vector of End-Effector orientation variables
Euler = [phi_euler; theta_euler; psi_euler];

% Compute orientation component of Analytical Jacobian: J_aPhi = d(Euler)/dq
JaPhi = simplify(jacobian(Euler, q));

% Stack linear and angular components to form complete Analytical Jacobian (Ja = [JaP; JaPhi])
Ja = simplify([JaP; JaPhi]);

% Display Analytical Jacobian matrix
fprintf("\n--- Analytical Jacobian JG (6x%d) ---\n", num_joints);
disp(Ja);
fprintf("============================================================\n\n");