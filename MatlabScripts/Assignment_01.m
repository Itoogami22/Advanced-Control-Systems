%% DH-Table, Direct Kinematics, Inverse Kinematics, Geometric and Analytical Jacobian
close all; clear; clc;

%% Display Robot and Check URDF Joint Limits
robot = importrobot('00001_PRR_zzx FILIPPO SEGALLA.urdf');
robot.DataFormat = 'row'; % Allow joint configurations as row vectors
showdetails(robot)

% Read the position limits directly from the imported URDF model.
joint_names = {};
joint_limits = zeros(0, 2);
for i = 1:numel(robot.Bodies)
    joint = robot.Bodies{i}.Joint;
    if ~strcmp(joint.Type, 'fixed')
        joint_names{end + 1} = joint.Name; %#ok<SAGROW>
        joint_limits(end + 1, :) = joint.PositionLimits; %#ok<SAGROW>
    end
end
if size(joint_limits, 1) ~= 3
    error('Expected three movable joints in the imported PRR robot.');
end

fprintf('\n--- Joint Position Limits from the URDF ---\n');
for i = 1:numel(joint_names)
    fprintf('%s: [%.4f, %.4f]\n', joint_names{i}, ...
        joint_limits(i, 1), joint_limits(i, 2));
end

% Figure 1: show the robot's home configuration.
home_config = homeConfiguration(robot);
figure;
show(robot, home_config);
title('Robot Home Configuration');
grid on; axis equal;
xlim([-1.0, 1.0]); ylim([-1.0, 1.0]); zlim([0.0, 1.4]);
rotate3d on;

fprintf("\n");

% Figure 2: compare two example configurations and confirm both are admissible.
example_configs = [-0.25, 5*pi/6, pi/3; ...
                   -0.40, pi/4, -pi/3];
example_names = {'Example configuration 1', 'Example configuration 2'};
limit_tolerance = 1e-9;
figure;
for k = 1:size(example_configs, 1)
    config = example_configs(k, :);
    within_limits = all(config >= joint_limits(:, 1).' - limit_tolerance) && ...
                    all(config <= joint_limits(:, 2).' + limit_tolerance);
    if within_limits
        limit_status = 'WITHIN LIMITS';
    else
        limit_status = 'OUTSIDE LIMITS';
    end
    fprintf('%s: q = [%.4f, %.4f, %.4f] -> %s\n', ...
        example_names{k}, config(1), config(2), config(3), ...
        limit_status);

    ax = subplot(1, 2, k);
    show(robot, config, 'Parent', ax);
    title(ax, example_names{k});
    grid(ax, 'on'); axis(ax, 'equal');
    xlim(ax, [-1.0, 1.0]); ylim(ax, [-1.0, 1.0]); zlim(ax, [0.0, 1.4]);
end
rotate3d on;
cameratoolbar('Show');

fprintf("------------------------------------------------------------\n");

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

% Cross-check the manually computed transformation with MATLAB's rigidBodyTree.
% Evaluate both matrices at the same numerical joint configuration.
q_check = [-0.2, 0.4, 0.6];
T_final_check = double(subs(T_final, q, q_check.'));
T_builtin = getTransform(robot, q_check, 'ee', 'world');
transform_equal = isequal(size(T_final_check), size(T_builtin)) && ...
    all(abs(T_final_check(:) - T_builtin(:)) < 1e-9);
if transform_equal
    fprintf("Transformation cross-check T_final vs getTransform(): MATCH\n");
else
    fprintf("Transformation cross-check T_final vs getTransform(): DIFFER\n");
end
fprintf("============================================================\n\n");

%% Inverse Kinematics (Closed-Form Analytical Solution for PRR)

% --- Numerical Example for Inverse Kinematics ---
px_val = 0.25;
py_val = 0.25;
pz_val = 0.9;
% -----------------------------------------------

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

% Display and filter analytical solutions against the URDF joint limits.
fprintf("\n--- Analytical Closed-Form Inverse Kinematics Solutions ---\n");
for k = 1:2
    fprintf("Solution Branch %d:\n", k);

    % Evaluate symbolic solutions with numerical values
    q1_num = double(subs(q1_ik_sol(k), [px py pz], [px_val py_val pz_val]));
    q2_num = double(subs(q2_ik_sol(k), [px py pz], [px_val py_val pz_val]));
    q3_num = double(subs(q3_ik_sol(k), [px py pz], [px_val py_val pz_val]));
    config_ik = [q1_num, q2_num, q3_num];
    within_limits = all(config_ik >= joint_limits(:, 1).' - limit_tolerance) && ...
                    all(config_ik <= joint_limits(:, 2).' + limit_tolerance);

    fprintf("\nNumerical Solution (px=%.2f, py=%.2f, pz=%.2f):\n", px_val, py_val, pz_val);
    fprintf("q1 = %.4f\n", q1_num);
    fprintf("q2 = %.4f\n", q2_num);
    fprintf("q3 = %.4f\n", q3_num);
    fprintf("\nSymbolic Form:\n");
    fprintf("q1 = %s\n", char(vpa(q1_ik_sol(k), 4)));
    fprintf("q2 = %s\n", char(vpa(q2_ik_sol(k), 4)));
    fprintf("q3 = %s\n", char(vpa(q3_ik_sol(k), 4)));
    if ~within_limits
        fprintf("\nJoint-limit check: REJECTED: this mathematical solution exceeds a URDF joint limit.\n");
        fprintf("------------------------------------------------------------\n");
        continue;
    end
    fprintf("\nJoint-limit check: ACCEPTED.\n");
    fprintf("------------------------------------------------------------\n");

    % --- Visualization of the Inverse Kinematics Solution ---
    figure;
    show(robot, config_ik);
    ax = gca;
    hold on;
    [sphere_x, sphere_y, sphere_z] = sphere(20);
    target_radius = 0.035;
    surf(px_val + target_radius * sphere_x, ...
         py_val + target_radius * sphere_y, ...
         pz_val + target_radius * sphere_z, ...
         'FaceColor', 'b', 'EdgeColor', 'none');
    text(px_val + 0.08, py_val + 0.03, pz_val, 'Target Point', 'Color', 'r', 'FontWeight', 'bold');
    title(sprintf('Inverse Kinematics - Branch %d', k));
    grid on; axis equal;
    xlim([-1.0, 1.0]); ylim([-1.0, 1.0]); zlim([0.0, 1.4]);
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    view(3);

    % Enable independent 3D rotation for this branch's axes
    rotate3d(ax, 'on');

    hold off;
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

% MATLAB returns [angular; linear], while JG is arranged as [linear; angular].
JG_check = double(subs(JG, q, q_check.'));
JG_builtin = geometricJacobian(robot, q_check, 'ee');
JG_builtin = [JG_builtin(4:6, :); JG_builtin(1:3, :)];
geometric_jacobian_equal = isequal(size(JG_check), size(JG_builtin)) && ...
    all(abs(JG_check(:) - JG_builtin(:)) < 1e-9);
if geometric_jacobian_equal
    fprintf("Geometric Jacobian cross-check JG vs geometricJacobian(): MATCH\n");
else
    fprintf("Geometric Jacobian cross-check JG vs geometricJacobian(): DIFFER\n");
end
fprintf("============================================================\n\n");

%% Analytical Jacobian
% Compute linear position component: JaP = d(p_e)/dq
% Manual computation of each element to avoid using the jacobian() function
JaP = sym(zeros(3, num_joints));
for i = 1:3      % Coordinates: px, py, pz
    for j = 1:num_joints  % Joints: q1, q2, q3
        JaP(i,j) = simplify(diff(pe(i), q(j)));
    end
end

% Extract End-Effector rotation matrix R (3x3)
R = T_final(1:3, 1:3);

% Extract orientation parameters using Euler angles representation (ZYZ sequence)
phi_euler   = simplify(atan2(R(2,3), R(1,3)));
theta_euler = simplify(atan2(sqrt(R(1,3)^2 + R(2,3)^2), R(3,3)));
psi_euler   = simplify(atan2(R(3,2), -R(3,1)));

Euler = [phi_euler; theta_euler; psi_euler];

% Compute orientation component: JaPhi = d(Euler)/dq
% Manual computation of each element to avoid using the jacobian() function
JaPhi = sym(zeros(3, num_joints));
for i = 1:3      % Angles: phi, theta, psi
    for j = 1:num_joints  % Joints: q1, q2, q3
        JaPhi(i,j) = simplify(diff(Euler(i), q(j)));
    end
end

% Complete Analytical Jacobian: Ja = [JaP; JaPhi]
Ja = simplify([JaP; JaPhi]);

% Remove abs() and dirac() artifacts from symbolic differentiation
Ja = subs(Ja, abs(sin(q3)), sin(q3));
Ja = subs(Ja, dirac(sin(q3)), 0);

% Display Analytical Jacobian matrix in simplified decimal format
fprintf("\n--- Analytical Jacobian Ja (6x%d) ---\n", num_joints);
disp(vpa(Ja, 4));
fprintf("============================================================\n\n");

% Convert the built-in angular velocity to ZYZ Euler-angle rates, then
% compare all six rows of the analytical and built-in Jacobians.
phi_check = double(subs(phi_euler, q, q_check.'));
theta_check = double(subs(theta_euler, q, q_check.'));
if abs(sin(theta_check)) < 1e-9
    fprintf("Full analytical Jacobian cross-check Ja: SKIPPED (ZYZ singularity)\n");
else
    EulerRateToAngular = [0, -sin(phi_check),  sin(theta_check) * cos(phi_check); ...
                          0,  cos(phi_check),  sin(theta_check) * sin(phi_check); ...
                          1,  0,               cos(theta_check)];
    Ja_builtin = [JG_builtin(1:3, :); ...
                  EulerRateToAngular \ JG_builtin(4:6, :)];
    Ja_check = double(subs(Ja, q, q_check.'));
    analytical_jacobian_equal = isequal(size(Ja_check), size(Ja_builtin)) && ...
        all(abs(Ja_check(:) - Ja_builtin(:)) < 1e-9);
    if analytical_jacobian_equal
        fprintf("Full analytical Jacobian cross-check Ja: MATCH\n");
    else
        fprintf("Full analytical Jacobian cross-check Ja: DIFFER\n");
    end
end
fprintf("The built-in angular velocity is converted to ZYZ Euler-angle rates for comparison.\n");
fprintf("============================================================\n\n");
