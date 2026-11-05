%% ========================================================================
%  ADVANCED CONTROL SYSTEMS - SECOND ASSIGNMENT (LIVE SYNCHRONIZED VERSION)
% ========================================================================

close all; clear; clc;

%% 1. Robot Import and Setup
% Load the robot model from the URDF file to handle kinematics and visualization.
robot = importrobot('00001_PRR_zzx FILIPPO SEGALLA.urdf');
% Set the data format to 'column' to ensure joint configurations are handled as vectors.
robot.DataFormat = 'column';

%% 2. Symbolic Variables and Energy Formulas
% Define the number of degrees of freedom (DOF).
num_joints = 3;

% Define symbolic variables for DH parameters and robot geometry.
syms a alpha d theta real

% Define joint positions (q) and velocities (dq).
syms q1 q2 q3 real
syms dq1 dq2 dq3 real

% Define mass (m) and gravity (g).
syms m1 m2 m3 real
syms g real

% Define the moment of inertia tensors for each link (diagonal matrices).
syms I1xx I1yy I1zz real
syms I2xx I2yy I2zz real
syms I3xx I3yy I3zz real

% Vectorize joint positions, velocities, and masses for cleaner matrix operations.
q  = [q1; q2; q3];
dq = [dq1; dq2; dq3];
m  = [m1, m2, m3];
g0 = [0; 0; -g]; % Gravity vector acting in the -Z direction.

% Assign the local baricentric inertia tensors as diagonal matrices for each link.
I_local = cell(num_joints, 1);
I_local{1} = diag([I1xx, I1yy, I1zz]);
I_local{2} = diag([I2xx, I2yy, I2zz]);
I_local{3} = diag([I3xx, I3yy, I3zz]);

% Define the general Denavit-Hartenberg (DH) transformation matrix.
A = [cos(theta), -sin(theta)*cos(alpha),  sin(theta)*sin(alpha), a*cos(theta); ...
     sin(theta),  cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta); ...
     0,           sin(alpha),             cos(alpha),            d; ...
     0,           0,                      0,                     1];

% Base transformation matrix to position the robot correctly in the world frame.
TB0 = [ 0, -1,  0,  0; 
        1,  0,  0,  0; 
        0,  0,  1,  0.65; 
        0,  0,  0,  1 ];

% DH table: [a, alpha, d, theta] for each link.
DH_table = [ 0,   0,    q1,  0; 
             0,   pi/2, 0.4, q2; 
             0.4, 0,    0,   q3 ];
% Specify joint types: P = Prismatic, R = Revolute.
joint_types = ['P'; 'R'; 'R'];

% Compute the transformation matrices T and rotation matrices R for each link.
T = cell(num_joints, 1);
R = cell(num_joints, 1);
T_curr = sym(TB0);

for i = 1:num_joints
    % Substitute the DH parameters for the i-th link into the general matrix A.
    A_i = subs(A, [a, alpha, d, theta], DH_table(i,:));
    % Accumulate transformations to get the pose of each link relative to the base.
    T_curr = simplify(T_curr * A_i);
    T{i} = T_curr;
    % Extract the 3x3 rotation matrix.
    R{i} = T{i}(1:3, 1:3);
end

% Define the center of mass (COM) for each link in its local frame.
p_l_local = {[0; 0; -0.25; 1], [0; 0; -0.20; 1], [-0.20; 0; 0; 1]};
p_l = cell(num_joints, 1);

for i = 1:num_joints
    % Transform the COM position to the world frame.
    p_l{i} = simplify(T{i} * p_l_local{i});
    p_l{i} = p_l{i}(1:3);
end

%% --- COMPUTATION OF STEINER INERTIA TENSORS FOR REVOLUTE JOINTS ---
I_steiner = cell(num_joints, 1);
for i = 1:num_joints
    if joint_types(i) == 'R'
        c_i = p_l_local{i}(1:3); % CoM position relative to the joint origin
        % Steiner parallel axis tensor: S = m * ( (c'*c)*I - c*c' )
        S_matrix = m(i) * ((c_i' * c_i) * eye(3) - (c_i * c_i'));
        I_steiner{i} = simplify(I_local{i} + S_matrix);
    else
        I_steiner{i} = I_local{i}; % Keep baricentric inertia for prismatic link
    end
end

% Derive the linear velocity (JP) and angular velocity (JO) Jacobians for each COM.
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

% Compute Jacobian columns based on joint type
for i = 1:num_joints
    JP_i = sym(zeros(3, num_joints));
    JO_i = sym(zeros(3, num_joints));
    
    for j = 1:i
        z_prev = z_axes{j};
        p_prev = p_origins{j};
        
        if joint_types(j) == 'P'
            % Prismatic joint: linear velocity contribution, no angular contribution.
            JP_i(:, j) = z_prev; JO_i(:, j) = [0; 0; 0];
        else
            % Revolute joint: linear velocity (cross product) and angular velocity contribution.
            JP_i(:, j) = cross(z_prev, (p_l{i} - p_prev)); JO_i(:, j) = z_prev;
        end
    end
    
    JP_l{i} = simplify(JP_i); 
    JO_l{i} = simplify(JO_i);
end

% Calculate the Symbolic Inertia Matrix B(q) using I_steiner for revolute joints.
B_symbolic = sym(zeros(num_joints, num_joints));

for i = 1:num_joints
    B_symbolic = simplify(B_symbolic + m(i) * (JP_l{i}' * JP_l{i}) + JO_l{i}' * R{i} * I_steiner{i} * R{i}' * JO_l{i}); % I_steiner used here
end

% Define the total Kinetic Energy T = 0.5 * dq' * B * dq.
T_kinetic_symbolic = simplify(0.5 * dq' * B_symbolic * dq);

% Define the total Potential Energy U = sum(-m_i * g' * p_com_i).
U_potential_symbolic = sym(0);
for i = 1:num_joints
    U_potential_symbolic = U_potential_symbolic - m(i) * (g0' * p_l{i});
end

% --- PRINT SYMBOLIC OUTPUTS TO CONSOLE ---
fprintf('\n============================================================\n');
fprintf('SYMBOLIC DERIVATIONS\n');
fprintf('============================================================\n');
fprintf('\n1. Inertia Matrix B(q):\n');
disp(simplify(B_symbolic));
fprintf('\n2. Kinetic Energy T(q, dq):\n');
disp(T_kinetic_symbolic);
fprintf('\n3. Potential Energy U(q):\n');
disp(simplify(U_potential_symbolic));
fprintf('============================================================\n\n');

% Numerical values for masses and inertia tensors for simulation.
m1_val = sym(1); m2_val = sym(1); m3_val = sym(1); 
g_val = sym('981/100');
I1xx_val = sym('1/12')*m1_val*(sym('36/1000')^2+sym('1/2')^2); 
I1yy_val = I1xx_val; 
I1zz_val = sym('1/12')*m1_val*(sym('36/1000')^2+sym('36/1000')^2);

I2xx_val = sym('1/12')*m2_val*(3*sym('19/1000')^2+sym('4/10')^2); 
I2yy_val = I2xx_val; 
I2zz_val = sym('1/2')*m2_val*(sym('19/1000')^2);

I3xx_val = sym('1/2')*m3_val*(sym('19/1000')^2); 
I3yy_val = sym('1/12')*m3_val*(3*sym('19/1000')^2+sym('4/10')^2); 
I3zz_val = I3yy_val;

% Map symbolic variables to their numerical values for substitutions.
sym_vars = [m1, m2, m3, g, I1xx, I1yy, I1zz, I2xx, I2yy, I2zz, I3xx, I3yy, I3zz];
num_vals = [m1_val, m2_val, m3_val, g_val, I1xx_val, I1yy_val, I1zz_val, I2xx_val, I2yy_val, I2zz_val, I3xx_val, I3yy_val, I3zz_val];

% Create numerically-defined symbolic expressions for T and U.
T_num_exact = subs(T_kinetic_symbolic, sym_vars, num_vals);
U_num_exact = subs(U_potential_symbolic, sym_vars, num_vals);

%% 3. Trajectory Setup
% Simulation time parameters.
t_start = 0; t_end = 2.0;
t_vec = linspace(t_start, t_end, 100);

% Boundary conditions for the quintic polynomial trajectory.
q1_0 = -0.5; q1_f = 0.0; 
q2_0 = 0; q2_f = pi/4; 
q3_0 = 0; q3_f = pi/6;

get_quintic = @(q0, qf, tf, t) [ ...
    q0 + 10*(qf-q0)/tf^3*t^3 - 15*(qf-q0)/tf^4*t^4 + 6*(qf-q0)/tf^5*t^5, ... % Position
    15*(qf-q0)/tf^3*t^2 - 24*(qf-q0)/tf^4*t^3 + 15*(qf-q0)/tf^5*t^4, ...     % Velocity
    30*(qf-q0)/tf^3*t - 36*(qf-q0)/tf^4*t^2 + 15*(qf-q0)/tf^5*t^3 ...        % Acceleration
];

%% 4. Unified Window Setup (Dark Mode)
fig = figure('Name', 'PRR Robot: Live Energy Evolution', 'NumberTitle', 'off', 'Color', [0.15, 0.15, 0.15]);

% Subplot 1: Robot Animation (Sinistra - Righe 1 e 2)
ax_robot = subplot(4, 2, [1, 3]);
grid on; axis equal;
xlim([-1.5, 1.5]); ylim([-1.5, 1.5]); zlim([-0.2, 2.0]);
set(ax_robot, 'Color', [0.2, 0.2, 0.2], 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w', 'GridColor', [0.4, 0.4, 0.4]);
axis off;
title('Live Robot Pose', 'Color', 'w');

% Subplot 2: Traiettoria Giunto Prismatico q1 [m] (Sinistra - Riga 3)
ax_q_p = subplot(4, 2, 5);
h_q1 = plot(NaN, NaN, 'Color', [0, 1, 0], 'LineWidth', 1.5);
legend(h_q1, {'q_1 [m]'}, 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none', 'Location', 'northwest');
grid on; set(ax_q_p, 'Color', [0.2, 0.2, 0.2], 'XColor', 'w', 'YColor', 'w', 'GridColor', [0.4, 0.4, 0.4]);
ylabel('Position [m]', 'Color', 'w'); 
title('Prismatic Joint (q_1)', 'Color', 'w');
xlim([t_start, t_end]);
ylim([-0.6, 0.1]);
set(gca, 'FontSize', 8);

% Subplot 3: Traiettorie Giunti Rotoidali q2, q3 [rad] (Sinistra - Riga 4)
ax_q_r = subplot(4, 2, 7);
h_q2 = plot(NaN, NaN, 'Color', [1, 1, 0], 'LineWidth', 1.5); hold on;
h_q3 = plot(NaN, NaN, 'Color', [0, 1, 1], 'LineWidth', 1.5);
legend([h_q2, h_q3], {'q_2 [rad]', 'q_3 [rad]'}, 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none', 'Location', 'northwest');
grid on; set(ax_q_r, 'Color', [0.2, 0.2, 0.2], 'XColor', 'w', 'YColor', 'w', 'GridColor', [0.4, 0.4, 0.4]);
xlabel('Time [s]', 'Color', 'w'); 
ylabel('Angle [rad]', 'Color', 'w'); 
title('Revolute Joints (q_2, q_3)', 'Color', 'w');
xlim([t_start, t_end]);
ylim([-0.2, 1.0]);
set(gca, 'FontSize', 8);

% Subplot 4: Energia Cinetica (Destra Alto)
ax_t = subplot(3, 2, 2);
h_t = plot(NaN, NaN, 'Color', [0, 0.8, 1], 'LineWidth', 2);
grid on; set(ax_t, 'Color', [0.2, 0.2, 0.2], 'XColor', 'w', 'YColor', 'w', 'GridColor', [0.4, 0.4, 0.4]);
ylabel('T [J]', 'Color', 'w'); title('Kinetic Energy', 'Color', 'w');
xlim([t_start, t_end]);
ylim([0, 4]);

% Subplot 5: Energia Potenziale (Destra Centro)
ax_u = subplot(3, 2, 4);
h_u = plot(NaN, NaN, 'Color', [1, 0.4, 0.4], 'LineWidth', 2);
grid on; set(ax_u, 'Color', [0.2, 0.2, 0.2], 'XColor', 'w', 'YColor', 'w', 'GridColor', [0.4, 0.4, 0.4]);
ylabel('U [J]', 'Color', 'w'); title('Potential Energy', 'Color', 'w');
xlim([t_start, t_end]);
ylim([9.81, 26]);

% Subplot 6: Energia Totale (Destra Basso)
ax_e = subplot(3, 2, 6);
h_e = plot(NaN, NaN, 'Color', [0.9, 0.9, 0.9], 'LineWidth', 2);
grid on; set(ax_e, 'Color', [0.2, 0.2, 0.2], 'XColor', 'w', 'YColor', 'w', 'GridColor', [0.4, 0.4, 0.4]);
xlabel('Time [s]', 'Color', 'w'); ylabel('E [J]', 'Color', 'w'); title('Total Energy', 'Color', 'w');
xlim([t_start, t_end]);
ylim([9.81, 30]);

%% 5. Live Execution Loop
t_hist = zeros(1, length(t_vec));
q1_hist = zeros(1, length(t_vec));
q2_hist = zeros(1, length(t_vec));
q3_hist = zeros(1, length(t_vec));
T_hist = zeros(1, length(t_vec));
U_hist = zeros(1, length(t_vec));
E_hist = zeros(1, length(t_vec));

fprintf('\nCompiling symbolic expressions into numeric functions for faster execution...\n');
T_func = matlabFunction(T_num_exact, 'Vars', [q1; q2; q3; dq1; dq2; dq3]);
U_func = matlabFunction(U_num_exact, 'Vars', [q1; q2; q3]);
fprintf('Compilation complete! Starting synchronized animation...\n');

q_init = [q1_0; q2_0; q3_0];
show(robot, q_init, 'Parent', ax_robot, 'PreservePlot', false);
h_title = title(ax_robot, 'Live Robot Pose - Time: 0.00 s', 'Color', 'w');

sync_timer = tic;

for k = 1:length(t_vec)
    t = t_vec(k);

    r1 = get_quintic(q1_0, q1_f, t_end, t);
    r2 = get_quintic(q2_0, q2_f, t_end, t);
    r3 = get_quintic(q3_0, q3_f, t_end, t);
    
    q_curr = [r1(1); r2(1); r3(1)];
    dq_curr = [r1(2); r2(2); r3(2)];

    show(robot, q_curr, 'Parent', ax_robot, 'FastUpdate', true, 'PreservePlot', false);
    set(h_title, 'String', sprintf('Live Robot Pose - Time: %.2f s', t));

    T_val = T_func(q_curr(1), q_curr(2), q_curr(3), dq_curr(1), dq_curr(2), dq_curr(3));
    U_val = U_func(q_curr(1), q_curr(2), q_curr(3));
    E_val = T_val + U_val;

    t_hist(k) = t;
    q1_hist(k) = q_curr(1);
    q2_hist(k) = q_curr(2);
    q3_hist(k) = q_curr(3);
    T_hist(k) = T_val;
    U_hist(k) = U_val;
    E_hist(k) = E_val;

    set(h_q1, 'XData', t_hist(1:k), 'YData', q1_hist(1:k));
    set(h_q2, 'XData', t_hist(1:k), 'YData', q2_hist(1:k));
    set(h_q3, 'XData', t_hist(1:k), 'YData', q3_hist(1:k));
    set(h_t, 'XData', t_hist(1:k), 'YData', T_hist(1:k));
    set(h_u, 'XData', t_hist(1:k), 'YData', U_hist(1:k));
    set(h_e, 'XData', t_hist(1:k), 'YData', E_hist(1:k));

    drawnow limitrate;

    elapsed_time = toc(sync_timer);
    time_to_wait = t - elapsed_time;
    if time_to_wait > 0
        pause(time_to_wait);
    end
end

%% 6. Console Output
fprintf('\n------------------------------------------------------------\n');
fprintf('FINAL DYNAMIC ANALYSIS SUMMARY\n');
fprintf('------------------------------------------------------------\n');
fprintf('Sample Time t = %.2f s (End of trajectory)\n', t_end);
fprintf('Final Configuration q_final: [%.4f, %.4f, %.4f]^T\n', q_curr);
fprintf('Final Velocities dq_final: [%.4f, %.4f, %.4f]^T\n\n', dq_curr);

fprintf('Peak Values during motion:\n');
fprintf('Max Kinetic Energy T:   %.6f J\n', max(T_hist));
fprintf('Max Potential Energy U: %.6f J\n', max(U_hist));
fprintf('Max Total Energy E:     %.6f J\n', max(E_hist));

fprintf('\nEvaluated Inertia Matrix B(q_final):\n');
B_eval_final_num = double(subs(subs(B_symbolic, sym_vars, num_vals), [q1, q2, q3], q_curr'));
disp(B_eval_final_num);

fprintf('Final Energy State:\n');
fprintf('T: %.6f J | U: %.6f J | E: %.6f J\n', T_val, U_val, E_val);
fprintf('============================================================\n');