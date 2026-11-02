%% ========================================================================
%  ADVANCED CONTROL SYSTEMS - SECOND ASSIGNMENT (DYNAMIC ANALYSIS VERSION)
% ========================================================================
% This script calculates the kinetic and potential energies, as well as the
% generalized inertia matrix B(q), for a 3-DOF PRR manipulator.
% Instead of a single sample point, it analyzes the evolution of energy
% over the entire trajectory and plots the results.
% ========================================================================

close all; clear; clc;

%% 1. Robot Import and Setup
robot = importrobot('00001_PRR_zzx FILIPPO SEGALLA.urdf');
robot.DataFormat = 'column';

%% 2. Symbolic Variables Definition
num_joints = 3;
syms a alpha d theta real
syms q1 q2 q3 real
syms dq1 dq2 dq3 real
syms m1 m2 m3 real
syms g real
syms I1xx I1yy I1zz real
syms I2xx I2yy I2zz real
syms I3xx I3yy I3zz real

q  = [q1; q2; q3];
dq = [dq1; dq2; dq3];
m  = [m1, m2, m3];
g0 = [0; 0; -g];

I_local = cell(num_joints, 1);
I_local{1} = diag([I1xx, I1yy, I1zz]);
I_local{2} = diag([I2xx, I2yy, I2zz]);
I_local{3} = diag([I3xx, I3yy, I3zz]);

%% 3. Direct Kinematics and DH Transformation Matrices
A = [cos(theta), -sin(theta)*cos(alpha),  sin(theta)*sin(alpha), a*cos(theta); ...
     sin(theta),  cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta); ...
     0,           sin(alpha),             cos(alpha),            d; ...
     0,           0,                      0,                     1];

TB0 = [ 0, -1,  0,  0;
        1,  0,  0,  0;
        0,  0,  1,  0.65;
        0,  0,  0,  1 ];

DH_table = [ 0,    0,     q1,   0;       % Link 1: Prismatic
             0,    pi/2,  0.4,  q2;     % Link 2: Revolute
             0.4,  0,     0,    q3 ];    % Link 3: Revolute

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

%% 4. Center of Mass (CoM) Positions in Base Frame
p_l_local = {[0; 0; -0.25; 1], [0; 0; -0.20; 1], [-0.20; 0; 0; 1]};
p_l = cell(num_joints, 1);

for i = 1:num_joints
    p_temp = simplify(T{i} * p_l_local{i});
    p_l{i} = p_temp(1:3);
end

%% 5. Partial Jacobians for Link CoMs
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
        else
            JP_i(:, j) = cross(z_prev, (p_l{i} - p_prev));
            JO_i(:, j) = z_prev;
        end
    end
    JP_l{i} = simplify(JP_i);
    JO_l{i} = simplify(JO_i);
end

%% 6. Energy Formulas (Symbolic)
B_symbolic = sym(zeros(num_joints, num_joints));
for i = 1:num_joints
    B_symbolic = B_symbolic + m(i) * (JP_l{i}' * JP_l{i}) + JO_l{i}' * R{i} * I_local{i} * R{i}' * JO_l{i};
end
B_symbolic = simplify(B_symbolic);
T_kinetic_symbolic = simplify(0.5 * dq' * B_symbolic * dq);
U_potential_symbolic = sym(0);
for i = 1:num_joints
    U_potential_symbolic = U_potential_symbolic - m(i) * (g0' * p_l{i});
end
U_potential_symbolic = simplify(U_potential_symbolic);

%% 7. Numerical Substitution (Physical Constants)
m1_val = sym(1); m2_val = sym(1); m3_val = sym(1);
g_val  = sym('981/100');
I1xx_val = sym('1/12') * m1_val * (sym('36/1000')^2 + sym('1/2')^2);
I1yy_val = sym('1/12') * m1_val * (sym('36/1000')^2 + sym('1/2')^2);
I1zz_val = sym('1/12') * m1_val * (sym('36/1000')^2 + sym('36/1000')^2);
I2xx_val = sym('1/12') * m2_val * (3*sym('19/1000')^2 + sym('4/10')^2);
I2yy_val = sym('1/12') * m2_val * (3*sym('19/1000')^2 + sym('4/10')^2);
I2zz_val = sym('1/2')  * m2_val * (sym('19/1000')^2);
I3xx_val = sym('1/2')  * m3_val * (sym('19/1000')^2);
I3yy_val = sym('1/12') * m3_val * (3*sym('19/1000')^2 + sym('4/10')^2);
I3zz_val = sym('1/12') * m3_val * (3*sym('19/1000')^2 + sym('4/10')^2);

sym_vars = [m1, m2, m3, g, I1xx, I1yy, I1zz, I2xx, I2yy, I2zz, I3xx, I3yy, I3zz];
num_vals = [m1_val, m2_val, m3_val, g_val, I1xx_val, I1yy_val, I1zz_val, I2xx_val, I2yy_val, I2zz_val, I3xx_val, I3yy_val, I3zz_val];

T_num_exact = subs(T_kinetic_symbolic, sym_vars, num_vals);
U_num_exact = subs(U_potential_symbolic, sym_vars, num_vals);

%% 8. Trajectory Setup
t_start = 0;
t_end = 2;
t_vec = linspace(t_start, t_end, 100); % 100 points for smooth plots

q1_0 = -0.5; q1_f = 0.0;
q2_0 = 0; q2_f = pi/4;
q3_0 = 0; q3_f = pi/6;

get_quintic = @(q0, qf, tf, t) [ ...
    q0 + 10*(qf-q0)/tf^3*t^3 - 15*(qf-q0)/tf^4*t^4 + 6*(qf-q0)/tf^5*t^5, ...
    15*(qf-q0)/tf^3*t^2 - 24*(qf-q0)/tf^4*t^3 + 15*(qf-q0)/tf^5*t^4 ...
];

%% 9. Dynamic Evaluation over the Trajectory
T_history = zeros(size(t_vec));
U_history = zeros(size(t_vec));
E_history = zeros(size(t_vec));

fprintf('Calculating energy evolution... please wait.\n');

for k = 1:length(t_vec)
    t = t_vec(k);

    % Current joint state
    r1 = get_quintic(q1_0, q1_f, t_end, t);
    r2 = get_quintic(q2_0, q2_f, t_end, t);
    r3 = get_quintic(q3_0, q3_f, t_end, t);

    q_curr = [r1(1); r2(1); r3(1)];
    dq_curr = [r1(2); r2(2); r3(2)];

    % Numerical evaluation
    T_history(k) = double(subs(T_num_exact, [q1; q2; q3; dq1; dq2; dq3], [q_curr; dq_curr]));
    U_history(k) = double(subs(U_num_exact, [q1; q2; q3], q_curr));
    E_history(k) = T_history(k) + U_history(k);
end

%% 10. Plotting the Results
figure('Name', 'Dynamic Energy Analysis of PRR Robot', 'NumberTitle', 'off', 'Color', [0.15, 0.15, 0.15]);

subplot(3,1,1);
plot(t_vec, T_history, 'Color', [0, 0.8, 1], 'LineWidth', 2);
grid on;
set(gca, 'Color', [0.2, 0.2, 0.2], 'XColor', 'w', 'YColor', 'w', 'GridColor', [0.4, 0.4, 0.4]);
ylabel('Kinetic Energy T [J]', 'Color', 'w');
title('Energy Evolution along the Quintic Trajectory', 'Color', 'w');

subplot(3,1,2);
plot(t_vec, U_history, 'Color', [1, 0.4, 0.4], 'LineWidth', 2);
grid on;
set(gca, 'Color', [0.2, 0.2, 0.2], 'XColor', 'w', 'YColor', 'w', 'GridColor', [0.4, 0.4, 0.4]);
ylabel('Potential Energy U [J]', 'Color', 'w');

subplot(3,1,3);
plot(t_vec, E_history, 'Color', [0.9, 0.9, 0.9], 'LineWidth', 2);
grid on;
set(gca, 'Color', [0.2, 0.2, 0.2], 'XColor', 'w', 'YColor', 'w', 'GridColor', [0.4, 0.4, 0.4]);
xlabel('Time [s]', 'Color', 'w'); ylabel('Total Energy E [J]', 'Color', 'w');

% Final Summary
fprintf('\n------------------------------------------------------------\n');
fprintf('DYNAMIC ANALYSIS COMPLETE\n');
fprintf('------------------------------------------------------------\n');
fprintf('Max Kinetic Energy:   %.6f J\n', max(T_history));
fprintf('Max Potential Energy: %.6f J\n', max(U_history));
fprintf('Total Energy Range:   [%.4f, %.4f] J\n', min(E_history), max(E_history));
fprintf('============================================================\n');
