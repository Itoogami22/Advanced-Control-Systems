%% ========================================================================
%  ADVANCED CONTROL SYSTEMS - ASSIGNMENT 3 (MATLAB DESKTOP VERSION)
%  Euler-Lagrange Dynamic Equations of Motion for 3-DOF PRR Manipulator
% ========================================================================

close all; clear; clc;

%% 1. Robot Import and Interactive 3D Visualization
robot = importrobot('00001_PRR_zzx FILIPPO SEGALLA.urdf');
robot.DataFormat = 'column';
showdetails(robot)

fig = figure('Name', 'PRR Robot Visualization - Assignment 3', 'NumberTitle', 'off');
config = homeConfiguration(robot);
show(robot, config);

rotate3d on;              
cameratoolbar('Show');    
grid on;                  
axis equal;              

xlim([-1.5, 1.5]);        
ylim([-1.5, 1.5]);        
zlim([-0.2, 2.0]);        

%% 2. Symbolic Variables Definition
num_joints = 3;           

% Symbolic variables for DH parameters
syms a alpha d theta real

% Generalized coordinates, velocities, and accelerations
syms q1 q2 q3 real                                  
syms dq1 dq2 dq3 real                               
syms ddq1 ddq2 ddq3 real                            

% Symbolic dynamic parameters (letters)
syms m1 m2 m3 real                                  % Link masses (kg)
syms g real                                         % Gravitational acceleration (m/s^2)
syms I1xx I1yy I1zz real                            % Link 1 inertia moments
syms I2xx I2yy I2zz real                            % Link 2 inertia moments
syms I3xx I3yy I3zz real                            % Link 3 inertia moments

q   = [q1; q2; q3];
dq  = [dq1; dq2; dq3];
ddq = [ddq1; ddq2; ddq3];
m   = [m1, m2, m3];

% Gravity acceleration vector expressed in Base Frame Sigma_0
g0 = [0; 0; -g];                                    

% Local inertia tensors relative to CoMs in local link frames
I_local = cell(num_joints, 1);
I_local{1} = diag([I1xx, I1yy, I1zz]);
I_local{2} = diag([I2xx, I2yy, I2zz]);
I_local{3} = diag([I3xx, I3yy, I3zz]);

%% 3. Forward Kinematics and DH Transformation Matrices (Exact Rationals)
A = [cos(theta), -sin(theta)*cos(alpha),  sin(theta)*sin(alpha), a*cos(theta); ...
     sin(theta),  cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta); ...
     0,           sin(alpha),             cos(alpha),            d; ...
     0,           0,                      0,                     1];

TB0 = [ 0, -1,  0,  0;
        1,  0,  0,  0;
        0,  0,  1,  sym('65/100');
        0,  0,  0,  1 ];

% Direct use of exact rational fractions to prevent floating point noise
DH_table = [ 0,               0,             q1,             0;
             0,               sym(pi)/2,     sym('4/10'),    q2;
             sym('4/10'),     0,             0,              q3 ];

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

%% 4. Center of Mass (CoM) Positions & Partial Jacobians
p_l_local = {[0; 0; -sym('1/4'); 1], [0; 0; -sym('1/5'); 1], [-sym('1/5'); 0; 0; 1]};
p_l = cell(num_joints, 1);

for i = 1:num_joints
    p_temp = simplify(T{i} * p_l_local{i});
    p_l{i} = p_temp(1:3);
end

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

%% 5. Computation of Steiner Inertia Tensors & Generalized Inertia Matrix B(q)
I_steiner = cell(num_joints, 1);
for i = 1:num_joints
    if joint_types(i) == 'R'
        c_i = p_l_local{i}(1:3); % CoM offset relative to joint frame
        % Steiner's Parallel Axis Theorem: S = m * ( (c'*c)*I - c*c' )
        S_matrix = m(i) * ((c_i' * c_i) * eye(3) - (c_i * c_i'));
        I_steiner{i} = simplify(I_local{i} + S_matrix);
    else
        I_steiner{i} = I_local{i}; % Prismatic link maintains local inertia
    end
end

B_symbolic = sym(zeros(num_joints, num_joints));
for i = 1:num_joints
    B_symbolic = B_symbolic + m(i) * (JP_l{i}' * JP_l{i}) + JO_l{i}' * R{i} * I_steiner{i} * R{i}' * JO_l{i};
end

B_symbolic = simplify(B_symbolic);

%% 6. Explicit Christoffel Symbols c_ijk(q) and Coriolis Matrix C(q, dq)
c_ijk_sym = sym(zeros(num_joints, num_joints, num_joints));
C_symbolic = sym(zeros(num_joints, num_joints));

for i = 1:num_joints
    for j = 1:num_joints
        c_ij = sym(0);
        for k = 1:num_joints
            c_ijk = 0.5 * ( diff(B_symbolic(i,j), q(k)) + ...
                            diff(B_symbolic(i,k), q(j)) - ...
                            diff(B_symbolic(j,k), q(i)) );
            c_ijk_sym(i,j,k) = simplify(c_ijk);
            c_ij = c_ij + c_ijk_sym(i,j,k) * dq(k);
        end
        C_symbolic(i,j) = simplify(c_ij);
    end
end

%% 7. Gravity Vector g(q)
U_potential_symbolic = sym(0);
for i = 1:num_joints
    U_potential_symbolic = U_potential_symbolic - m(i) * (g0' * p_l{i});
end

g_symbolic = sym(zeros(num_joints, 1));
for i = 1:num_joints
    g_symbolic(i) = simplify(diff(U_potential_symbolic, q(i)));
end

%% 8. Purely Symbolic Equations of Motion Tau(q, dq, ddq)
tau_symbolic = simplify(B_symbolic * ddq + C_symbolic * dq + g_symbolic);

tau1_sym = tau_symbolic(1);
tau2_sym = tau_symbolic(2);
tau3_sym = tau_symbolic(3);

%% 9. Skew-Symmetry Verification of N(q, dq) = B_dot - 2*C
B_dot_symbolic = sym(zeros(num_joints, num_joints));
for k = 1:num_joints
    B_dot_symbolic = B_dot_symbolic + diff(B_symbolic, q(k)) * dq(k);
end

N_symbolic = simplify(B_dot_symbolic - 2 * C_symbolic);
is_skew_symmetric = simplify(N_symbolic + N_symbolic'); % Should be zeros(3,3)

%% 10. Numerical Substitutions (Exact Symbolic Rationals)
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

c_ijk_num_exact = subs(c_ijk_sym, sym_vars, num_vals);
B_num_exact     = subs(B_symbolic, sym_vars, num_vals);
C_num_exact     = subs(C_symbolic, sym_vars, num_vals);
g_num_exact     = subs(g_symbolic, sym_vars, num_vals);
tau_num_exact   = subs(tau_symbolic, sym_vars, num_vals);
N_num_exact     = subs(N_symbolic, sym_vars, num_vals);

% Format clean numerical outputs with 4 significant digits (vpa)
B_numeric   = vpa(expand(B_num_exact), 4);
C_numeric   = vpa(expand(C_num_exact), 4);
g_numeric   = vpa(expand(g_num_exact), 4);
tau_numeric = vpa(expand(tau_num_exact), 4);
N_numeric   = vpa(expand(N_num_exact), 4);

tau1_num = tau_numeric(1);
tau2_num = tau_numeric(2);
tau3_num = tau_numeric(3);

%% 11. Display Console Output
fprintf('------------------------------------------------------------\n');
fprintf('PURELY SYMBOLIC FORMULAS (LETTERS & FRACTIONS)\n');
fprintf('------------------------------------------------------------\n');

fprintf('\n--- Inertia Matrix B(q) (3x3) ---\n');
disp(B_symbolic);

fprintf('\n------------------------------------------------------------\n');
fprintf('DETAILED SYMBOLIC CHRISTOFFEL SYMBOLS c_ijk(q)\n');
fprintf('------------------------------------------------------------\n');
count_nonzero_sym = 0;
for i = 1:num_joints
    for j = 1:num_joints
        for k = 1:num_joints
            val = c_ijk_sym(i,j,k);
            if ~isequal(val, sym(0))
                count_nonzero_sym = count_nonzero_sym + 1;
                fprintf('  c_%d%d%d = %s\n', i, j, k, char(val));
            end
        end
    end
end
fprintf('Total non-zero symbolic Christoffel terms: %d / 27\n', count_nonzero_sym);

fprintf('\n--- Coriolis Matrix C(q, dq) (3x3) ---\n');
disp(C_symbolic);

fprintf('\n--- Gravity Vector g(q) (3x1) ---\n');
disp(g_symbolic);

fprintf('\n------------------------------------------------------------\n');
fprintf('EXPLICIT JOINT TORQUE EQUATIONS Tau_i(q, dq, ddq) (PURE SYMBOLIC)\n');
fprintf('------------------------------------------------------------\n');
fprintf('  tau_1 = %s\n\n', char(tau1_sym));
fprintf('  tau_2 = %s\n\n', char(tau2_sym));
fprintf('  tau_3 = %s\n\n', char(tau3_sym));

fprintf('\n--- Skew-Symmetry Verification (N = B_dot - 2*C) ---\n');
disp(N_symbolic);
fprintf('Verification of N + N^T == 0:\n');
disp(is_skew_symmetric);

fprintf('\n------------------------------------------------------------\n');
fprintf('EVALUATED NUMERICAL VALUES\n');
fprintf('------------------------------------------------------------\n');

fprintf('\n------------------------------------------------------------\n');
fprintf('EVALUATED NUMERICAL CHRISTOFFEL SYMBOLS c_ijk(q)\n');
fprintf('------------------------------------------------------------\n');
count_nonzero_num = 0;
for i = 1:num_joints
    for j = 1:num_joints
        for k = 1:num_joints
            val_exact = c_ijk_num_exact(i,j,k);
            if ~isequal(val_exact, sym(0))
                count_nonzero_num = count_nonzero_num + 1;
                val_vpa = vpa(expand(val_exact), 4);
                fprintf('  c_%d%d%d = %s\n', i, j, k, char(val_vpa));
            end
        end
    end
end
if count_nonzero_num == 0
    fprintf('  All evaluated Christoffel symbols are zero.\n');
else
    fprintf('  Total non-zero numerical Christoffel terms: %d / 27\n', count_nonzero_num);
end

fprintf('\n------------------------------------------------------------\n');
fprintf('EXPLICIT JOINT TORQUE EQUATIONS Tau_i (NUMERICAL)\n');
fprintf('------------------------------------------------------------\n');
fprintf('  tau_1 = %s\n\n', char(tau1_num));
fprintf('  tau_2 = %s\n\n', char(tau2_num));
fprintf('  tau_3 = %s\n\n', char(tau3_num));

fprintf('\n--- Evaluated Inertia Matrix B(q) (3x3) ---\n');
disp(B_numeric);

fprintf('\n--- Evaluated Coriolis Matrix C(q, dq) (3x3) ---\n');
disp(C_numeric);

fprintf('\n--- Evaluated Gravity Vector g(q) (3x1) ---\n');
disp(g_numeric);

fprintf('\n--- Evaluated Vector Tau(q, dq, ddq) ---\n');
disp(tau_numeric);

fprintf('\n--- Evaluated Matrix N(q, dq) = B_dot - 2*C ---\n');
disp(N_numeric);
fprintf('============================================================\n\n');