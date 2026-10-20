%%  Kinetic and Potential Energy
close all; clear; clc;

%% Robot Import and Interactive 3D Visualization
% Import the robot model from URDF description file
robot = importrobot('00001_PRR_zzx FILIPPO SEGALLA.urdf');
showdetails(robot)

% Create figure window and display robot in home configuration
fig = figure;
config = homeConfiguration(robot);
show(robot, config);

% Enable interactive 3D controls (rotation, pan, zoom)
rotate3d on;              % Dynamic 3D rotation with mouse left-click drag
cameratoolbar('Show');    % Display top camera controls toolbar
grid on;                  % Enable grid lines for spatial reference
axis equal;               % Ensure 1:1 scale for all spatial axes

% Define extended workspace limits to center the robot visual grid
xlim([-1.5, 1.5]);        % Expanded X-axis limits
ylim([-1.5, 1.5]);        % Expanded Y-axis limits
zlim([-0.2, 2.0]);        % Expanded Z-axis limits

%% Symbolic Variables Definition (Generic Letters)
num_joints = 3;           % Number of degrees of freedom (DOF)

% Symbolic variables for DH parameters
syms a alpha d theta real

% Symbolic generalized coordinates and velocities
syms q1 q2 q3 real                                  
syms dq1 dq2 dq3 real                               

% Symbolic dynamic parameters (expressed as letters)
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
% Standard Denavit-Hartenberg homogeneous transformation matrix
A = [cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha) a*cos(theta); ...
     sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta); ...
     0           sin(alpha)             cos(alpha)            d; ...
     0           0                      0                     1];

% Fixed Base Frame transformation matrix TB0 (World -> Frame 0)
TB0 = [ 0 -1  0  0;
        1  0  0  0;
        0  0  1  0.65;
        0  0  0  1 ];

% Clean DH parameters table [a, alpha, d, theta]
DH_table = [ 0    0     q1   0;
             0    pi/2  0.4  q2;
             0.4  0     0    q3 ];

joint_types = ['P'; 'R'; 'R'];                     % Joint classification

% Initialize storage cell arrays
T = cell(num_joints, 1);
R = cell(num_joints, 1);
T_curr = sym(TB0);

% Compute cumulative homogeneous transformations T_0^i and rotation matrices R_i
for i = 1:num_joints
    A_i = subs(A, [a alpha d theta], DH_table(i,:));
    T_curr = simplify(T_curr * A_i);
    T{i} = T_curr;                                % Transformation T_0^i
    R{i} = T{i}(1:3, 1:3);                         % Rotation matrix R_i
end

%% Center of Mass (CoM) Positions in Base Frame Sigma_0
% Local CoM positions in link frames (homogeneous coordinates)
p_l_local = {[0; 0; -0.25; 1], [0; 0; -0.20; 1], [-0.20; 0; 0; 1]};
p_l = cell(num_joints, 1);

% Compute CoM positions in Base Frame Sigma_0: p_li = T_0^i * p_li^i
for i = 1:num_joints
    p_temp = simplify(T{i} * p_l_local{i});
    p_l{i} = p_temp(1:3);
end

%% Partial Jacobians for Link CoMs (JP_li, JO_li)
JP_l = cell(num_joints, 1);
JO_l = cell(num_joints, 1);

% Store joint origins and Z-axes in Base Frame Sigma_0
z_axes = cell(num_joints+1, 1);
p_origins = cell(num_joints+1, 1);

z_axes{1} = sym(TB0(1:3, 3));
p_origins{1} = sym(TB0(1:3, 4));

for i = 1:num_joints
    z_axes{i+1} = T{i}(1:3, 3);
    p_origins{i+1} = T{i}(1:3, 4);
end

% Construct partial Jacobians via column-by-column evaluation
for i = 1:num_joints
    JP_i = sym(zeros(3, num_joints));
    JO_i = sym(zeros(3, num_joints));
    
    for j = 1:i
        z_prev = z_axes{j};
        p_prev = p_origins{j};
        
        if joint_types(j) == 'P'
            JP_i(:, j) = z_prev;                  % Prismatic joint linear term
            JO_i(:, j) = [0; 0; 0];               % Prismatic joint angular term
        else % Revolute 'R'
            JP_i(:, j) = cross(z_prev, (p_l{i} - p_prev)); % Revolute joint linear term
            JO_i(:, j) = z_prev;                  % Revolute joint angular term
        end
    end
    
    JP_l{i} = simplify(JP_i);
    JO_l{i} = simplify(JO_i);
end

%% Inertia Matrix B(q) Computation (Symbolic Letters)
B_symbolic = sym(zeros(num_joints, num_joints));

for i = 1:num_joints
    B_symbolic = B_symbolic + m(i) * (JP_l{i}' * JP_l{i}) + JO_l{i}' * R{i} * I_local{i} * R{i}' * JO_l{i};
end

B_symbolic = simplify(B_symbolic);

%% Total Kinetic Energy T(q, dq) and Potential Energy U(q) (Symbolic Letters)
% Kinetic Energy: T(q, dq) = 0.5 * dq' * B(q) * dq
T_kinetic_symbolic = simplify(0.5 * dq' * B_symbolic * dq);

% Potential Energy: U(q) = - sum( m_i * g0' * p_li )
U_potential_symbolic = sym(0);
for i = 1:num_joints
    U_potential_symbolic = U_potential_symbolic - m(i) * (g0' * p_l{i});
end
U_potential_symbolic = simplify(U_potential_symbolic);

%% Numerical Value Substitution
% Numerical values for link masses (kg)
m1_val = 1.0; m2_val = 1.0; m3_val = 1.0;
g_val  = 9.81;

% Numerical inertia values derived from link geometries
I1xx_val = (1/12) * m1_val * (0.036^2 + 0.5^2);
I1yy_val = (1/12) * m1_val * (0.036^2 + 0.5^2);
I1zz_val = (1/12) * m1_val * (0.036^2 + 0.036^2);

I2xx_val = (1/12) * m2_val * (3*0.019^2 + 0.4^2);
I2yy_val = (1/12) * m2_val * (3*0.019^2 + 0.4^2);
I2zz_val = (1/2)  * m2_val * (0.019^2);

I3xx_val = (1/2)  * m3_val * (0.019^2);
I3yy_val = (1/12) * m3_val * (3*0.019^2 + 0.4^2);
I3zz_val = (1/12) * m3_val * (3*0.019^2 + 0.4^2);

% Vectors of symbolic variables and corresponding numerical values
sym_vars = [m1, m2, m3, g, ...
            I1xx, I1yy, I1zz, ...
            I2xx, I2yy, I2zz, ...
            I3xx, I3yy, I3zz];

num_vals = [m1_val, m2_val, m3_val, g_val, ...
            I1xx_val, I1yy_val, I1zz_val, ...
            I2xx_val, I2yy_val, I2zz_val, ...
            I3xx_val, I3yy_val, I3zz_val];

% Evaluate numerical matrices and energies
B_numeric           = vpa(subs(B_symbolic, sym_vars, num_vals), 4);
T_kinetic_numeric   = vpa(subs(T_kinetic_symbolic, sym_vars, num_vals), 4);
U_potential_numeric = vpa(subs(U_potential_symbolic, sym_vars, num_vals), 4);

%% Display Results (1. Purely Symbolic, 2. Numerical Decimal)
fprintf('------------------------------------------------------------\n');
fprintf('PURELY SYMBOLIC FORMULAS (LETTERS)\n');
fprintf('------------------------------------------------------------\n');

fprintf('\n--- CoM Positions p_l in Base Frame Sigma_0 ---\n');
for i = 1:num_joints
    fprintf('p_l%d =\n', i);
    disp(p_l{i});
end

fprintf('\n--- Symbolic Inertia Matrix B(q) (3x3) ---\n');
disp(B_symbolic);

fprintf('\n--- Symbolic Kinetic Energy T(q, dq) ---\n');
disp(T_kinetic_symbolic);

fprintf('\n--- Symbolic Potential Energy U(q) ---\n');
disp(U_potential_symbolic);

fprintf('\n------------------------------------------------------------\n');
fprintf('EVALUATED NUMERICAL VALUES (DECIMALS)\n');
fprintf('------------------------------------------------------------\n');

fprintf('\n--- Evaluated Inertia Matrix B(q) (3x3) ---\n');
disp(B_numeric);

fprintf('\n--- Evaluated Kinetic Energy T(q, dq) ---\n');
disp(T_kinetic_numeric);

fprintf('\n--- Evaluated Potential Energy U(q) (with g = 9.81) ---\n');
disp(U_potential_numeric);
fprintf('============================================================\n\n');