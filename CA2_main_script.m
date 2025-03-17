clear all;
clc;

% Pre-Processing

% Load mesh data
load Mesh_data_2024_coarse.mat
thick = 1;  % Out-of-plane thickness of element

% Input parameters
T_air = 20;
T_water = 30;
kt = 1.0;    % Tile conductivity
kc = 2.0;    % Concrete conductivity
alpha_wc = 1000; % Transfer coefficient water to concrete
alpha_ta = 10;   % Transfer coefficient tiles to air
c_pt = 800;  % Heat capacity of tiles
c_pc = 1000; % Heat capacity of concrete
rho_t = 1900; % Density of tiles
rho_c = 2300; % Density of concrete

% Conductivity matrices
D_concrete = [2 0; 0 2]; % Conductivity matrix for concrete
D_tile = [1 0; 0 1]; % Conductivity matrix for tile

NoDofs = size(Dof, 1); % Number of global DOFs per node

% ============================================================
% Solving the Problem
%============================================================

% ------------------------------------------------------------
% Calculation of the Convective Stiffness Matrix
% ------------------------------------------------------------

Kc = zeros(NoDofs);  % Sparse conductivity matrix for surface 1 (Air)
fc = zeros(NoDofs, 1);  % Sparse load vector for surface 1 (Air)

for i = 1:size(boundaryEdof, 1)
    
    % Identify the type of boundary condition (Air or Water convection)
    if boundaryMaterial(i, 2) == 1  % Air convection (Γ₁)
        alpha = alpha_ta;
        Tamb = T_air;
    elseif boundaryMaterial(i, 2) == 2  % Water convection (Γ₅)
        alpha = alpha_wc;
        Tamb = T_water;
    end

    % Compute convective element contributions
    [Kc_e, fc_e] = convecte(boundaryEx(i, :), boundaryEy(i, :), alpha, thick, Tamb);
    
    % Assemble into the air convection system
    edof = boundaryEdof(i, 2:3);
    Kc(edof, edof) = Kc(edof, edof) + Kc_e;
    fc(edof, 1) = fc(edof, 1) + fc_e;
    
end

%------------------------------------------------------------
% Calculation of the Conductive Stiffness Matrix
% ------------------------------------------------------------

K_stiffness = zeros(NoDofs);  % Sparse conductivity matrix for Concrete
fl = zeros(NoDofs, 1);  % Sparse load vector for Concrete

for i = 1:size(Edof, 1)
    
    % Extract element coordinates
    ex = Ex(i, :);
    ey = Ey(i, :);
    ep = thick;
    eq = 0;

    % Identify the type of material (Tile or Concrete)
    if matrlIndex(i, 1) == 1  % Tiles
        D = D_tile;
    elseif matrlIndex(i, 1) == 2  % Concrete
        D = D_concrete;
    end

    % Compute element matrix and load vector
    [K_e, f_e] = flw2te(ex, ey, ep, D, eq);
    
    % Assemble into the tile conduction system
    edof = Edof;
    [K_stiffness, fl] = assem(edof(i, :), K_stiffness, K_e, fl, f_e);

end

% Combine contributions into the final system
K = K_stiffness + Kc;
f = fl + fc;

% Solve system of equations
[a, Q] = solveq(K, f);


% ------------------------------------------------------------
%Plotting temperature distribution
% ------------------------------------------------------------

% Convert node positions to mm
x = Coord(:,1) * 1000;  
y = Coord(:,2) * 1000;  

% Plot temperature distribution
figure(1);
trisurf(Edof(:,2:4), x, y, a,'FaceColor', 'interp', 'EdgeColor', 'none');
colorbar;
colormap('jet');
xlabel('X Position (mm)');
ylabel('Y Position (mm)');
zlabel('Temperature (°C)');
title('Temperature Distribution Through the Concrete Slab (T_w = 30°C)');
view(2); % 2D view

% Overlay filtered mesh
hold on;
triplot(Edof(:,2:4), x, y, 'k'); 
hold off;



%------------------------------------------------------------
% Compute Heat Flow Vectors for Each Element
% ------------------------------------------------------------

nel = size(Edof, 1);  % Number of elements
es = zeros(nel, 2);   % Initialize heat flow matrix

for i = 1:nel
    % Get element nodes
    nodes = Edof(i, 2:4);
    x = Coord(nodes, 1);
    y = Coord(nodes, 2);
    
    % Compute B matrix 
    A = polyarea(x, y);  % Element area
    B = [y(2)-y(3), y(3)-y(1), y(1)-y(2);
         x(3)-x(2), x(1)-x(3), x(2)-x(1)] / (2 * A);
    
    % Identify material type (Tile or Concrete)
    if matrlIndex(i, 1) == 1  % Tiles
        D = D_tile;
    elseif matrlIndex(i, 1) == 2  % Concrete
        D = D_concrete;
    end 
    
    % Compute heat flow vector (q = -D * B * a)
    a_e = a(nodes);
    q = -D * B * a_e;
    
    % Store heat flow values
    es(i, :) = q';
end

% ------------------------------------------------------------
% Plot Heat Flow Arrows using elflux2
%------------------------------------------------------------

figure(2);
hold on;
eldraw2(Ex, Ey);  % Draw the mesh
elflux2(Ex, Ey, es, [1, 4]);  % Plot heat flow arrows (solid red)
hold off;
title('Element Heat Flow Vectors (Tw = 30°C)');
xlabel('X Position (m)');
ylabel('Y Position (m)');
axis equal;


%------------------------------------------------------------
%Floor surface temperature plot
%------------------------------------------------------------
for i = 1:7
    T(i)=a(boundaryEdof(i,2));
    x(i)=boundaryEx(i,1);
end 
figure(3)
plot(x,T)
grid on 
xlabel('X Position (m)');
ylabel('Temperature (°C)')
title('Temperature distribution across the boundary surface 1')

Q_e=zeros(1,7);
for i=1:7
    alpha=10;
    ex=boundaryEx(i,:);
    ey=boundaryEy(i,:);
    Le=sqrt((ey(2)-ey(1))^2+(ex(2)-ex(1))^2);
    Ti(i)=a(boundaryEdof(i,2));
    Tj(i)=a(boundaryEdof(i,3));
    Q_e(i)=alpha*Le*thick*((Ti(i)-Tj(i))/2-T_air);
end 
Q=sum(Q_e)