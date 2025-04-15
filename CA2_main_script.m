clear all;
clc;


%--------------------------------------------------
% Pre-Processing
%---------------------------------------------------


% Load mesh data
load Mesh_data_2024.mat
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
xlabel('X Position (mm)','FontSize', 24);
ylabel('Y Position (mm)','FontSize', 24);
zlabel('Temperature (°C)','FontSize', 24);
title('Temperature Distribution Through the Concrete Slab (T_w = 30°C)','FontSize', 24);
set(gca, 'FontSize', 20);
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
elflux2(Ex, Ey, es, [1, 4],0.000005);  % Plot heat flow arrows (solid red)
hold off;
title('Element Heat Flow Vectors (Tw = 30°C)');
xlabel('X Position (m)');
ylabel('Y Position (m)');
axis equal;


%------------------------------------------------------------
%Floor surface temperature plot
%------------------------------------------------------------
S1_edof=find((boundaryMaterial(:,2)==1));
boundary_row=size((S1_edof),1); %rows in the edof matrix which corresponds to surface 1

for i = 1:boundary_row
    T(i)=a(boundaryEdof(i,2));
    x(i)=boundaryEx(i,1);   
end 
T(end)=a(boundaryEdof(boundary_row,3));
x(end)=boundaryEx(boundary_row,2);
figure(3)
plot(x,T)
grid on 
xlabel('X Position (m)');
ylabel('Temperature (°C)')
title('Temperature distribution across the boundary surface 1')

for i=1:boundary_row
    alpha=10;
    ex=boundaryEx(i,:);
    ey=boundaryEy(i,:);
    Le=sqrt((ey(2)-ey(1))^2+(ex(2)-ex(1))^2);
    Ti(i)=a(boundaryEdof(i,2));
    Tj(i)=a(boundaryEdof(i,3));
    Q_e(i)=alpha*Le*thick*((Ti(i)+Tj(i))/2-T_air);
end 
Q=sum(Q_e);
disp(Q)




%-------------------------------------------------------------------------
%Comparision of abaqus results and matlab results
%------------------------------------------------------------------------

figure(4);
hold on; 

% Plot simulation results
plot(x, T, 'b-', 'LineWidth', 2); % Plot T vs. x in blue

% Abaqus data
x_abaqus_coarse = [0.0000, 0.0054, 0.0103, 0.0151, 0.0201, 0.0250, 0.0300, 0.0350, 0.0400, 0.0450, ...
            0.0500, 0.0550, 0.0600, 0.0650, 0.0700, 0.0750, 0.0799, 0.0849, 0.0897, 0.0946, 0.1000];

T_abaqus_coarse = [27.5556, 27.5056, 27.3905, 27.2315, 27.0510, 26.8676, 26.6923, 26.5309, 26.3854, 26.2561, ...
            26.1423, 26.0431, 25.9573, 25.8839, 25.8220, 25.7709, 25.7300, 25.6988, 25.6769, 25.6635, 25.6583];


x_abaqus_fine = [0, 0.003070363076403737, 0.005973363760858774, 0.00887497328221798, 0.011789732612669468, ...
                 0.014715086668729782, 0.017647285014390945, 0.020583538338541985, 0.02352208085358143, 0.026461869478225708, ...
            0.029402321204543114, 0.032343123108148575, 0.03528410196304321, 0.03822517395019531, 0.041166290640830994, 0.044107429683208466, 0.04704857990145683, 0.05293089151382446, 0.05587204545736313, 0.0588131844997406, 0.06175430491566658, 0.06469538062810898, 0.06763637065887451, 0.07057719677686691, 0.07351769506931305, 0.07645757496356964, 0.07939630001783371, 0.08233289420604706, 0.08526572585105896, 0.08819224685430527, 0.09110905230045319, 0.09401410073041916, 0.0969223827123642, 0.10000000149011612];

T_abaqus_fine = [27.549598693847656, 27.53305435180664, 27.490398406982422, 27.425003051757812, 27.341243743896484, 27.244388580322266, 27.139535903930664, 27.03101348876953, 26.922182083129883, 26.81546974182129, 26.712514877319336, 26.614349365234375, 26.521554946899414, 26.434391021728516, 26.352914810180664, 26.27704620361328, 26.206632614135742, 26.081350326538086, 26.026033401489258, 25.975303649902344, 25.9289493560791, 25.886775970458984, 25.848602294921875, 25.814266204833984, 25.78361701965332, 25.756528854370117, 25.732891082763672, 25.712608337402344, 25.69559669494629, 25.681785583496094, 25.671106338500977, 25.66348648071289, 25.658823013305664, 25.65707778930664];

x_abaqus_quad = [0, 0.004419696982949972, 0.009083370678126812, 0.014041165821254253, 0.0191547442227602, 0.024311354383826256, 0.02945711649954319, 0.0345756821334362, 0.03966720774769783, 0.044736891984939575, 0.04979023337364197, 0.05483165755867958, 0.05986436828970909, 0.06489064544439316, 0.06991210579872131, 0.07492997497320175, 0.07994525134563446, 0.08495887368917465, 0.08997182548046112, 0.0949852392077446, 0.10000000149011612];

T_abaqus_quad = [27.549333572387695, 27.5172176361084, 27.421316146850586, 27.2708683013916, 27.088258743286133, 26.897403717041016, 26.71435546875, 26.546619415283203, 26.39654541015625, 26.26404571533203, 26.14806365966797, 26.04728126525879, 25.96039390563965, 25.886234283447266, 25.82379913330078, 25.772258758544922, 25.730939865112305, 25.699317932128906, 25.67699432373047, 25.66370391845703, 25.659284591674805];


plot(x_abaqus_coarse, T_abaqus_coarse, 'ro-', 'MarkerFaceColor', 'c')
xlabel('X Position (m)', 'FontSize', 24);
ylabel('Temperature (°C)', 'FontSize', 24);
title('Temperature Distribution Comparison', 'FontSize', 24);
legend('Matlab Results', 'Abaqus coarse mesh','Location', 'Best');
set(gca, 'FontSize', 20);
grid on; 
hold off; 


figure(5)
hold on; 
plot(x_abaqus_coarse, T_abaqus_coarse, 'ro-', 'MarkerFaceColor', 'c')
plot(x_abaqus_fine, T_abaqus_fine, 'ro-', 'MarkerFaceColor', 'y')
xlabel('X Position (m)', 'FontSize', 24);
ylabel('Temperature (°C)', 'FontSize', 24);
title('Temperature Distribution Comparison', 'FontSize', 24);
legend('Abaqus coarse mesh', 'Abaqus fine mesh','Location', 'Best');
grid on; 
set(gca, 'FontSize', 20);
hold off;


figure(6)
hold on; 
plot(x_abaqus_coarse, T_abaqus_coarse, 'ro-', 'MarkerFaceColor', 'c')
plot(x_abaqus_fine, T_abaqus_fine, 'ro-', 'MarkerFaceColor', 'y')
plot(x_abaqus_quad, T_abaqus_quad, 'ro-', 'MarkerFaceColor', 'r')
xlabel('X Position (m)', 'FontSize', 24);
ylabel('Temperature (°C)', 'FontSize', 24);
title('Temperature Distribution Comparison', 'FontSize', 24);
legend('Abaqus coarse mesh', 'Abaqus fine mesh','Abaqus quad mesh','Location', 'Best');
set(gca, 'FontSize', 20);
grid on; 
hold off;









