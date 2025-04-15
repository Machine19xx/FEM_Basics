clear all;
clc;



Edof=[1 1 2;
      2 2 3;
      3 3 4;
      4 4 5;
      5 5 6;
      6 6 7;
      7 7 8;
      8 8 9;
      9 9 10;
      10 10 11;
      11 11 12;
      12 12 13;
      13 13 14;
      14 14 15];

Boundary_edof=[1 1;
               14 15];
Dof=[1;2;3;4;5;6;7;8;9;10;11;12;13;14;15];
T_in=21;
T_out=2.3;
bc=[1 T_out; 15 T_in];
h_inside = 8; % Convective coefficient inside (W/m2C)
h_outside = 25; % Convective coefficient outside (W/m2C)
Nodof=size(Dof,1);
K_stiffness=zeros(Nodof);
F=zeros(Nodof,1);
A=1;
Ex=10^-3*[0 20;
        20 40;
        40 60;
        60 80;
        80 100;
        100 120;
        120 140; 
        140 160;
        160 180;
        180 200;
        200 220;
        220 222;
        222 247;
        247 260];
k_cond=[14 0.034 5.50 0.034 0.20];
for i = 1:size(Edof,1)
    ex = Ex(i,:);
    
    % Assign correct conductivity values
    if i == 1
        k_conduction = k_cond(1);
    elseif i >= 2 && i <= 11
        k_conduction = k_cond(2);
    elseif i == 12
        k_conduction = k_cond(3);
    elseif i == 13
        k_conduction = k_cond(4);
    elseif i == 14
        k_conduction = k_cond(5);
    end
    
    % Element stiffness matrix
    Ke = (A * k_conduction / (ex(2) - ex(1))) * [1 -1; -1 1];
    edof = Edof(i, 2:3);
    K_stiffness(edof, edof) = K_stiffness(edof, edof) + Ke;
end 
K=K_stiffness;

free = [2 3 4 5 6 7 8 9 10 11 12 13 14];  % Free DOFs
const = [1 15];  % Constrained DOFs
% Applying Dirichlet BCs for Problem 1.1
K_dir = K;
F_dir = F;
T_dir= zeros(Nodof,1);
T_dir(const) = bc(:,2);  % Apply known temperatures
T_dir(free) = K_dir(free, free) \ (F_dir(free)-K_dir(free, const) * T_dir(const)); 

% Applying Robin BCs for Problem 1.2
K_robin = K;
F_robin = F;
K_robin(1,1) = K_robin(1,1) + h_outside;
F_robin(1) = F_robin(1) + h_outside * T_out;
K_robin(end,end) = K_robin(end,end) + h_inside;
F_robin(end) = F_robin(end) + h_inside * T_in;
T_robin = K_robin \ F_robin;


Length=[0 20 40 60 80 100 120 140 160 180 200 220 222 247 260];


% Plot results

figure;
plot(Length, T_dir, 'b', 'LineWidth', 5); 
hold on;
plot(Length, T_robin, 'r', 'LineWidth', 5);
%plot(Length, T, 'k--', 'LineWidth', 2); % Corrected Analytical Solution
legend('Dirichlet BCs', 'Robin BCs','Location','south');
xlabel('Wall thickness (mm)','FontSize',24); ylabel('Temperature (°C)','FontSize',24);
title('Temperature Distribution Across the Wall','FontSize',24); grid on;
set(gca, 'FontSize', 20)

Q_dir_total = 0; % Initialize total heat loss
Q_robin_total = 0;
for  i = 1:size(Edof,1)
    edof = Edof(i, 2:3); % Get the DOFs for this element
    if i == 1
        k_conduction = k_cond(1);
    elseif i >= 2 && i <= 11
        k_conduction = k_cond(2);
    elseif i == 12
        k_conduction = k_cond(3);
    elseif i == 13
        k_conduction = k_cond(4);
    elseif i == 14
        k_conduction = k_cond(5);
    end
 
    % Compute the heat loss for the element
    Q_dir = -k_conduction * (T_dir(edof(2)) - T_dir(edof(1))) / (Ex(i,2) - Ex(i,1));
    Q_robin = -k_conduction * (T_robin(edof(2)) - T_robin(edof(1))) / (Ex(i,2) - Ex(i,1));
    % Sum up the total heat loss
    q_dir(i)=Q_dir;
    q_robin(i)=Q_robin;
    Q_dir_total = Q_dir_total + Q_dir; 
    Q_robin_total = Q_robin_total + Q_robin;
end
fprintf('Heat loss (Dirichlet): %.2f W/m²\n', Q_dir);
fprintf('Heat loss (Robin): %.2f W/m²\n', Q_robin);


