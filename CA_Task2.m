clear all;
clc;

F=5*10^6; %Applied force
w0=2;
w1=3;
t=2;
h=20;
g=9.81;
rho=2500;
E=37.3*10^9;
L=h;

n_elem =10;         % Number of elements
n_nodes = n_elem + 1;

Edof = zeros(n_elem, 3);

for i = 1:n_elem
    Edof(i, :) = [i, i, i+1];
end

Dof=(1:n_nodes)';

% Generate node coordinates
x_coords = linspace(0, L, n_nodes); 

% Generate Ex matrix: x-coordinates of each element [x_start, x_end]
Ex = [x_coords(1:end-1)', x_coords(2:end)'];


Nodof=size(Dof,1);
K_stiffness=zeros(Nodof);
fb=zeros(Nodof,1);
fl=zeros(Nodof,1);

for i = 1:size(Edof,1)
    ex = Ex(i,:);
    z1=Ex(i,1); 
    z2=Ex(i,2);
    x1=w0+(w1-w0)*z1/h; %Calculation of the element width 1
    x2=w0+(w1-w0)*z2/h; %Calculation of the element width 2
    A1=((x1*t)+(x2*t))/2;
    A_1(i)=A1;%average cross-sectional area of the element 
    Le=Ex(i,2)-Ex(i,1);%Calculation of the element length
    A2=(x1+x2)/2;
    A_2(i)=A2; %element area calculation
    b=A2*t*rho*g; %calculation of the body force
    %b=0;
    % Element stiffness matrix
    Ke = ( A1*E/ Le) * [1 -1; -1 1];
    % Element volume load matrix
    fl_e=b*Le*[0.5;0.5];
    edof = Edof(i, 2:3);
    %Assembly operation
    K_stiffness(edof, edof) = K_stiffness(edof, edof) + Ke;
    fl(edof)=fl(edof)+fl_e;
end 
K=K_stiffness;
free =Dof(1:end-1);%Free DOFs
const = Dof(end);  % Constrained DOFs
fb(1)=fb(1)+F; % Assigning the externl force
F=fb+fl;
U=zeros(Nodof,1);
U(const)=0;  % Apply known displacement
U(free) = K(free, free) \ (F(free)-K(free, const) * U(const)); %Calculation of the unknown displacement

eps= zeros(size(Edof,1),1);
sigma= zeros(size(Edof,1),1);

for  i = 1:size(Edof,1)
    edof = Edof(i, 2:3); % Get the DOFs for this element
    eps(i)=(U(edof(1)) - U(edof(2))) / (Ex(i,2) - Ex(i,1));
    % Compute the stress induced for the element
    sigma(i) = E * eps(i);
end


figure(1)
plot(Dof(:,1),U*10^3,'r-s','LineWidth', 3)
grid on
xlabel('DOFs (each grid is an element)','FontSize',24)
ylabel('Displacement [mm]','FontSize',24)
title('Displacment Distribution Across the Elements','FontSize',24)
set(gca, 'FontSize', 20)


figure(2)
plot(Edof(:,1),sigma*10^-6,'b-o','LineWidth', 3)
grid on
xlabel('Elements','FontSize',24)
ylabel('Stress [MPa]','FontSize',24)
title('Stress Distribution Across the Elements','FontSize',24)
set(gca, 'FontSize', 20)
fprintf('Stress at top (node 1): %.4f MPa\n', sigma(1)*1e-6)
fprintf('Stress at bottom (node %d): %.4f MPa\n', 11, sigma(end)*1e-6)




listelement_list = [1000 800 600 400 200 150 120 100 80 50 40 30 20 15 10];
sigma_max_comprn = [1.25 1.25 1.25 1.25 1.25 1.25 1.25 1.25 1.25 1.248 1.248 1.247 1.246 1.245 1.244];
sigma_x_h = [1.24 1.24 1.24 1.24 1.24 1.24 1.24 1.24 1.24 1.24 1.24 1.240 1.240 1.239 1.238];

sigma_compress_hand = 1.25;
sigma_x_h_hand = 1.24;

figure(3)
plot(listelement_list, sigma_max_comprn, 'r-s', 'LineWidth', 3);
hold on;
plot(listelement_list, sigma_x_h, 'b-o', 'LineWidth', 3);

% Add horizontal lines for hand calculations
yline(sigma_compress_hand, 'r--', 'LineWidth', 2, 'Label', 'Hand calc: max compressive', 'LabelVerticalAlignment', 'bottom');
yline(sigma_x_h_hand, 'b--', 'LineWidth', 2, 'Label', 'Hand calc: stress at x=h', 'LabelVerticalAlignment', 'top');

grid on
xlabel('Number of elements','FontSize',24)
ylabel('Stress [MPa]','FontSize',24)
title('Mesh Convergence Study','FontSize',24)
legend('Maximum compressive stress (FEA)', ...
       'Stress at x=h (FEA)', ...
       'Hand calc: max compressive', ...
       'Hand calc: stress at x=h', ...
       'Location','south');
set(gca, 'FontSize', 20)
