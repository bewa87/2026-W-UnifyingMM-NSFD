%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                %
% NSFD for Irreversible and Reversible Michaelis-Menten Kinetics %
%                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1: Define Problem Parameter and Setup Solution Vectors

% Definition of Problem Parameters

k_01 = 0.15;
k_02 = 0.10;
k_03 = 0.05;
k_04 = 0.5;

% Setup of Time Interval

T = 1000;
h = 1;
t = 0:h:T;

% Setup of Solution Vectors

S = zeros(1,length(t));
E = zeros(1,length(t));
C = zeros(1,length(t));
P = zeros(1,length(t));

% Setup of Initial Conditions

S_null = 0.5;
E_null = 0.2;
C_null = 0.0;
P_null = 0.2;

% Put Initial Conditions into Solution Vectors

S(1) = S_null;
E(1) = E_null;
C(1) = C_null;
P(1) = P_null;

% Step 2: NSFD Solution Loop

for j = 1:1:(length(t)-1)
  F_j    = 1 + h*((k_01 + k_04)*E(j) + k_02 + k_03) + h^2*(k_01*E(j)*(k_03 + k_04*E(j)) + k_02*k_04*E(j));
  G_j    = 1 + h*((k_01 + k_04)*E(j) + k_02 + k_03) + h^2*(k_01*E(j)*(k_02 + k_03 + k_04*E(j)) + k_02*k_04*E(j)) + h^3*(k_01*k_02*k_04*E(j)*E(j));
  S(j+1) = (G_j*S(j))/(F_j*(1 + h*k_01*E(j))) + (G_j*h*k_02*C(j))/(F_j*(1 + h*k_01*E(j))*(1 + h*(k_02 + k_03)))...
           + (G_j*h^3*k_02*k_03*k_04*E(j)*C(j))/(F_j*(1 + h*k_01*E(j))*(1 + h*(k_02 + k_03 + k_04*E(j)) + h^2*(k_02*k_04*E(j)))*(1 + h*(k_02 + k_03)))...
           + (G_j*h^2*k_02*k_04*E(j)*P(j))/(F_j*(1 + h*k_01*E(j))*(1 + h*(k_02 + k_03 + k_04*E(j)) + h^2*(k_02*k_04*E(j))));
  P(j+1) = ((1 + h*(k_02 + k_03))*P(j))/((1 + h*(k_02 + k_03 + k_04*E(j)) + h^2*(k_02*k_04*E(j))))...
           + (h*k_03*C(j))/((1 + h*(k_02 + k_03 + k_04*E(j)) + h^2*(k_02*k_04*E(j))))...
           + (h^2*k_01*k_03*E(j)*S(j+1))/((1 + h*(k_02 + k_03 + k_04*E(j)) + h^2*(k_02*k_04*E(j))));
  C(j+1) = (C(j) + h*E(j)*(k_01*S(j+1) + k_04*P(j+1)))/(1 + h*(k_02 + k_03));
  E(j+1) = E(1) + C(1) - C(j+1);
endfor

% Step 3: Function S*(k_04), E*(k_04), C*(k_04), P*(k_04)

A_01 = E_null + C_null;
A_02 = S_null + C_null + P_null;

S_star_irr = 0;
E_star_irr = E_null + C_null;
C_star_irr = 0;
P_star_irr = S_null + C_null + P_null;

C_star = @(x) (2*A_01*A_02)./(A_01 + A_02 + k_03./x + k_02/k_01 + sqrt((A_01 + A_02 + k_03./x + k_02/k_01).^2 - 4*A_01*A_02));
S_star = @(x) (k_02.*C_star(x))./(k_01*(A_01 - C_star(x)));
E_star = @(x) A_01 - C_star(x);
P_star = @(x) A_02 - S_star(x) - C_star(x);

% Step 4: Plotting of Solutions

figure(1)
plot(t,S,'linewidth',1,'linestyle','-','color','red')
hold on
plot(t,E,'linewidth',1,'linestyle','-','color','blue')
hold on
plot(t,C,'linewidth',1,'linestyle','-','color','magenta')
hold on
plot(t,P,'linewidth',1,'linestyle','-','color','black')
hold on
plot(t,S_star(0.5)*ones(1,length(t)),'linewidth',1,'linestyle','--','color','red')
hold on
plot(t,E_star(0.5)*ones(1,length(t)),'linewidth',1,'linestyle','--','color','blue')
hold on
plot(t,C_star(0.5)*ones(1,length(t)),'linewidth',1,'linestyle','--','color','magenta')
hold on
plot(t,P_star(0.5)*ones(1,length(t)),'linewidth',1,'linestyle','--','color','black')
hold on
title('Concentrations of enzymatic reaction (k_{4} = 0.5)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'S(t)','E(t)','C(t)','P(t)','S^{*}','E^{*}','C^{*}','P^{*}'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 150])
hold off
