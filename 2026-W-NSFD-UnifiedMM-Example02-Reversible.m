%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                %
% NSFD for Irreversible and Reversible Michaelis-Menten Kinetics %
%                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1: Define Problem Parameter and Setup Solution Vectors

% Definition of Problem Parameters

k_01   = 0.15;
k_02   = 0.10;
k_03   = 0.05;
k_04   = 0.01;
k_04_2 = 0.1;
k_04_3 = 1;
k_04_4 = 10;
k_04_5 = 100;

% Setup of Time Interval

T = 1000;
h = 1;
t = 0:h:T;

% Setup of Solution Vectors

S1 = zeros(1,length(t));
E1 = zeros(1,length(t));
C1 = zeros(1,length(t));
P1 = zeros(1,length(t));

S2 = zeros(1,length(t));
E2 = zeros(1,length(t));
C2 = zeros(1,length(t));
P2 = zeros(1,length(t));

S3 = zeros(1,length(t));
E3 = zeros(1,length(t));
C3 = zeros(1,length(t));
P3 = zeros(1,length(t));

S4 = zeros(1,length(t));
E4 = zeros(1,length(t));
C4 = zeros(1,length(t));
P4 = zeros(1,length(t));

S5 = zeros(1,length(t));
E5 = zeros(1,length(t));
C5 = zeros(1,length(t));
P5 = zeros(1,length(t));

% Setup of Initial Conditions

S_null = 0.5;
E_null = 0.2;
C_null = 0.0;
P_null = 0.2;

% Put Initial Conditions into Solution Vectors

S1(1) = S_null;
E1(1) = E_null;
C1(1) = C_null;
P1(1) = P_null;

S2(1) = S_null;
E2(1) = E_null;
C2(1) = C_null;
P2(1) = P_null;

S3(1) = S_null;
E3(1) = E_null;
C3(1) = C_null;
P3(1) = P_null;

S4(1) = S_null;
E4(1) = E_null;
C4(1) = C_null;
P4(1) = P_null;

S5(1) = S_null;
E5(1) = E_null;
C5(1) = C_null;
P5(1) = P_null;

% Step 2: NSFD Solution Loop

% Step 2.1: NSFD Loop for k_04 = 0.01

for j = 1:1:(length(t)-1)
  F_j     = 1 + h*((k_01 + k_04)*E1(j) + k_02 + k_03) + h^2*(k_01*E1(j)*(k_03 + k_04*E1(j)) + k_02*k_04*E1(j));
  G_j     = 1 + h*((k_01 + k_04)*E1(j) + k_02 + k_03) + h^2*(k_01*E1(j)*(k_02 + k_03 + k_04*E1(j)) + k_02*k_04*E1(j)) + h^3*(k_01*k_02*k_04*E1(j)*E1(j));
  S1(j+1) = (G_j*S1(j))/(F_j*(1 + h*k_01*E1(j))) + (G_j*h*k_02*C1(j))/(F_j*(1 + h*k_01*E1(j))*(1 + h*(k_02 + k_03)))...
           + (G_j*h^3*k_02*k_03*k_04*E1(j)*C1(j))/(F_j*(1 + h*k_01*E1(j))*(1 + h*(k_02 + k_03 + k_04*E1(j)) + h^2*(k_02*k_04*E1(j)))*(1 + h*(k_02 + k_03)))...
           + (G_j*h^2*k_02*k_04*E1(j)*P1(j))/(F_j*(1 + h*k_01*E1(j))*(1 + h*(k_02 + k_03 + k_04*E1(j)) + h^2*(k_02*k_04*E1(j))));
  P1(j+1) = ((1 + h*(k_02 + k_03))*P1(j))/((1 + h*(k_02 + k_03 + k_04*E1(j)) + h^2*(k_02*k_04*E1(j))))...
           + (h*k_03*C1(j))/((1 + h*(k_02 + k_03 + k_04*E1(j)) + h^2*(k_02*k_04*E1(j))))...
           + (h^2*k_01*k_03*E1(j)*S1(j+1))/((1 + h*(k_02 + k_03 + k_04*E1(j)) + h^2*(k_02*k_04*E1(j))));
  C1(j+1) = (C1(j) + h*E1(j)*(k_01*S1(j+1) + k_04*P1(j+1)))/(1 + h*(k_02 + k_03));
  E1(j+1) = E1(1) + C1(1) - C1(j+1);
endfor

% Step 2.2: NSFD Loop for k_04_2 = 0.1

for j = 1:1:(length(t)-1)
  F_j     = 1 + h*((k_01 + k_04_2)*E2(j) + k_02 + k_03) + h^2*(k_01*E2(j)*(k_03 + k_04_2*E2(j)) + k_02*k_04_2*E2(j));
  G_j     = 1 + h*((k_01 + k_04_2)*E2(j) + k_02 + k_03) + h^2*(k_01*E2(j)*(k_02 + k_03 + k_04_2*E2(j)) + k_02*k_04_2*E2(j)) + h^3*(k_01*k_02*k_04_2*E2(j)*E2(j));
  S2(j+1) = (G_j*S2(j))/(F_j*(1 + h*k_01*E2(j))) + (G_j*h*k_02*C2(j))/(F_j*(1 + h*k_01*E2(j))*(1 + h*(k_02 + k_03)))...
           + (G_j*h^3*k_02*k_03*k_04_2*E2(j)*C2(j))/(F_j*(1 + h*k_01*E2(j))*(1 + h*(k_02 + k_03 + k_04_2*E2(j)) + h^2*(k_02*k_04_2*E2(j)))*(1 + h*(k_02 + k_03)))...
           + (G_j*h^2*k_02*k_04_2*E2(j)*P2(j))/(F_j*(1 + h*k_01*E2(j))*(1 + h*(k_02 + k_03 + k_04_2*E2(j)) + h^2*(k_02*k_04_2*E2(j))));
  P2(j+1) = ((1 + h*(k_02 + k_03))*P2(j))/((1 + h*(k_02 + k_03 + k_04_2*E2(j)) + h^2*(k_02*k_04_2*E2(j))))...
           + (h*k_03*C2(j))/((1 + h*(k_02 + k_03 + k_04_2*E2(j)) + h^2*(k_02*k_04_2*E2(j))))...
           + (h^2*k_01*k_03*E2(j)*S2(j+1))/((1 + h*(k_02 + k_03 + k_04_2*E2(j)) + h^2*(k_02*k_04_2*E2(j))));
  C2(j+1) = (C2(j) + h*E2(j)*(k_01*S2(j+1) + k_04_2*P2(j+1)))/(1 + h*(k_02 + k_03));
  E2(j+1) = E2(1) + C2(1) - C2(j+1);
endfor

% Step 2.3: NSFD Loop for k_04_3 = 1

for j = 1:1:(length(t)-1)
  F_j     = 1 + h*((k_01 + k_04_3)*E3(j) + k_02 + k_03) + h^2*(k_01*E3(j)*(k_03 + k_04_3*E3(j)) + k_02*k_04_3*E3(j));
  G_j     = 1 + h*((k_01 + k_04_3)*E3(j) + k_02 + k_03) + h^2*(k_01*E3(j)*(k_02 + k_03 + k_04_3*E3(j)) + k_02*k_04_3*E3(j)) + h^3*(k_01*k_02*k_04_3*E3(j)*E3(j));
  S3(j+1) = (G_j*S3(j))/(F_j*(1 + h*k_01*E3(j))) + (G_j*h*k_02*C3(j))/(F_j*(1 + h*k_01*E3(j))*(1 + h*(k_02 + k_03)))...
           + (G_j*h^3*k_02*k_03*k_04_3*E3(j)*C3(j))/(F_j*(1 + h*k_01*E3(j))*(1 + h*(k_02 + k_03 + k_04_3*E3(j)) + h^2*(k_02*k_04_3*E3(j)))*(1 + h*(k_02 + k_03)))...
           + (G_j*h^2*k_02*k_04_3*E3(j)*P3(j))/(F_j*(1 + h*k_01*E3(j))*(1 + h*(k_02 + k_03 + k_04_3*E3(j)) + h^2*(k_02*k_04_3*E3(j))));
  P3(j+1) = ((1 + h*(k_02 + k_03))*P3(j))/((1 + h*(k_02 + k_03 + k_04_3*E3(j)) + h^2*(k_02*k_04_3*E3(j))))...
           + (h*k_03*C3(j))/((1 + h*(k_02 + k_03 + k_04_3*E3(j)) + h^2*(k_02*k_04_3*E3(j))))...
           + (h^2*k_01*k_03*E3(j)*S3(j+1))/((1 + h*(k_02 + k_03 + k_04_3*E3(j)) + h^2*(k_02*k_04_3*E3(j))));
  C3(j+1) = (C3(j) + h*E3(j)*(k_01*S3(j+1) + k_04_3*P3(j+1)))/(1 + h*(k_02 + k_03));
  E3(j+1) = E3(1) + C3(1) - C3(j+1);
endfor

% Step 2.4: NSFD Loop for k_04_4 = 10

for j = 1:1:(length(t)-1)
  F_j     = 1 + h*((k_01 + k_04_4)*E4(j) + k_02 + k_03) + h^2*(k_01*E4(j)*(k_03 + k_04_4*E4(j)) + k_02*k_04_4*E4(j));
  G_j     = 1 + h*((k_01 + k_04_4)*E4(j) + k_02 + k_03) + h^2*(k_01*E4(j)*(k_02 + k_03 + k_04_4*E4(j)) + k_02*k_04_4*E4(j)) + h^3*(k_01*k_02*k_04_4*E4(j)*E4(j));
  S4(j+1) = (G_j*S4(j))/(F_j*(1 + h*k_01*E4(j))) + (G_j*h*k_02*C4(j))/(F_j*(1 + h*k_01*E4(j))*(1 + h*(k_02 + k_03)))...
           + (G_j*h^3*k_02*k_03*k_04_4*E4(j)*C4(j))/(F_j*(1 + h*k_01*E4(j))*(1 + h*(k_02 + k_03 + k_04_4*E4(j)) + h^2*(k_02*k_04_4*E4(j)))*(1 + h*(k_02 + k_03)))...
           + (G_j*h^2*k_02*k_04_4*E4(j)*P4(j))/(F_j*(1 + h*k_01*E4(j))*(1 + h*(k_02 + k_03 + k_04_4*E4(j)) + h^2*(k_02*k_04_4*E4(j))));
  P4(j+1) = ((1 + h*(k_02 + k_03))*P4(j))/((1 + h*(k_02 + k_03 + k_04_4*E4(j)) + h^2*(k_02*k_04_4*E4(j))))...
           + (h*k_03*C4(j))/((1 + h*(k_02 + k_03 + k_04_4*E4(j)) + h^2*(k_02*k_04_4*E4(j))))...
           + (h^2*k_01*k_03*E4(j)*S4(j+1))/((1 + h*(k_02 + k_03 + k_04_4*E4(j)) + h^2*(k_02*k_04_4*E4(j))));
  C4(j+1) = (C4(j) + h*E4(j)*(k_01*S4(j+1) + k_04_4*P4(j+1)))/(1 + h*(k_02 + k_03));
  E4(j+1) = E4(1) + C4(1) - C4(j+1);
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

% Step 4.1: Numerical Solution for k_04 = 0.01

figure(1)
plot(t,S1,'linewidth',1,'linestyle','-','color','red')
hold on
plot(t,E1,'linewidth',1,'linestyle','-','color','blue')
hold on
plot(t,C1,'linewidth',1,'linestyle','-','color','magenta')
hold on
plot(t,P1,'linewidth',1,'linestyle','-','color','black')
hold on
plot(t,S_star(0.01)*ones(1,length(t)),'linewidth',1,'linestyle','--','color','red')
hold on
plot(t,E_star(0.01)*ones(1,length(t)),'linewidth',1,'linestyle','--','color','blue')
hold on
plot(t,C_star(0.01)*ones(1,length(t)),'linewidth',1,'linestyle','--','color','magenta')
hold on
plot(t,P_star(0.01)*ones(1,length(t)),'linewidth',1,'linestyle','--','color','black')
hold on
title('Concentrations of enzymatic reaction (k_{4} = 0.01)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'S(t)','E(t)','C(t)','P(t)','S^{*}','E^{*}','C^{*}','P^{*}'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 T])
hold off

% Step 4.2: Numerical Solution for k_04 = 0.1

figure(2)
plot(t,S2,'linewidth',1,'linestyle','-','color','red')
hold on
plot(t,E2,'linewidth',1,'linestyle','-','color','blue')
hold on
plot(t,C2,'linewidth',1,'linestyle','-','color','magenta')
hold on
plot(t,P2,'linewidth',1,'linestyle','-','color','black')
hold on
plot(t,S_star(0.1)*ones(1,length(t)),'linewidth',1,'linestyle','--','color','red')
hold on
plot(t,E_star(0.1)*ones(1,length(t)),'linewidth',1,'linestyle','--','color','blue')
hold on
plot(t,C_star(0.1)*ones(1,length(t)),'linewidth',1,'linestyle','--','color','magenta')
hold on
plot(t,P_star(0.1)*ones(1,length(t)),'linewidth',1,'linestyle','--','color','black')
hold on
title('Concentrations of enzymatic reaction (k_{4} = 0.1)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'S(t)','E(t)','C(t)','P(t)','S^{*}','E^{*}','C^{*}','P^{*}'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 500])
hold off

% Step 4.3: Numerical Solution for k_04_3 = 1

figure(3)
plot(t,S3,'linewidth',1,'linestyle','-','color','red')
hold on
plot(t,E3,'linewidth',1,'linestyle','-','color','blue')
hold on
plot(t,C3,'linewidth',1,'linestyle','-','color','magenta')
hold on
plot(t,P3,'linewidth',1,'linestyle','-','color','black')
hold on
plot(t,S_star(1)*ones(1,length(t)),'linewidth',1,'linestyle','--','color','red')
hold on
plot(t,E_star(1)*ones(1,length(t)),'linewidth',1,'linestyle','--','color','blue')
hold on
plot(t,C_star(1)*ones(1,length(t)),'linewidth',1,'linestyle','--','color','magenta')
hold on
plot(t,P_star(1)*ones(1,length(t)),'linewidth',1,'linestyle','--','color','black')
hold on
title('Concentrations of enzymatic reaction (k_{4} = 1)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'S(t)','E(t)','C(t)','P(t)','S^{*}','E^{*}','C^{*}','P^{*}'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 300])
hold off

% Step 4.4: Numerical Solution for k_04_4 = 10

figure(4)
plot(t,S4,'linewidth',1,'linestyle','-','color','red')
hold on
plot(t,E4,'linewidth',1,'linestyle','-','color','blue')
hold on
plot(t,C4,'linewidth',1,'linestyle','-','color','magenta')
hold on
plot(t,P4,'linewidth',1,'linestyle','-','color','black')
hold on
plot(t,S_star(10)*ones(1,length(t)),'linewidth',1,'linestyle','--','color','red')
hold on
plot(t,E_star(10)*ones(1,length(t)),'linewidth',1,'linestyle','--','color','blue')
hold on
plot(t,C_star(10)*ones(1,length(t)),'linewidth',1,'linestyle','--','color','magenta')
hold on
plot(t,P_star(10)*ones(1,length(t)),'linewidth',1,'linestyle','--','color','black')
hold on
title('Concentrations of enzymatic reaction (k_{4} = 10)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'S(t)','E(t)','C(t)','P(t)','S^{*}','E^{*}','C^{*}','P^{*}'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 200])
hold off

% Step 4.5: Numerical Equilibrium States vs. Theoretical Ones for S*(k_{4})

figure(5)
plot(0:0.1:10,S_star(0:0.1:10),'linewidth',1,'linestyle','-','color','red')
hold on
plot(0.01,S1(end),'marker','*','red')
hold on
plot(0.1,S2(end),'marker','*','red')
hold on
plot(1,S3(end),'marker','*','red')
hold on
plot(10,S4(end),'marker','*','red')
hold on
title('Theoretical S*(k_{4})-Values and Numerical Equilibrium Values','interpreter','tex','fontsize',14,'fontweight','normal')
xlabel('k_{4}','interpreter','tex','fontsize',12)
ylabel('S*(k_{4})','interpreter','tex','fontsize',12)
xlim([0 10])
hold off

% Step 4.6: Numerical Equilibrium States vs. Theoretical Ones for E*(k_{4})

figure(6)
plot(0:0.1:10,E_star(0:0.1:10),'linewidth',1,'linestyle','-','color','blue')
hold on
plot(0.01,E1(end),'marker','*','blue')
hold on
plot(0.1,E2(end),'marker','*','blue')
hold on
plot(1,E3(end),'marker','*','blue')
hold on
plot(10,E4(end),'marker','*','blue')
hold on
title('Theoretical E*(k_{4})-Values and Numerical Equilibrium Values','interpreter','tex','fontsize',14,'fontweight','normal')
xlabel('k_{4}','interpreter','tex','fontsize',12)
ylabel('E*(k_{4})','interpreter','tex','fontsize',12)
xlim([0 10])
hold off

% Step 4.7: Numerical Equilibrium States vs. Theoretical Ones for C*(k_{4})

figure(7)
plot(0:0.1:10,C_star(0:0.1:10),'linewidth',1,'linestyle','-','color','magenta')
hold on
plot(0.01,C1(end),'marker','*','magenta')
hold on
plot(0.1,C2(end),'marker','*','magenta')
hold on
plot(1,C3(end),'marker','*','magenta')
hold on
plot(10,C4(end),'marker','*','magenta')
hold on
title('Theoretical C*(k_{4})-Values and Numerical Equilibrium Values','interpreter','tex','fontsize',14,'fontweight','normal')
xlabel('k_{4}','interpreter','tex','fontsize',12)
ylabel('C*(k_{4})','interpreter','tex','fontsize',12)
xlim([0 10])
hold off

% Step 4.8: Numerical Equilibrium States vs. Theoretical Ones for P*(k_{4})

figure(8)
plot(0:0.1:10,P_star(0:0.1:10),'linewidth',1,'linestyle','-','color','black')
hold on
plot(0.01,P1(end),'marker','*','black')
hold on
plot(0.1,P2(end),'marker','*','black')
hold on
plot(1,P3(end),'marker','*','black')
hold on
plot(10,P4(end),'marker','*','black')
hold on
title('Theoretical P*(k_{4})-Values and Numerical Equilibrium Values','interpreter','tex','fontsize',14,'fontweight','normal')
xlabel('k_{4}','interpreter','tex','fontsize',12)
ylabel('P*(k_{4})','interpreter','tex','fontsize',12)
xlim([0 10])
hold off
