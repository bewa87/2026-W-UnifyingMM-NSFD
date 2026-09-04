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
k_04 = 0.0;

% Setup of Time Intervals

T1_exp = 200;
T1_imp = 200000;
T2     = 1000;

% Time Interval 1: Large Time Step

h1_exp = 20;
h1_imp = 1000;
t1_exp = 0:h1_exp:T1_exp;
t1_imp = 0:h1_imp:T1_imp;

% Time Interval 2: Small Time Step

h2 = 0.1;
t2 = 0:h2:T2;

% Setup of Different Solution Vectors For Time-Stepping Methods

% Solution Vectors 1: Large Time Step

S1      = zeros(1,length(t1_exp));
E1      = zeros(1,length(t1_exp));
C1      = zeros(1,length(t1_exp));
P1      = zeros(1,length(t1_exp));

S1_imp  = zeros(1,length(t1_imp));
E1_imp  = zeros(1,length(t1_imp));
C1_imp  = zeros(1,length(t1_imp));
P1_imp  = zeros(1,length(t1_imp));

S1_RK2  = zeros(1,length(t1_exp));
E1_RK2  = zeros(1,length(t1_exp));
C1_RK2  = zeros(1,length(t1_exp));
P1_RK2  = zeros(1,length(t1_exp));

S1_RK4  = zeros(1,length(t1_exp));
E1_RK4  = zeros(1,length(t1_exp));
C1_RK4  = zeros(1,length(t1_exp));
P1_RK4  = zeros(1,length(t1_exp));

S1_IE   = zeros(1,length(t1_imp));
E1_IE   = zeros(1,length(t1_imp));
C1_IE   = zeros(1,length(t1_imp));
P1_IE   = zeros(1,length(t1_imp));

S1_BDF2 = zeros(1,length(t1_imp));
E1_BDF2 = zeros(1,length(t1_imp));
C1_BDF2 = zeros(1,length(t1_imp));
P1_BDF2 = zeros(1,length(t1_imp));

% Solution Vectors 2: Small Time Step

S2      = zeros(1,length(t2));
E2      = zeros(1,length(t2));
C2      = zeros(1,length(t2));
P2      = zeros(1,length(t2));

S2_RK2  = zeros(1,length(t2));
E2_RK2  = zeros(1,length(t2));
C2_RK2  = zeros(1,length(t2));
P2_RK2  = zeros(1,length(t2));

S2_RK4  = zeros(1,length(t2));
E2_RK4  = zeros(1,length(t2));
C2_RK4  = zeros(1,length(t2));
P2_RK4  = zeros(1,length(t2));

S2_IE   = zeros(1,length(t2));
E2_IE   = zeros(1,length(t2));
C2_IE   = zeros(1,length(t2));
P2_IE   = zeros(1,length(t2));

S2_BDF2 = zeros(1,length(t2));
E2_BDF2 = zeros(1,length(t2));
C2_BDF2 = zeros(1,length(t2));
P2_BDF2 = zeros(1,length(t2));

% Setup of Initial Conditions

S_null = 0.5;
E_null = 0.2;
C_null = 0.0;
P_null = 0.2;

% Put Initial Conditions into Solution Vectors 1: Large Time Step

S1(1)      = S_null;
E1(1)      = E_null;
C1(1)      = C_null;
P1(1)      = P_null;

S1_imp(1)  = S_null;
E1_imp(1)  = E_null;
C1_imp(1)  = C_null;
P1_imp(1)  = P_null;

S1_RK2(1)  = S_null;
E1_RK2(1)  = E_null;
C1_RK2(1)  = C_null;
P1_RK2(1)  = P_null;

S1_RK4(1)  = S_null;
E1_RK4(1)  = E_null;
C1_RK4(1)  = C_null;
P1_RK4(1)  = P_null;

S1_IE(1)   = S_null;
E1_IE(1)   = E_null;
C1_IE(1)   = C_null;
P1_IE(1)   = P_null;

S1_BDF2(1) = S_null;
E1_BDF2(1) = E_null;
C1_BDF2(1) = C_null;
P1_BDF2(1) = P_null;

% Put Initial Conditions into Solution Vectors 2: Small Time Step

S2(1)      = S_null;
E2(1)      = E_null;
C2(1)      = C_null;
P2(1)      = P_null;

S2_RK2(1)  = S_null;
E2_RK2(1)  = E_null;
C2_RK2(1)  = C_null;
P2_RK2(1)  = P_null;

S2_RK4(1)  = S_null;
E2_RK4(1)  = E_null;
C2_RK4(1)  = C_null;
P2_RK4(1)  = P_null;

S2_IE(1)   = S_null;
E2_IE(1)   = E_null;
C2_IE(1)   = C_null;
P2_IE(1)   = P_null;

S2_BDF2(1) = S_null;
E2_BDF2(1) = E_null;
C2_BDF2(1) = C_null;
P2_BDF2(1) = P_null;

% Step 2: Solution Loops

% Step 2.1: NSFD Solution Loop 1: Large Time Step

for j = 1:1:(length(t1_exp)-1)
  F_j    = 1 + h1_exp*((k_01 + k_04)*E1(j) + k_02 + k_03) + (h1_exp)^2*(k_01*E1(j)*(k_03 + k_04*E1(j)) + k_02*k_04*E1(j));
  G_j    = 1 + h1_exp*((k_01 + k_04)*E1(j) + k_02 + k_03) + (h1_exp)^2*(k_01*E1(j)*(k_02 + k_03 + k_04*E1(j)) + k_02*k_04*E1(j)) + (h1_exp)^3*(k_01*k_02*k_04*E1(j)*E1(j));
  S1(j+1) = (G_j*S1(j))/(F_j*(1 + h1_exp*k_01*E1(j))) + (G_j*h1_exp*k_02*C1(j))/(F_j*(1 + h1_exp*k_01*E1(j))*(1 + h1_exp*(k_02 + k_03)))...
           + (G_j*(h1_exp)^3*k_02*k_03*k_04*E1(j)*C1(j))/(F_j*(1 + h1_exp*k_01*E1(j))*(1 + h1_exp*(k_02 + k_03 + k_04*E1(j)) + (h1_exp)^2*(k_02*k_04*E1(j)))*(1 + h1_exp*(k_02 + k_03)))...
           + (G_j*(h1_exp)^2*k_02*k_04*E1(j)*P1(j))/(F_j*(1 + h1_exp*k_01*E1(j))*(1 + h1_exp*(k_02 + k_03 + k_04*E1(j)) + (h1_exp)^2*(k_02*k_04*E1(j))));
  P1(j+1) = ((1 + h1_exp*(k_02 + k_03))*P1(j))/((1 + h1_exp*(k_02 + k_03 + k_04*E1(j)) + (h1_exp)^2*(k_02*k_04*E1(j))))...
           + (h1_exp*k_03*C1(j))/((1 + h1_exp*(k_02 + k_03 + k_04*E1(j)) + (h1_exp)^2*(k_02*k_04*E1(j))))...
           + ((h1_exp)^2*k_01*k_03*E1(j)*S1(j+1))/((1 + h1_exp*(k_02 + k_03 + k_04*E1(j)) + (h1_exp)^2*(k_02*k_04*E1(j))));
  C1(j+1) = (C1(j) + h1_exp*E1(j)*(k_01*S1(j+1) + k_04*P1(j+1)))/(1 + h1_exp*(k_02 + k_03));
  E1(j+1) = E1(1) + C1(1) - C1(j+1);
endfor

% Step 2.2: RK2 Solution Loop 1: Large Time Step

f_S = @(S,C,P) -k_01*S*(E_null + C_null - C) + k_02*C;
f_C = @(S,C,P) k_01*S*(E_null + C_null - C) - (k_02 + k_03)*C + k_04*P*(E_null + C_null - C);
f_P = @(S,C,P) k_03*C - k_04*P*(E_null + C_null - C);

for j = 1:1:(length(t1_exp)-1)
  k1_S   = f_S(S1_RK2(j),C1_RK2(j),P1_RK2(j));
  k1_C   = f_C(S1_RK2(j),C1_RK2(j),P1_RK2(j));
  k1_P   = f_P(S1_RK2(j),C1_RK2(j),P1_RK2(j));

  k2_S   = f_S(S1_RK2(j)+h1_exp*k1_S,C1_RK2(j)+h1_exp*k1_C,P1_RK2(j)+h1_exp*k1_P);
  k2_C   = f_C(S1_RK2(j)+h1_exp*k1_S,C1_RK2(j)+h1_exp*k1_C,P1_RK2(j)+h1_exp*k1_P);
  k2_P   = f_P(S1_RK2(j)+h1_exp*k1_S,C1_RK2(j)+h1_exp*k1_C,P1_RK2(j)+h1_exp*k1_P);

  S1_RK2(j+1) = S1_RK2(j) + (h1_exp/2)*(k1_S + k2_S);
  C1_RK2(j+1) = C1_RK2(j) + (h1_exp/2)*(k1_C + k2_C);
  P1_RK2(j+1) = P1_RK2(j) + (h1_exp/2)*(k1_P + k2_P);
  E1_RK2(j+1) = E_null + C_null - C1_RK2(j+1);
endfor

% Step 2.3: RK4 Solution Loop 1: Large Time Step

for j = 1:1:(length(t1_exp)-1)
  k1_S   = f_S(S1_RK4(j),C1_RK4(j),P1_RK4(j));
  k1_C   = f_C(S1_RK4(j),C1_RK4(j),P1_RK4(j));
  k1_P   = f_P(S1_RK4(j),C1_RK4(j),P1_RK4(j));

  k2_S   = f_S(S1_RK4(j)+(h1_exp/2)*k1_S,C1_RK4(j)+(h1_exp/2)*k1_C,P1_RK4(j)+(h1_exp/2)*k1_P);
  k2_C   = f_C(S1_RK4(j)+(h1_exp/2)*k1_S,C1_RK4(j)+(h1_exp/2)*k1_C,P1_RK4(j)+(h1_exp/2)*k1_P);
  k2_P   = f_P(S1_RK4(j)+(h1_exp/2)*k1_S,C1_RK4(j)+(h1_exp/2)*k1_C,P1_RK4(j)+(h1_exp/2)*k1_P);

  k3_S   = f_S(S1_RK4(j)+(h1_exp/2)*k2_S,C1_RK4(j)+(h1_exp/2)*k2_C,P1_RK4(j)+(h1_exp/2)*k2_P);
  k3_C   = f_C(S1_RK4(j)+(h1_exp/2)*k2_S,C1_RK4(j)+(h1_exp/2)*k2_C,P1_RK4(j)+(h1_exp/2)*k2_P);
  k3_P   = f_P(S1_RK4(j)+(h1_exp/2)*k2_S,C1_RK4(j)+(h1_exp/2)*k2_C,P1_RK4(j)+(h1_exp/2)*k2_P);

  k4_S   = f_S(S1_RK4(j)+h1_exp*k3_S,C1_RK4(j)+h1_exp*k3_C,P1_RK4(j)+h1_exp*k3_P);
  k4_C   = f_C(S1_RK4(j)+h1_exp*k3_S,C1_RK4(j)+h1_exp*k3_C,P1_RK4(j)+h1_exp*k3_P);
  k4_P   = f_P(S1_RK4(j)+h1_exp*k3_S,C1_RK4(j)+h1_exp*k3_C,P1_RK4(j)+h1_exp*k3_P);

  S1_RK4(j+1) = S1_RK4(j) + (h1_exp/6)*(k1_S + 2*k2_S + 2*k3_S + k4_S);
  C1_RK4(j+1) = C1_RK4(j) + (h1_exp/6)*(k1_C + 2*k2_C + 2*k3_C + k4_C);
  P1_RK4(j+1) = P1_RK4(j) + (h1_exp/6)*(k1_P + 2*k2_P + 2*k3_P + k4_P);
  E1_RK4(j+1) = E_null + C_null - C1_RK4(j+1);
endfor

% Step 2.4A: NSFD Solution Loop 1: Large Time Step

for j = 1:1:(length(t1_imp)-1)
  F_j         = 1 + h1_imp*((k_01 + k_04)*E1_imp(j) + k_02 + k_03) + (h1_imp)^2*(k_01*E1_imp(j)*(k_03 + k_04*E1_imp(j)) + k_02*k_04*E1_imp(j));
  G_j         = 1 + h1_imp*((k_01 + k_04)*E1_imp(j) + k_02 + k_03) + (h1_imp)^2*(k_01*E1_imp(j)*(k_02 + k_03 + k_04*E1_imp(j)) + k_02*k_04*E1_imp(j)) + (h1_imp)^3*(k_01*k_02*k_04*E1_imp(j)*E1_imp(j));
  S1_imp(j+1) = (G_j*S1_imp(j))/(F_j*(1 + h1_imp*k_01*E1_imp(j))) + (G_j*h1_imp*k_02*C1_imp(j))/(F_j*(1 + h1_imp*k_01*E1_imp(j))*(1 + h1_imp*(k_02 + k_03)))...
              + (G_j*(h1_imp)^3*k_02*k_03*k_04*E1_imp(j)*C1_imp(j))/(F_j*(1 + h1_imp*k_01*E1_imp(j))*(1 + h1_imp*(k_02 + k_03 + k_04*E1_imp(j)) + (h1_imp)^2*(k_02*k_04*E1_imp(j)))*(1 + h1_imp*(k_02 + k_03)))...
              + (G_j*(h1_imp)^2*k_02*k_04*E1_imp(j)*P1_imp(j))/(F_j*(1 + h1_imp*k_01*E1_imp(j))*(1 + h1_imp*(k_02 + k_03 + k_04*E1_imp(j)) + (h1_imp)^2*(k_02*k_04*E1_imp(j))));
  P1_imp(j+1) = ((1 + h1_imp*(k_02 + k_03))*P1_imp(j))/((1 + h1_imp*(k_02 + k_03 + k_04*E1_imp(j)) + (h1_imp)^2*(k_02*k_04*E1_imp(j))))...
              + (h1_imp*k_03*C1_imp(j))/((1 + h1_imp*(k_02 + k_03 + k_04*E1_imp(j)) + (h1_imp)^2*(k_02*k_04*E1_imp(j))))...
              + ((h1_imp)^2*k_01*k_03*E1_imp(j)*S1_imp(j+1))/((1 + h1_imp*(k_02 + k_03 + k_04*E1_imp(j)) + (h1_imp)^2*(k_02*k_04*E1_imp(j))));
  C1_imp(j+1) = (C1_imp(j) + h1_imp*E1_imp(j)*(k_01*S1_imp(j+1) + k_04*P1_imp(j+1)))/(1 + h1_imp*(k_02 + k_03));
  E1_imp(j+1) = E1_imp(1) + C1_imp(1) - C1_imp(j+1);
endfor

% Step 2.4B: IE Solution Loop 1: Large Time Step

F1    = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) S_new + h1_imp*k_01*S_new*(E_null+C_null-C_new) - h1_imp*k_02*C_new - S_prev;
F2    = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) -h1_imp*k_01*S_new*(E_null+C_null-C_new) + (1+h1_imp*k_02+h1_imp*k_03)*C_new - h1_imp*k_04*P_new*(E_null+C_null-C_new) - C_prev;
F3    = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) P_new -h1_imp*k_03*C_new + h1_imp*k_04*P_new*(E_null+C_null-C_new) - P_prev;
F_vec = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) [F1(S_new,C_new,P_new,S_prev,C_prev,P_prev);F2(S_new,C_new,P_new,S_prev,C_prev,P_prev);F3(S_new,C_new,P_new,S_prev,C_prev,P_prev)];

F1_dS = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) 1 + h1_imp*k_01*(E_null+C_null-C_new);
F1_dC = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) -h1_imp*k_01*S_new - h1_imp*k_02;
F1_dP = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) 0;
F2_dS = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) -h1_imp*k_01*(E_null+C_null-C_new);
F2_dC = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) 1 + h1_imp*k_01*S_new + h1_imp*k_02 + h1_imp*k_03 + h1_imp*k_04*P_new;
F2_dP = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) -h1_imp*k_04*(E_null+C_null-C_new);
F3_dS = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) 0;
F3_dC = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) -h1_imp*k_03 - h1_imp*k_04*P_new;
F3_dP = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) 1 + h1_imp*k_04*(E_null+C_null-C_new);

JF_mat  = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) [F1_dS(S_new,C_new,P_new,S_prev,C_prev,P_prev) F1_dC(S_new,C_new,P_new,S_prev,C_prev,P_prev) F1_dP(S_new,C_new,P_new,S_prev,C_prev,P_prev); F2_dS(S_new,C_new,P_new,S_prev,C_prev,P_prev) F2_dC(S_new,C_new,P_new,S_prev,C_prev,P_prev) F2_dP(S_new,C_new,P_new,S_prev,C_prev,P_prev); F3_dS(S_new,C_new,P_new,S_prev,C_prev,P_prev) F3_dC(S_new,C_new,P_new,S_prev,C_prev,P_prev) F3_dP(S_new,C_new,P_new,S_prev,C_prev,P_prev)];
X_curr  = [S1_IE(1);C1_IE(1);P1_IE(1)];

for l = 1:1:(length(t1_imp)-1)
  for j = 1:1:50
    JF_step = JF_mat(X_curr(1),X_curr(2),X_curr(3),S1_IE(l),C1_IE(l),P1_IE(l));
    H       = -inv(JF_step)*F_vec(X_curr(1),X_curr(2),X_curr(3),S1_IE(l),C1_IE(l),P1_IE(l));
    X_curr  = X_curr + H;
  endfor
  S1_IE(l+1) = X_curr(1);
  C1_IE(l+1) = X_curr(2);
  P1_IE(l+1) = X_curr(3);
  E1_IE(l+1) = E_null+C_null-C1_IE(l+1);
endfor

% Step 2.4C: BDF2 Solution Loop 1: Large Time Step

JF_mat  = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) [F1_dS(S_new,C_new,P_new,S_prev,C_prev,P_prev) F1_dC(S_new,C_new,P_new,S_prev,C_prev,P_prev) F1_dP(S_new,C_new,P_new,S_prev,C_prev,P_prev); F2_dS(S_new,C_new,P_new,S_prev,C_prev,P_prev) F2_dC(S_new,C_new,P_new,S_prev,C_prev,P_prev) F2_dP(S_new,C_new,P_new,S_prev,C_prev,P_prev); F3_dS(S_new,C_new,P_new,S_prev,C_prev,P_prev) F3_dC(S_new,C_new,P_new,S_prev,C_prev,P_prev) F3_dP(S_new,C_new,P_new,S_prev,C_prev,P_prev)];
X_curr  = [S1_BDF2(1);C1_BDF2(1);P1_BDF2(1)];

for j = 1:1:100
  JF_step = JF_mat(X_curr(1),X_curr(2),X_curr(3),S1_BDF2(1),C1_BDF2(1),P1_BDF2(1));
  H       = -inv(JF_step)*F_vec(X_curr(1),X_curr(2),X_curr(3),S1_BDF2(1),C1_BDF2(1),P1_BDF2(1));
  X_curr  = X_curr + H;
endfor

S1_BDF2(2) = X_curr(1);
C1_BDF2(2) = X_curr(2);
P1_BDF2(2) = X_curr(3);
E1_BDF2(2) = E_null + C_null - C1_BDF2(2);

F1_BDF2    = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) -S_new + (4/3)*S_curr - (1/3)*S_old + (2*h1_imp/3)*(-k_01*S_new*(E_null+C_null-C_new) + k_02*C_new);
F2_BDF2    = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) -C_new + (4/3)*C_curr - (1/3)*C_old + (2*h1_imp/3)*(k_01*S_new*(E_null+C_null-C_new) - (k_02+k_03)*C_new + k_04*P_new*(E_null+C_null-C_new));
F3_BDF2    = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) -P_new + (4/3)*P_curr - (1/3)*P_old + (2*h1_imp/3)*(k_03*C_new - k_04*P_new*(E_null+C_null-C_new));
F_vec_BDF2 = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) [F1_BDF2(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old);F2_BDF2(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old);F3_BDF2(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old)];

F1_BDF2_dS = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) -1 - (2/3)*h1_imp*k_01*(E_null+C_null-C_new);
F1_BDF2_dC = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) (2/3)*h1_imp*k_01*S_new + (2/3)*h1_imp*k_02;
F1_BDF2_dP = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) 0;
F2_BDF2_dS = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) (2/3)*h1_imp*k_01*(E_null+C_null-C_new);
F2_BDF2_dC = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) -1 - (2/3)*h1_imp*k_01*S_new - (2/3)*h1_imp*(k_02+k_03) - (2/3)*h1_imp*k_04*P_new;
F2_BDF2_dP = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) (2/3)*h1_imp*k_04*(E_null+C_null-C_new);
F3_BDF2_dS = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) 0;
F3_BDF2_dC = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) (2/3)*h1_imp*k_03 + (2/3)*h1_imp*k_04*P_new;
F3_BDF2_dP = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) -1 - (2/3)*h1_imp*k_04*(E_null+C_null-C_new);

JF_mat_BDF2 = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) [F1_BDF2_dS(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) F1_BDF2_dC(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) F1_BDF2_dP(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old);F2_BDF2_dS(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) F2_BDF2_dC(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) F2_BDF2_dP(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old);F3_BDF2_dS(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) F3_BDF2_dC(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) F3_BDF2_dP(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old)];

for l = 1:1:(length(t1_imp)-2)
  X_curr  = [S1_BDF2(l+1);C1_BDF2(l+1);P1_BDF2(l+1)];
  for j = 1:1:100
    JF_step = JF_mat_BDF2(X_curr(1),X_curr(2),X_curr(3),S1_BDF2(l+1),C1_BDF2(l+1),P1_BDF2(l+1),S1_BDF2(l),C1_BDF2(l),P1_BDF2(l));
    H       = -inv(JF_step)*F_vec_BDF2(X_curr(1),X_curr(2),X_curr(3),S1_BDF2(l+1),C1_BDF2(l+1),P1_BDF2(l+1),S1_BDF2(l),C1_BDF2(l),P1_BDF2(l));
    X_curr  = X_curr + H;
  endfor
  S1_BDF2(l+2) = X_curr(1);
  C1_BDF2(l+2) = X_curr(2);
  P1_BDF2(l+2) = X_curr(3);
  E1_BDF2(l+2) = E_null+C_null-C1_BDF2(l+2);
endfor

% Step 2.5: NSFD Solution Loop 1: Small Time Step

for j = 1:1:(length(t2)-1)
  F_j    = 1 + h2*((k_01 + k_04)*E2(j) + k_02 + k_03) + h2^2*(k_01*E2(j)*(k_03 + k_04*E2(j)) + k_02*k_04*E2(j));
  G_j    = 1 + h2*((k_01 + k_04)*E2(j) + k_02 + k_03) + h2^2*(k_01*E2(j)*(k_02 + k_03 + k_04*E2(j)) + k_02*k_04*E2(j)) + h2^3*(k_01*k_02*k_04*E2(j)*E2(j));
  S2(j+1) = (G_j*S2(j))/(F_j*(1 + h2*k_01*E2(j))) + (G_j*h2*k_02*C2(j))/(F_j*(1 + h2*k_01*E2(j))*(1 + h2*(k_02 + k_03)))...
           + (G_j*h2^3*k_02*k_03*k_04*E2(j)*C2(j))/(F_j*(1 + h2*k_01*E2(j))*(1 + h2*(k_02 + k_03 + k_04*E2(j)) + h2^2*(k_02*k_04*E2(j)))*(1 + h2*(k_02 + k_03)))...
           + (G_j*h2^2*k_02*k_04*E2(j)*P2(j))/(F_j*(1 + h2*k_01*E2(j))*(1 + h2*(k_02 + k_03 + k_04*E2(j)) + h2^2*(k_02*k_04*E2(j))));
  P2(j+1) = ((1 + h2*(k_02 + k_03))*P2(j))/((1 + h2*(k_02 + k_03 + k_04*E2(j)) + h2^2*(k_02*k_04*E2(j))))...
           + (h2*k_03*C2(j))/((1 + h2*(k_02 + k_03 + k_04*E2(j)) + h2^2*(k_02*k_04*E2(j))))...
           + (h2^2*k_01*k_03*E2(j)*S2(j+1))/((1 + h2*(k_02 + k_03 + k_04*E2(j)) + h2^2*(k_02*k_04*E2(j))));
  C2(j+1) = (C2(j) + h2*E2(j)*(k_01*S2(j+1) + k_04*P2(j+1)))/(1 + h2*(k_02 + k_03));
  E2(j+1) = E2(1) + C2(1) - C2(j+1);
endfor

% Step 2.6: RK2 Solution Loop 2: Small Time Step

for j = 1:1:(length(t2)-1)
  k1_S   = f_S(S2_RK2(j),C2_RK2(j),P2_RK2(j));
  k1_C   = f_C(S2_RK2(j),C2_RK2(j),P2_RK2(j));
  k1_P   = f_P(S2_RK2(j),C2_RK2(j),P2_RK2(j));

  k2_S   = f_S(S2_RK2(j)+h2*k1_S,C2_RK2(j)+h2*k1_C,P2_RK2(j)+h2*k1_P);
  k2_C   = f_C(S2_RK2(j)+h2*k1_S,C2_RK2(j)+h2*k1_C,P2_RK2(j)+h2*k1_P);
  k2_P   = f_P(S2_RK2(j)+h2*k1_S,C2_RK2(j)+h2*k1_C,P2_RK2(j)+h2*k1_P);

  S2_RK2(j+1) = S2_RK2(j) + (h2/2)*(k1_S + k2_S);
  C2_RK2(j+1) = C2_RK2(j) + (h2/2)*(k1_C + k2_C);
  P2_RK2(j+1) = P2_RK2(j) + (h2/2)*(k1_P + k2_P);
  E2_RK2(j+1) = E_null + C_null - C2_RK2(j+1);
endfor

% Step 2.7: RK4 Solution Loop 2: Small Time Step

for j = 1:1:(length(t2)-1)
  k1_S   = f_S(S2_RK4(j),C2_RK4(j),P2_RK4(j));
  k1_C   = f_C(S2_RK4(j),C2_RK4(j),P2_RK4(j));
  k1_P   = f_P(S2_RK4(j),C2_RK4(j),P2_RK4(j));

  k2_S   = f_S(S2_RK4(j)+(h2/2)*k1_S,C2_RK4(j)+(h2/2)*k1_C,P2_RK4(j)+(h2/2)*k1_P);
  k2_C   = f_C(S2_RK4(j)+(h2/2)*k1_S,C2_RK4(j)+(h2/2)*k1_C,P2_RK4(j)+(h2/2)*k1_P);
  k2_P   = f_P(S2_RK4(j)+(h2/2)*k1_S,C2_RK4(j)+(h2/2)*k1_C,P2_RK4(j)+(h2/2)*k1_P);

  k3_S   = f_S(S2_RK4(j)+(h2/2)*k2_S,C2_RK4(j)+(h2/2)*k2_C,P2_RK4(j)+(h2/2)*k2_P);
  k3_C   = f_C(S2_RK4(j)+(h2/2)*k2_S,C2_RK4(j)+(h2/2)*k2_C,P2_RK4(j)+(h2/2)*k2_P);
  k3_P   = f_P(S2_RK4(j)+(h2/2)*k2_S,C2_RK4(j)+(h2/2)*k2_C,P2_RK4(j)+(h2/2)*k2_P);

  k4_S   = f_S(S2_RK4(j)+h2*k3_S,C2_RK4(j)+h2*k3_C,P2_RK4(j)+h2*k3_P);
  k4_C   = f_C(S2_RK4(j)+h2*k3_S,C2_RK4(j)+h2*k3_C,P2_RK4(j)+h2*k3_P);
  k4_P   = f_P(S2_RK4(j)+h2*k3_S,C2_RK4(j)+h2*k3_C,P2_RK4(j)+h2*k3_P);

  S2_RK4(j+1) = S2_RK4(j) + (h2/6)*(k1_S + 2*k2_S + 2*k3_S + k4_S);
  C2_RK4(j+1) = C2_RK4(j) + (h2/6)*(k1_C + 2*k2_C + 2*k3_C + k4_C);
  P2_RK4(j+1) = P2_RK4(j) + (h2/6)*(k1_P + 2*k2_P + 2*k3_P + k4_P);
  E2_RK4(j+1) = E_null + C_null - C2_RK4(j+1);
endfor

% Step 2.8: IE Solution Loop 2: Small Time Step

F1    = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) S_new + h2*k_01*S_new*(E_null+C_null-C_new) - h2*k_02*C_new - S_prev;
F2    = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) -h2*k_01*S_new*(E_null+C_null-C_new) + (1+h2*k_02+h2*k_03)*C_new - h2*k_04*P_new*(E_null+C_null-C_new) - C_prev;
F3    = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) P_new - h2*k_03*C_new + h2*k_04*P_new*(E_null+C_null-C_new) - P_prev;
F_vec = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) [F1(S_new,C_new,P_new,S_prev,C_prev,P_prev);F2(S_new,C_new,P_new,S_prev,C_prev,P_prev);F3(S_new,C_new,P_new,S_prev,C_prev,P_prev)];

F1_dS = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) 1 + h2*k_01*(E_null+C_null-C_new);
F1_dC = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) -h2*k_01*S_new - h1_imp*k_02;
F1_dP = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) 0;
F2_dS = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) -h2*k_01*(E_null+C_null-C_new);
F2_dC = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) 1 + h2*k_01*S_new + h2*k_02 + h2*k_03 + h2*k_04*P_new;
F2_dP = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) -h2*k_04*(E_null+C_null-C_new);
F3_dS = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) 0;
F3_dC = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) -h2*k_03 - h1_imp*k_04*P_new;
F3_dP = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) 1 + h2*k_04*(E_null+C_null-C_new);

JF_mat  = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) [F1_dS(S_new,C_new,P_new,S_prev,C_prev,P_prev) F1_dC(S_new,C_new,P_new,S_prev,C_prev,P_prev) F1_dP(S_new,C_new,P_new,S_prev,C_prev,P_prev); F2_dS(S_new,C_new,P_new,S_prev,C_prev,P_prev) F2_dC(S_new,C_new,P_new,S_prev,C_prev,P_prev) F2_dP(S_new,C_new,P_new,S_prev,C_prev,P_prev); F3_dS(S_new,C_new,P_new,S_prev,C_prev,P_prev) F3_dC(S_new,C_new,P_new,S_prev,C_prev,P_prev) F3_dP(S_new,C_new,P_new,S_prev,C_prev,P_prev)];
X_curr  = [S2_IE(1);C2_IE(1);P2_IE(1)];

for l = 1:1:(length(t2)-1)
  for j = 1:1:50
    JF_step = JF_mat(X_curr(1),X_curr(2),X_curr(3),S2_IE(l),C2_IE(l),P2_IE(l));
    H       = -inv(JF_step)*F_vec(X_curr(1),X_curr(2),X_curr(3),S2_IE(l),C2_IE(l),P2_IE(l));
    X_curr  = X_curr + H;
  endfor
  S2_IE(l+1) = X_curr(1);
  C2_IE(l+1) = X_curr(2);
  P2_IE(l+1) = X_curr(3);
  E2_IE(l+1) = E_null+C_null-C2_IE(l+1);
endfor

% Step 2.9: BDF2 Solution Loop 2: Small Time Step

JF_mat  = @(S_new,C_new,P_new,S_prev,C_prev,P_prev) [F1_dS(S_new,C_new,P_new,S_prev,C_prev,P_prev) F1_dC(S_new,C_new,P_new,S_prev,C_prev,P_prev) F1_dP(S_new,C_new,P_new,S_prev,C_prev,P_prev); F2_dS(S_new,C_new,P_new,S_prev,C_prev,P_prev) F2_dC(S_new,C_new,P_new,S_prev,C_prev,P_prev) F2_dP(S_new,C_new,P_new,S_prev,C_prev,P_prev); F3_dS(S_new,C_new,P_new,S_prev,C_prev,P_prev) F3_dC(S_new,C_new,P_new,S_prev,C_prev,P_prev) F3_dP(S_new,C_new,P_new,S_prev,C_prev,P_prev)];
X_curr  = [S2_BDF2(1);C2_BDF2(1);P2_BDF2(1)];

for j = 1:1:100
  JF_step = JF_mat(X_curr(1),X_curr(2),X_curr(3),S2_BDF2(1),C2_BDF2(1),P2_BDF2(1));
  H       = -inv(JF_step)*F_vec(X_curr(1),X_curr(2),X_curr(3),S2_BDF2(1),C2_BDF2(1),P2_BDF2(1));
  X_curr  = X_curr + H;
endfor

S2_BDF2(2) = X_curr(1);
C2_BDF2(2) = X_curr(2);
P2_BDF2(2) = X_curr(3);
E2_BDF2(2) = E_null + C_null - C2_BDF2(2);

F1_BDF2    = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) -S_new + (4/3)*S_curr - (1/3)*S_old + (2*h2/3)*(-k_01*S_new*(E_null+C_null-C_new) + k_02*C_new);
F2_BDF2    = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) -C_new + (4/3)*C_curr - (1/3)*C_old + (2*h2/3)*(k_01*S_new*(E_null+C_null-C_new) - (k_02+k_03)*C_new + k_04*P_new*(E_null+C_null-C_new));
F3_BDF2    = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) -P_new + (4/3)*P_curr - (1/3)*P_old + (2*h2/3)*(k_03*C_new - k_04*P_new*(E_null+C_null-C_new));
F_vec_BDF2 = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) [F1_BDF2(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old);F2_BDF2(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old);F3_BDF2(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old)];

F1_BDF2_dS = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) -1 - (2/3)*h2*k_01*(E_null+C_null-C_new);
F1_BDF2_dC = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) (2/3)*h2*k_01*S_new + (2/3)*h2*k_02;
F1_BDF2_dP = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) 0;
F2_BDF2_dS = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) (2/3)*h2*k_01*(E_null+C_null-C_new);
F2_BDF2_dC = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) -1 - (2/3)*h2*k_01*S_new - (2/3)*h2*(k_02+k_03) - (2/3)*h1_imp*k_04*P_new;
F2_BDF2_dP = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) (2/3)*h2*k_04*(E_null+C_null-C_new);
F3_BDF2_dS = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) 0;
F3_BDF2_dC = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) (2/3)*h2*k_03 + (2/3)*h1_imp*k_04*P_new;
F3_BDF2_dP = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) -1 - (2/3)*h2*k_04*(E_null+C_null-C_new);

JF_mat_BDF2 = @(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) [F1_BDF2_dS(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) F1_BDF2_dC(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) F1_BDF2_dP(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old);F2_BDF2_dS(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) F2_BDF2_dC(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) F2_BDF2_dP(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old);F3_BDF2_dS(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) F3_BDF2_dC(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old) F3_BDF2_dP(S_new,C_new,P_new,S_curr,C_curr,P_curr,S_old,C_old,P_old)];

for l = 1:1:(length(t2)-2)
  X_curr  = [S2_BDF2(l+1);C2_BDF2(l+1);P2_BDF2(l+1)];
  for j = 1:1:100
    JF_step = JF_mat_BDF2(X_curr(1),X_curr(2),X_curr(3),S2_BDF2(l+1),C2_BDF2(l+1),P2_BDF2(l+1),S2_BDF2(l),C2_BDF2(l),P2_BDF2(l));
    H       = -inv(JF_step)*F_vec_BDF2(X_curr(1),X_curr(2),X_curr(3),S2_BDF2(l+1),C2_BDF2(l+1),P2_BDF2(l+1),S2_BDF2(l),C2_BDF2(l),P2_BDF2(l));
    X_curr  = X_curr + H;
  endfor
  S2_BDF2(l+2) = X_curr(1);
  C2_BDF2(l+2) = X_curr(2);
  P2_BDF2(l+2) = X_curr(3);
  E2_BDF2(l+2) = E_null+C_null-C2_BDF2(l+2);
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

% Step 4.1: Numerical Solutions of NSFD 1: Large Time Step

figure(1)
plot(t1_exp,S1,'linewidth',1,'linestyle','-','color','red')
hold on
plot(t1_exp,E1,'linewidth',1,'linestyle','-','color','blue')
hold on
plot(t1_exp,C1,'linewidth',1,'linestyle','-','color','magenta')
hold on
plot(t1_exp,P1,'linewidth',1,'linestyle','-','color','black')
hold on
plot(t1_exp,S_star_irr*ones(1,length(t1_exp)),'linewidth',1,'linestyle','--','color','red')
hold on
plot(t1_exp,E_star_irr*ones(1,length(t1_exp)),'linewidth',1,'linestyle','--','color','blue')
hold on
plot(t1_exp,C_star_irr*ones(1,length(t1_exp)),'linewidth',1,'linestyle','--','color','magenta')
hold on
plot(t1_exp,P_star_irr*ones(1,length(t1_exp)),'linewidth',1,'linestyle','--','color','black')
hold on
title('Concentrations of enzymatic reaction (k_{4} = 0.0, h = 20, NSFD)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'S(t)','E(t)','C(t)','P(t)','S^{*}','E^{*}','C^{*}','P^{*}'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 T1_exp])
hold off

% Step 4.2: Numerical Solutions of RK2 1: Large Time Step

figure(2)
plot(t1_exp,S1_RK2,'linewidth',1,'linestyle','-','color','red')
hold on
plot(t1_exp,E1_RK2,'linewidth',1,'linestyle','-','color','blue')
hold on
plot(t1_exp,C1_RK2,'linewidth',1,'linestyle','-','color','magenta')
hold on
plot(t1_exp,P1_RK2,'linewidth',1,'linestyle','-','color','black')
hold on
plot(t1_exp,S_star_irr*ones(1,length(t1_exp)),'linewidth',1,'linestyle','--','color','red')
hold on
plot(t1_exp,E_star_irr*ones(1,length(t1_exp)),'linewidth',1,'linestyle','--','color','blue')
hold on
plot(t1_exp,C_star_irr*ones(1,length(t1_exp)),'linewidth',1,'linestyle','--','color','magenta')
hold on
plot(t1_exp,P_star_irr*ones(1,length(t1_exp)),'linewidth',1,'linestyle','--','color','black')
hold on
title('Concentrations of enzymatic reaction (k_{4} = 0.0, h = 20, RK2)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'S(t)','E(t)','C(t)','P(t)','S^{*}','E^{*}','C^{*}','P^{*}'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 T1_exp])
hold off

% Step 4.3: Numerical Solutions of RK4 1: Large Time Step

figure(3)
plot(t1_exp,S1_RK4,'linewidth',1,'linestyle','-','color','red')
hold on
plot(t1_exp,E1_RK4,'linewidth',1,'linestyle','-','color','blue')
hold on
plot(t1_exp,C1_RK4,'linewidth',1,'linestyle','-','color','magenta')
hold on
plot(t1_exp,P1_RK4,'linewidth',1,'linestyle','-','color','black')
hold on
plot(t1_exp,S_star_irr*ones(1,length(t1_exp)),'linewidth',1,'linestyle','--','color','red')
hold on
plot(t1_exp,E_star_irr*ones(1,length(t1_exp)),'linewidth',1,'linestyle','--','color','blue')
hold on
plot(t1_exp,C_star_irr*ones(1,length(t1_exp)),'linewidth',1,'linestyle','--','color','magenta')
hold on
plot(t1_exp,P_star_irr*ones(1,length(t1_exp)),'linewidth',1,'linestyle','--','color','black')
hold on
title('Concentrations of enzymatic reaction (k_{4} = 0.0, h = 20, RK4)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'S(t)','E(t)','C(t)','P(t)','S^{*}','E^{*}','C^{*}','P^{*}'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 T1_exp])
hold off

% Step 4.4: Numerical Solution of NSFD 1: Large Time Step

figure(4)
plot(t1_imp,S1_imp,'linewidth',1,'linestyle','-','color','red')
hold on
plot(t1_imp,E1_imp,'linewidth',1,'linestyle','-','color','blue')
hold on
plot(t1_imp,C1_imp,'linewidth',1,'linestyle','-','color','magenta')
hold on
plot(t1_imp,P1_imp,'linewidth',1,'linestyle','-','color','black')
hold on
plot(t1_imp,S_star_irr*ones(1,length(t1_imp)),'linewidth',1,'linestyle','--','color','red')
hold on
plot(t1_imp,E_star_irr*ones(1,length(t1_imp)),'linewidth',1,'linestyle','--','color','blue')
hold on
plot(t1_imp,C_star_irr*ones(1,length(t1_imp)),'linewidth',1,'linestyle','--','color','magenta')
hold on
plot(t1_imp,P_star_irr*ones(1,length(t1_imp)),'linewidth',1,'linestyle','--','color','black')
hold on
title('Concentrations of enzymatic reaction (k_{4} = 0.0, h = 1000, NSFD)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'S(t)','E(t)','C(t)','P(t)','S^{*}','E^{*}','C^{*}','P^{*}'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 4000])
hold off

% Step 4.5: Numerical Solutions of IE 1: Large Time Step

figure(5)
plot(t1_imp,S1_IE,'linewidth',1,'linestyle','-','color','red')
hold on
plot(t1_imp,E1_IE,'linewidth',1,'linestyle','-','color','blue')
hold on
plot(t1_imp,C1_IE,'linewidth',1,'linestyle','-','color','magenta')
hold on
plot(t1_imp,P1_IE,'linewidth',1,'linestyle','-','color','black')
hold on
plot(t1_imp,S_star_irr*ones(1,length(t1_imp)),'linewidth',1,'linestyle','--','color','red')
hold on
plot(t1_imp,E_star_irr*ones(1,length(t1_imp)),'linewidth',1,'linestyle','--','color','blue')
hold on
plot(t1_imp,C_star_irr*ones(1,length(t1_imp)),'linewidth',1,'linestyle','--','color','magenta')
hold on
plot(t1_imp,P_star_irr*ones(1,length(t1_imp)),'linewidth',1,'linestyle','--','color','black')
hold on
title('Concentrations of enzymatic reaction (k_{4} = 0.0, h = 1000, IE)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'S(t)','E(t)','C(t)','P(t)','S^{*}','E^{*}','C^{*}','P^{*}'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 4000])
hold off

% Step 4.6: Numerical Solutions of BDF2: Large Time Step

figure(6)
plot(t1_imp,S1_BDF2,'linewidth',1,'linestyle','-','color','red')
hold on
plot(t1_imp,E1_BDF2,'linewidth',1,'linestyle','-','color','blue')
hold on
plot(t1_imp,C1_BDF2,'linewidth',1,'linestyle','-','color','magenta')
hold on
plot(t1_imp,P1_BDF2,'linewidth',1,'linestyle','-','color','black')
hold on
plot(t1_imp,S_star_irr*ones(1,length(t1_imp)),'linewidth',1,'linestyle','--','color','red')
hold on
plot(t1_imp,E_star_irr*ones(1,length(t1_imp)),'linewidth',1,'linestyle','--','color','blue')
hold on
plot(t1_imp,C_star_irr*ones(1,length(t1_imp)),'linewidth',1,'linestyle','--','color','magenta')
hold on
plot(t1_imp,P_star_irr*ones(1,length(t1_imp)),'linewidth',1,'linestyle','--','color','black')
hold on
title('Concentrations of enzymatic reaction (k_{4} = 0.0, h = 1000, BDF2)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'S(t)','E(t)','C(t)','P(t)','S^{*}','E^{*}','C^{*}','P^{*}'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 4000])
hold off

% Step 4.7: Conservation Laws of NSFD 1: Large Time Step

figure(7)
plot(t1_exp,E1+C1,'linewidth',1,'linestyle','-','color','black')
hold on
plot(t1_exp,S1+C1+P1,'linewidth',1,'linestyle','-','color','red')
hold on
title('Conservation laws of enzymatic reaction (k_{4} = 0.0, h = 20, NSFD)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'E(t) + C(t)', 'S(t) + C(t) + P(t)'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 T1_exp])
ylim([0 0.8])
hold off

% Step 4.8: Conservation Laws of RK2 1: Large Time Step

figure(8)
plot(t1_exp,E1_RK2+C1_RK2,'linewidth',1,'linestyle','-','color','black')
hold on
plot(t1_exp,S1_RK2+C1_RK2+P1_RK2,'linewidth',1,'linestyle','-','color','red')
hold on
title('Conservation laws of enzymatic reaction (k_{4} = 0.0, h = 20, RK2)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'E(t) + C(t)', 'S(t) + C(t) + P(t)'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 T1_exp])
ylim([0 0.8])
hold off

% Step 4.9: Conservation Laws of RK4 1: Large Time Step

figure(9)
plot(t1_exp,E1_RK4+C1_RK4,'linewidth',1,'linestyle','-','color','black')
hold on
plot(t1_exp,S1_RK4+C1_RK4+P1_RK4,'linewidth',1,'linestyle','-','color','red')
hold on
title('Conservation laws of enzymatic reaction (k_{4} = 0.0, h = 20, RK4)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'E(t) + C(t)', 'S(t) + C(t) + P(t)'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 T1_exp])
ylim([0 0.8])
hold off

% Step 4.10: Conservation Laws of NSFD 1: Large Time Step

figure(10)
plot(t1_imp,E1_imp+C1_imp,'linewidth',1,'linestyle','-','color','black')
hold on
plot(t1_imp,S1_imp+C1_imp+P1_imp,'linewidth',1,'linestyle','-','color','red')
hold on
title('Conservation laws of enzymatic reaction (k_{4} = 0.0, h = 1000, NSFD)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'E(t) + C(t)', 'S(t) + C(t) + P(t)'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 T1_imp])
ylim([0 0.8])
hold off

% Step 4.11: Conservation Laws of IE 1: Large Time Step

figure(11)
plot(t1_imp,E1_IE+C1_IE,'linewidth',1,'linestyle','-','color','black')
hold on
plot(t1_imp,S1_IE+C1_IE+P1_IE,'linewidth',1,'linestyle','-','color','red')
hold on
title('Conservation laws of enzymatic reaction (k_{4} = 0.0, h = 1000, IE)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'E(t) + C(t)', 'S(t) + C(t) + P(t)'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 T1_imp])
ylim([0 0.8])
hold off

% Step 4.12: Conservation Laws of BDF2 1: Large Time Step

figure(12)
plot(t1_imp,E1_BDF2+C1_BDF2,'linewidth',1,'linestyle','-','color','black')
hold on
plot(t1_imp,S1_BDF2+C1_BDF2+P1_BDF2,'linewidth',1,'linestyle','-','color','red')
hold on
title('Conservation laws of enzymatic reaction (k_{4} = 0.0, h = 1000, BDF2)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'E(t) + C(t)', 'S(t) + C(t) + P(t)'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 T1_imp])
ylim([0 0.8])
hold off

% Step 4.13: Error of Conservation Law of All Explicit Methods: Large Time Step

figure(13)
plot(t1_exp,abs(S1+C1+P1-S_null-C_null-P_null),'linewidth',1,'linestyle','-','color','red')
hold on
plot(t1_exp,abs(S1_RK2+C1_RK2+P1_RK2-S_null-C_null-P_null),'linewidth',1,'linestyle','-','color','blue')
hold on
plot(t1_exp,abs(S1_RK4+C1_RK4+P1_RK4-S_null-C_null-P_null),'linewidth',1,'linestyle','-','color','black')
hold on
title('Absolute Errors of Conservation (k_{4} = 0.0, h = 20, All Explicit Methods)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'NSFD', 'RK2', 'RK4'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Absolute Error','fontsize',12)
xlim([0 T1_exp])
hold off

% Step 4.14: Error of Conservation Law of All Implicit Methods: Large Time Step

figure(14)
plot(t1_imp,abs(S1_imp+C1_imp+P1_imp-S_null-C_null-P_null),'linewidth',1,'linestyle','-','color','red')
hold on
plot(t1_imp,abs(S1_IE+C1_IE+P1_IE-S_null-C_null-P_null),'linewidth',1,'linestyle','-','color','blue')
hold on
plot(t1_imp,abs(S1_BDF2+C1_BDF2+P1_BDF2-S_null-C_null-P_null),'linewidth',1,'linestyle','-','color','black')
hold on
title('Absolute Errors of Conservation (k_{4} = 0.0, h = 1000, All Implicit Methods)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'NSFD', 'IE','BDF2'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Absolute Error','fontsize',12)
xlim([0 T1_imp])
hold off

% Step 4.15: Numerical Solutions of NSFD 2: Small Time Step

figure(15)
plot(t2(1:100:end),S2(1:100:end),'linewidth',1,'linestyle','-','color','red')
hold on
plot(t2(1:100:end),E2(1:100:end),'linewidth',1,'linestyle','-','color','blue')
hold on
plot(t2(1:100:end),C2(1:100:end),'linewidth',1,'linestyle','-','color','magenta')
hold on
plot(t2(1:100:end),P2(1:100:end),'linewidth',1,'linestyle','-','color','black')
hold on
plot(t2(1:100:end),S_star_irr*ones(1,length(t2(1:100:end))),'linewidth',1,'linestyle','--','color','red')
hold on
plot(t2(1:100:end),E_star_irr*ones(1,length(t2(1:100:end))),'linewidth',1,'linestyle','--','color','blue')
hold on
plot(t2(1:100:end),C_star_irr*ones(1,length(t2(1:100:end))),'linewidth',1,'linestyle','--','color','magenta')
hold on
plot(t2(1:100:end),P_star_irr*ones(1,length(t2(1:100:end))),'linewidth',1,'linestyle','--','color','black')
hold on
title('Concentrations of enzymatic reaction (k_{4} = 0.0, h = 0.1, NSFD)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'S(t)','E(t)','C(t)','P(t)','S^{*}','E^{*}','C^{*}','P^{*}'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 T2])
hold off

% Step 4.16: Numerical Solutions of RK2 2: Small Time Step

figure(16)
plot(t2(1:100:end),S2_RK2(1:100:end),'linewidth',1,'linestyle','-','color','red')
hold on
plot(t2(1:100:end),E2_RK2(1:100:end),'linewidth',1,'linestyle','-','color','blue')
hold on
plot(t2(1:100:end),C2_RK2(1:100:end),'linewidth',1,'linestyle','-','color','magenta')
hold on
plot(t2(1:100:end),P2_RK2(1:100:end),'linewidth',1,'linestyle','-','color','black')
hold on
plot(t2(1:100:end),S_star_irr*ones(1,length(t2(1:100:end))),'linewidth',1,'linestyle','--','color','red')
hold on
plot(t2(1:100:end),E_star_irr*ones(1,length(t2(1:100:end))),'linewidth',1,'linestyle','--','color','blue')
hold on
plot(t2(1:100:end),C_star_irr*ones(1,length(t2(1:100:end))),'linewidth',1,'linestyle','--','color','magenta')
hold on
plot(t2(1:100:end),P_star_irr*ones(1,length(t2(1:100:end))),'linewidth',1,'linestyle','--','color','black')
hold on
title('Concentrations of enzymatic reaction (k_{4} = 0.0, h = 0.1, RK2)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'S(t)','E(t)','C(t)','P(t)','S^{*}','E^{*}','C^{*}','P^{*}'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 T2])
hold off

% Step 4.17: Numerical Solutions of RK4 2: Small Time Step

figure(17)
plot(t2(1:100:end),S2_RK4(1:100:end),'linewidth',1,'linestyle','-','color','red')
hold on
plot(t2(1:100:end),E2_RK4(1:100:end),'linewidth',1,'linestyle','-','color','blue')
hold on
plot(t2(1:100:end),C2_RK4(1:100:end),'linewidth',1,'linestyle','-','color','magenta')
hold on
plot(t2(1:100:end),P2_RK4(1:100:end),'linewidth',1,'linestyle','-','color','black')
hold on
plot(t2(1:100:end),S_star_irr*ones(1,length(t2(1:100:end))),'linewidth',1,'linestyle','--','color','red')
hold on
plot(t2(1:100:end),E_star_irr*ones(1,length(t2(1:100:end))),'linewidth',1,'linestyle','--','color','blue')
hold on
plot(t2(1:100:end),C_star_irr*ones(1,length(t2(1:100:end))),'linewidth',1,'linestyle','--','color','magenta')
hold on
plot(t2(1:100:end),P_star_irr*ones(1,length(t2(1:100:end))),'linewidth',1,'linestyle','--','color','black')
hold on
title('Concentrations of enzymatic reaction (k_{4} = 0.0, h = 0.1, RK4)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'S(t)','E(t)','C(t)','P(t)','S^{*}','E^{*}','C^{*}','P^{*}'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 T2])
hold off

% Step 4.18: Numerical Solutions of IE 2: Small Time Step

figure(18)
plot(t2(1:100:end),S2_IE(1:100:end),'linewidth',1,'linestyle','-','color','red')
hold on
plot(t2(1:100:end),E2_IE(1:100:end),'linewidth',1,'linestyle','-','color','blue')
hold on
plot(t2(1:100:end),C2_IE(1:100:end),'linewidth',1,'linestyle','-','color','magenta')
hold on
plot(t2(1:100:end),P2_IE(1:100:end),'linewidth',1,'linestyle','-','color','black')
hold on
plot(t2(1:100:end),S_star_irr*ones(1,length(t2(1:100:end))),'linewidth',1,'linestyle','--','color','red')
hold on
plot(t2(1:100:end),E_star_irr*ones(1,length(t2(1:100:end))),'linewidth',1,'linestyle','--','color','blue')
hold on
plot(t2(1:100:end),C_star_irr*ones(1,length(t2(1:100:end))),'linewidth',1,'linestyle','--','color','magenta')
hold on
plot(t2(1:100:end),P_star_irr*ones(1,length(t2(1:100:end))),'linewidth',1,'linestyle','--','color','black')
hold on
title('Concentrations of enzymatic reaction (k_{4} = 0.0, h = 0.1, IE)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'S(t)','E(t)','C(t)','P(t)','S^{*}','E^{*}','C^{*}','P^{*}'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 T2])
hold off

% Step 4.19: Numerical Solutions of IE 2: Small Time Step

figure(19)
plot(t2(1:100:end),S2_BDF2(1:100:end),'linewidth',1,'linestyle','-','color','red')
hold on
plot(t2(1:100:end),E2_BDF2(1:100:end),'linewidth',1,'linestyle','-','color','blue')
hold on
plot(t2(1:100:end),C2_BDF2(1:100:end),'linewidth',1,'linestyle','-','color','magenta')
hold on
plot(t2(1:100:end),P2_BDF2(1:100:end),'linewidth',1,'linestyle','-','color','black')
hold on
plot(t2(1:100:end),S_star_irr*ones(1,length(t2(1:100:end))),'linewidth',1,'linestyle','--','color','red')
hold on
plot(t2(1:100:end),E_star_irr*ones(1,length(t2(1:100:end))),'linewidth',1,'linestyle','--','color','blue')
hold on
plot(t2(1:100:end),C_star_irr*ones(1,length(t2(1:100:end))),'linewidth',1,'linestyle','--','color','magenta')
hold on
plot(t2(1:100:end),P_star_irr*ones(1,length(t2(1:100:end))),'linewidth',1,'linestyle','--','color','black')
hold on
title('Concentrations of enzymatic reaction (k_{4} = 0.0, h = 0.1, BDF2)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'S(t)','E(t)','C(t)','P(t)','S^{*}','E^{*}','C^{*}','P^{*}'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 T2])
hold off

% Step 4.20: Conservation Laws of NSFD 2: Small Time Step

figure(20)
plot(t2,E2+C2,'linewidth',1,'linestyle','-','color','black')
hold on
plot(t2,S2+C2+P2,'linewidth',1,'linestyle','-','color','red')
hold on
title('Conservation laws of enzymatic reaction (k_{4} = 0.0, h = 0.1, NSFD)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'E(t) + C(t)', 'S(t) + C(t) + P(t)'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 T2])
ylim([0 0.8])
hold off

% Step 4.21: Conservation Laws of RK2 2: Small Time Step

figure(21)
plot(t2,E2_RK2+C2_RK2,'linewidth',1,'linestyle','-','color','black')
hold on
plot(t2,S2_RK2+C2_RK2+P2_RK2,'linewidth',1,'linestyle','-','color','red')
hold on
title('Conservation laws of enzymatic reaction (k_{4} = 0.0, h = 0.1, RK2)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'E(t) + C(t)', 'S(t) + C(t) + P(t)'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 T2])
ylim([0 0.8])
hold off

% Step 4.22: Conservation Laws of RK4 2: Small Time Step

figure(22)
plot(t2,E2_RK4+C2_RK4,'linewidth',1,'linestyle','-','color','black')
hold on
plot(t2,S2_RK4+C2_RK4+P2_RK4,'linewidth',1,'linestyle','-','color','red')
hold on
title('Conservation laws of enzymatic reaction (k_{4} = 0.0, h = 0.1, RK4)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'E(t) + C(t)', 'S(t) + C(t) + P(t)'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 T2])
ylim([0 0.8])
hold off

% Step 4.23: Conservation Laws of IE 2: Small Time Step

figure(23)
plot(t2,E2_IE+C2_IE,'linewidth',1,'linestyle','-','color','black')
hold on
plot(t2,S2_IE+C2_IE+P2_IE,'linewidth',1,'linestyle','-','color','red')
hold on
title('Conservation laws of enzymatic reaction (k_{4} = 0.0, h = 0.1, IE)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'E(t) + C(t)', 'S(t) + C(t) + P(t)'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 T2])
ylim([0 0.8])
hold off

% Step 4.24: Conservation Laws of IE 2: Small Time Step

figure(24)
plot(t2,E2_BDF2+C2_BDF2,'linewidth',1,'linestyle','-','color','black')
hold on
plot(t2,S2_BDF2+C2_BDF2+P2_BDF2,'linewidth',1,'linestyle','-','color','red')
hold on
title('Conservation laws of enzymatic reaction (k_{4} = 0.0, h = 0.1, BDF2)','interpreter','tex','fontsize',14,'fontweight','normal')
legend({'E(t) + C(t)', 'S(t) + C(t) + P(t)'},'interpreter','tex','location','eastoutside','fontsize',12)
xlabel('Time','fontsize',12)
ylabel('Concentrations','fontsize',12)
xlim([0 T2])
ylim([0 0.8])
hold off
