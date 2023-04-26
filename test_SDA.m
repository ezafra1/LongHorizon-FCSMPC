clear all
clc
%% Original FCS
% \V = \{ u_{\min},...,u_{\max} \} \subset \Z
u_min=0;
u_max=1;

%% Multistep FCS
% u \in \V^n_u \subset \Z^n_u
% \U_k = \V^n_u Np \subset \Z^p; p=n_u Np

n_u = 3;             % Number of inputs u=[Sa Sb Sc]'
Np = 5;              % Prediction Horizon length
p = n_u*Np;          % Length of switching sequences U_k
Kbest = 8;           % Kb parameter for K-best SDA
Niter_max = 32000;   % Maximum number of explored nodes for SDA


%% Basic Matrices
I_2 = eye(2);          % Identity matrices
I_3 = eye(3);
O_2 = zeros(2);        % Zero matrices
J = [0 -1;
   1  0];
Tr = sqrt(2/3)*[1     -0.5         -0.5;      % Clarke transform
              0   sqrt(3)/2   -sqrt(3)/2];
                 

%% Continuous-Time alpha-beta Model -> 2 lvl VSI with output LC filter as example

n_y = 2;               % Number of outputs for control y= [vo_alpha vo_beta]
Ts = 50e-6;          % Sample Time
Lf = 0.002;          % LC filter inductor
Cf = 50e-6;          % LC filter capacitor
w = 2*pi*50;         % Angular frequency (@50 Hz)
Vdc = 400;           % Dc link voltage
lambda = 60.0;       % Weighting factor: Switching effort penalty

Ac_xy = [O_2        -(1/Lf)*I_2         O_2     ;
      (1/Cf)*I_2      O_2       -(1/Cf)*I_2   ;
       O_2            O_2             w*J     ;];

Bc_xy = [(Vdc/Lf)*I_2 O_2 O_2]'*Tr ;
        
Cxy = [I_2 O_2 O_2;
      O_2 I_2 O_2];   % C for observer
 
ct_sys = ss(Ac_xy,Bc_xy,Cxy,[]);

%% Discrete-Time alpha-beta Model
dt_sys = c2d(ct_sys,Ts);
A = dt_sys.a;
B = dt_sys.b;


%% Multistep Matrices (H is lower triangular)
C = zeros(n_y,6);        % C for control, 2 rows because only voltage is controlled
C(1:n_y,3:4) = eye(2);

[UpsilonT,Gamma,lSTE,W_inv,H] = MPC_Matrices_l(A,B,C,Np,lambda); %Multistep matrices


%% Initialize variables
Uop_ESA = zeros(p,1); 
Uop_SDA = zeros(p,1);
Uop_Kbest = zeros(p,1);
Niter=0.0;          % Count nodes
Niter_ac=0.0;
Niter_avg=0.0;
sim = 0;            % Number of simulations
simMAX = 1;

while (sim < simMAX)

%% Unconstrained solution (randomly generated in this example)
U_unc = 4*(rand(p,1)-0.5);
U_unc = min(u_max+1, max(u_min-1, U_unc));
U_bar_unc = H*U_unc;          % Center of the sphere in the transformed space

%% Initial Vector U_ini -> Only for SDA
U_ini = zeros(p,1);               % Null vector
%U_vec_ini=round(U_vec_unc);    % Babai
% for k =1 : p
%    if (U_vec_ini(k,1)<u_min)
%        U_vec_ini(k,1)=u_min;
%    end
%    if (U_vec_ini(k,1)>u_max)
%        U_vec_ini(k,1)=u_max;
%    end
% end


%% Exhaustive Search Algorithm
Uop_ESA = ESA(U_bar_unc, H, u_min, u_max, p);

%% Sphere Decoding Algorithm
[Uop_SDA,Niter,Rs2] = SDA_l(U_bar_unc, H, u_min, u_max, p, U_ini,Niter_max);

%% K-best Sphere Decoding Algorithm
[Uop_Kbest,NiterKb,Rs2Kb] = SDAKbestl(U_bar_unc, H, u_min, u_max, p,Kbest);

Uop_ESA'
Uop_SDA'
Niter
Uop_Kbest'
NiterKb

Niter_ac=Niter_ac+Niter;
sim = sim+1;

end

Niter_avg=Niter_ac/sim;