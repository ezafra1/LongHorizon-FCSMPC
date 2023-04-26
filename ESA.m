function[U_op]=ESA(U_bar_unc, H, u_min, u_max, p)

%% u_i formation
% u_i \in \V = \{ u_{\min},...,u_{\max} \} \subset \Z
% u \in \V^n_u \subset \Z^n_u

j = u_min;
k = 1;
while (j<u_max || j==u_max)
    u_i(k) = j;
    j = j+1;
    k = k+1;
end

%% Input Combinations
% \U_k = \V^n_u Np \subset \Z^p; p=n_u Np

% U_cand is a matrix with all possible U_k

U_cand = combn(u_i,p); 

%% Initialization of squared radius
Rs2 = 1e20; 

%%%%%%%%%%%%%%% Exhaustive Search Algorithm %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:(length(U_cand))
    U_k = U_cand(i,:)'; 
    Ed = norm(H*U_k - U_bar_unc); % Eucledian distance of Hu to U_bar_unc

    if (Ed < Rs2) % Check Optimum distance
        Rs2 = Ed; % Store Optimum distance
        U_op = U_k;  % Store Optimum solution
    end
end

%%%%%%%%%%%%%%% End of Exhaustive Search Algorithm %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end