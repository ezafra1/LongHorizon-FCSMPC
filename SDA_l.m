function[U_op,Niter,Rs2]=SDA_l(U_bar_unc, H, u_min, u_max, p, U_ini,Niter_max) 

%% Input Parameter defination
%
% U_bar_unc is the Unconstrained solution in the transformed space
% (hypersphere center)
% H = Lattice generator matrix (Lower Triangular)
%
% x(k+1)=Ax(k)+Bu(k) 
% x\in \R^m; u\in \R^n_u
% 
% u_i \in \V = \{ u_{\min},...,u_{\max} \} \subset \Z
% u \in \V^n_u \subset \Z^n_u
%
% \U_k = \V^n_u Np \subset \Z^p; p=n_u Np
%
% min |Yc - H*U_k|^2
%
% U_ini is the initial vector from the lattice
% d = An one dimensional array to store the calculated distance; d=zeros(1,(p+1));
% U = Input vector; U=zeros(1,p);
% U_op = The optimal input vector; U_op=zeros(1,p);

Niter=0;
d=zeros((p+1),1);
U_k=zeros(p,1);
U_op=U_ini;

%% Initial Radius of Sphere (Squared value)
Rs2 = ceil((norm(U_bar_unc - H*U_ini))^2);
Rs2_ini = Rs2;


%%%%%%%%%%%%%%%%%% Sphere Decoding  Algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic
i = 0; C=0;
% k is the index of input vector dimension
% C is the index of input vector value

while (i >=0) &&(Niter < Niter_max) % Optimal or early termination
    i=i+1;
    j = u_min + C;
    
    if j <= u_max
        U_k(i) = j;
        d(i+1) = ((H(i,1:i)*U_k(1:i) - U_bar_unc(i))^2 + d(i));
        Niter=Niter+1;
        if d(i+1) > Rs2 % Outside the sphere
            i=i-1; C=C+1;
        else            % Inside the sphere
            if i==p           % If entire sequence 
               U_op = U_k;    % Update incumbent solution
               Rs2 = d(i+1);  % Save the current squared radius
                    
               i=i-1;
               C=C+1;
            else        % Keep exploring sequence
               C=0;                  
            end % if i==p         
        end % if d(i+1) > Rs2      
    else 
        i=i-2;
        if i>=0
            C=U_k(i+1) - u_min + 1;
        end %if i>=0
    end % j <= u_max  
end % While
%toc
%%%%%%%%%%%%%%%%%%%%%%%% End of Sphere Decoding  Algorithm %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end