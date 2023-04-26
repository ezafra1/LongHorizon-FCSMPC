function[U_op,Niter,Rs2_o]=SDAKbestl(U_bar_unc, H, u_min, u_max, p,Kbest) 

%% Input Parameter definition
%
% U_bar_unc is the Unconstrained solution in the transformed space
% (hypersphere center)
% H = Lattice generator matrix (Lower Triangular)
%
% x(k+1)=Ax(k)+Bu(k) 
% x\in \R^m; u\in \R^n_u
% 
% \V = \{ u_{\min},...,u_{\max} \} \subset \Z
% u \in \V^n_u \subset \Z^n_u
%
% \U_k = \V^n_u Np \subset \Z^p; p=n_u Np
%
% min |U_bar_unc - H*U_k|^2
%
% U_k = Input sequence vector;
% Ucand is a matrix that stores the K_best switching sequences with its cost
% U_op = The optimal input vector;
% Ucand = [Uk_1          Uk_2     ...    Uk_2Kbest;
%                     
%          cost(Uk_1)   cost(Uk_2)  ...  cost(Uk_2Kbest) ];

%K-best algorithm only searches forward
%In each layer, we store the K-best elements. From these K-best elements,
%we study 2*K-best possibilities in the lower layer.

Niter=0;
Ucand=zeros(p+1,2*Kbest);
Ucand_aux=zeros(p+1,2*Kbest);
U_op=zeros(p,1);

%%%%%%%%%%%%%%%%%% K-best Sphere Decoding  Algorithm %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic
i = 1; C=0;
% k is the index for the tree layer
% C is the index of input vector value

while (i <= p)
    
    minK = min(Kbest,2^(i-1));
    
    for j=1:2:minK*2                        % Calculate cost of the new 2 K-best candiates
        for C=u_min:1:u_max 
            Ucand(i,j+C) = C;
            Ucand(p+1,j+C) = (H(i,1:i)*Ucand(1:i,j+C) - U_bar_unc(i))^2 + Ucand(p+1,j+C);
            Niter=Niter+1;
        end
    end
    
    Ucand=sortrows(Ucand',p+1);             % Sort the 2 K-best candidates
    Ucand=Ucand';
    
    minK = min(Kbest,2^(i));                
    
    l=1; j=1;
    
    while(j<=minK*2)                        % Select K-best candidates
        if (Ucand(p+1,l)>0)
            Ucand_aux(:,j+0) = Ucand(:,l);
            Ucand_aux(:,j+1) = Ucand(:,l);
            j = j + 2;
        end
        
        l=l+1;
    end
    Ucand = Ucand_aux;

    i = i+1;
end % While (k <= p)

Rs2_o = Ucand(p+1,1);
U_op = Ucand(1:p,1);
%toc
%%%%%%%%%%%%%%%%%%%%%%%% End of K-best Sphere Decoder %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end