function[UpsilonT,Gamma,lSTE,W_inv,H]=MPC_Matrices_l(A,B,C,Np,lambda)

n=length(A);
aux=size(B);
n_u=aux(2);
aux=size(C);
n_y=aux(1);


%\Upsilon = [CB            0      ...     0       0;
%            CAB          CB      ...     0       0;
%                    
%            CA^(Np-1)  CA^(Np-1) ... CAB CB];

Upsilon=zeros(n_y*Np,n_u*Np);

% Create last row   CA^(Np-1) CA^(Np-1) ... CAB CB
 j=1;
 row_i=n_y*Np-1;
 row_f=n_y*Np;
for k=1:n_u:n_u*Np     
    col_i=k;
    col_f=k+n_u-1;
    Upsilon(row_i:row_f,col_i:col_f)=C*A^(Np-j)*B; 
    j=j+1;
end

%rest
for j=n_u*(Np-1):-n_u:1
    for k=n_y*(Np-1):-n_y:1
        Upsilon(k-1:k+n_y-2,j-n_u+1:j)=Upsilon(k+n_y-1:k+2*n_y-2,j+1:j+n_u);
    end
end

%Gamma
Gamma=zeros(n_y*Np,n);
cont=1;
for k=1:n_y:n_y*Np
    row_i=k;
    row_f=k+n_y-1;
   
    Gamma(row_i:row_f,:)=C*A^(cont);
    cont=cont+1;
end

%S
Im=eye(n_u);
S=zeros(n_u*Np,n_u*Np);
%first row
S(1:n_u,1:n_u)=Im;
%rest
if (Np>1)
    for k=Np:-1:2
        S(n_u*(k-1)+1:n_u*k,n_u*(k-1)+1:n_u*k)=Im;
        S(n_u*(k-1)+1:n_u*k,n_u*(k-2)+1:n_u*(k-1))=-Im;
    end
end


% E
E=zeros(n_u*Np,n_u);
E(1:n_u,1:n_u)=Im;

% W
UpsilonT=Upsilon';
lSTE=lambda*S'*E;

W=(UpsilonT*Upsilon + lambda*(S'*S));
W_inv=inv(W);
Hinv=chol(W_inv);
H=inv(Hinv)';

end