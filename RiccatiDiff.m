function [ tout, E ] = RiccatiDiff( A,B,C, D, beta, gamma, T,Nt)
%Computation of the solution of
%the solution to the Riccati Differential Equation
%for the control problem:
%min J(u)=\frac{1}{2}[\int_0^T|u(t)|^2dt+\int_0^T|Cx(t)|^2dt].
%x_t+Ax=Bu,    t\in (0,T)
%x(0)=x_0.
%In particular, the Riccati Differential Equation 
%for such problem reads as:
% \begin{equation}\label{RDE_1}
% \begin{cases}
% {\mathcal{E}}_t=\frac12 {C}^*{C}-({\mathcal{E}}{A}+{A}^*{\mathcal{E}})-2{\mathcal{E}}{B}{B}^*{\mathcal{E}}\qquad \forall t\in (0, +\infty )\\
% {\mathcal{E}}(0)=0.\\
% \end{cases}
% \end{equation}
%We employ the linear representation
%See, for instance, "Contr{\^o}le optimal: th{\'e}orie \& applications"
%by Emmanuel Tr{\'e}lat
%Proposition 4.3.5 page 58.
%WARNING: Different notation
%between this code and the above book!!!

%We define the matrix
M=[A,B*transpose(B);(beta)*transpose(C)*C,-transpose(A)];

%We compute the resolvent (matrix exponential) of M.
tout=linspace(0,T,Nt);
%We define the tensor where
%we will store the resolvent (matrix exponential) of M
%at the time instances in "tout".
R = repmat(0, [Nt, 2*size(A,1), 2*size(A,1)]);
for i=1:Nt
   R(i,:,:)=expm(tout(i)*M); 
end
%We define blocks of R
R1 = repmat(0, [Nt, size(A,1), size(A,1)]);
R2 = repmat(0, [Nt, size(A,1), size(A,1)]);
R3 = repmat(0, [Nt, size(A,1), size(A,1)]);
R4 = repmat(0, [Nt, size(A,1), size(A,1)]);
for i=1:Nt
for j=1:size(A,1)
    for l=1:size(A,1)
   R1(i,j,l)=R(i,j,l);
   R2(i,j,l)=R(i,j,size(A,1)+l);
   R3(i,j,l)=R(i,size(A,1)+j,l);
   R4(i,j,l)=R(i,size(A,1)+j,size(A,1)+l);
    end
end
end

%We define the Riccati operator
E = repmat(0, [Nt, size(A,1), size(A,1)]);
for i=1:Nt
    E(i,:,:)=(squeeze(R3(i,:,:))+(gamma)*squeeze(R4(i,:,:))*transpose(D)*D)*inv(squeeze(R1(i,:,:))+(gamma)*squeeze(R2(i,:,:))*transpose(D)*D);
end



end

