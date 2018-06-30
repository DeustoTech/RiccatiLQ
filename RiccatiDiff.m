function [ tout, E ] = RiccatiDiff( A,B,C, D, beta, gamma, T,Nt)
%Computation of the solution to the Riccati Differential Equation
%for the control problem:
%min J(u)=\frac{1}{2}[\int_0^T|u(t)|^2dt+\int_0^T|Cx(t)|^2dt].
%x_t+Ax=Bu,    t\in (0,T)
%x(0)=x_0.
%In particular, the Riccati Differential Equation 
%for such problem reads as:
% \begin{cases}
% {\mathcal{E}}_t=\frac12 {C}^*{C}-({\mathcal{E}}{A}+{A}^*{\mathcal{E}})-2{\mathcal{E}}{B}{B}^*{\mathcal{E}}\qquad \forall t\in (0, +\infty )\\
% {\mathcal{E}}(0)=0.\\
% \end{cases}
%The partition of the time interval in the discretization is:
%tout=linspace(0,T,Nt).
%For Riccati Differential Equation see, for instance,
%"Contr{\^o}le optimal: th{\'e}orie \& applications"
%by professor Emmanuel Tr{\'e}lat, section 4.3.

%STEP 1: We solve the Riccati Differential Equation
%in a vector notation, by employing "ode113".
%We impose strict restrictions
%both on the relative erro and absolute error.

options=odeset('RelTol',1e-12,'AbsTol',1e-14);
INIT=(0.5)*(gamma)*transpose(D)*D;
INITvect = INIT(:); %Convert from "n"-by-"n" to "n^2"-by-1
tout=linspace(0,T,Nt);
[tout, Evect] = ode113(@(t,Evect)RiccatiDyn(t, Evect, A,B,C, beta), tout, INITvect,options);

%STEP 2: We reshape "Evect" in the matrix form "E".

netg=length(tout);
E = repmat(0, [netg, size(A,1), size(A,1)]);
for i=1:length(tout)
E(i,:,:) = reshape(Evect(i,:), size(A)); %Convert from "n^2"-by-1 to "n"-by-"n"
end



end

