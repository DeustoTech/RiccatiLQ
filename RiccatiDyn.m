function [ dEdt ] = RiccatiDyn( t,E,A,B,C, beta)
%Dynamics of Riccati Differential Equation
%function [ dEdt ] = RiccatiDynamics( t,E,A=[control dynamics matrix],B=[control matrix],
%C=[observation matrix for the runnning cost],
%beta=[parameter for the state].
%Control problem
%min J(u)=\frac{1}{2}[\int_0^T|u(t)|^2dt+\int_0^T|Cx(t)|^2dt].
%x_t+Ax=Bu,    t\in (0,T)
%x(0)=x_0.



%We define the Dynamics related to the Riccati Defferential Equation.

E = reshape(E, size(A,1),size(A,1)); %Convert from "n^2"-by-1 to "n"-by-"n"
dEdt = (beta)*transpose(C)*C-(E*A+transpose(A)*E)-E*B*transpose(B)*E; %Dynamics in matrix notation.
dEdt = dEdt(:); %Convert from "n"-by-"n" to "n^2"-by-1


end

