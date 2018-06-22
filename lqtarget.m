function [ uopt, x] = lqtarget( A,B,C,D,beta, gamma, q, z, x0,T,Nt )
%We look for the optimal pair (uopt, xopt) for the optimal control problem:
%min
%J(u)=\frac{1}{2}[\int_0^T|u(t)-q(t)|^2dt+\beta\int_0^T|C(x(t)-z(t))|^2dt+\gamma\|D(x(T)-z(T))\|^2],
%where
%x_t+Ax=Bu,    t\in (0,T)
%x(0)=x_0.
%The partition of the time interval in the discretization is:
%tout=linspace(0,T,Nt).
%We employ Riccati's theory. We follow the approach of:
%"Contr{\^o}le optimal: th{\'e}orie \& applications",
%by professor Emmanuel Tr{\'e}lat, Proposition 4.4.1 page 60.
%WARNING(1):
%both the control target and the state target are function of time,
%beta\geq 0,
%gamma\geq 0,
%x0\in\mathbb{R}^{size(A,1)}
%T>0
%Nt\in\mathbb{N}.
%WARNING(2):
%Different notation between this code and the aforementioned book.

%Compatibility conditions.

if (size(A,1)~=size(A,2))
    error('the matrix A must be a square matrix.');
end
if (size(A,1)~=size(B,1))
    error('the matrix B must have the same number of rows of A.');
end
if (size(C,2)~=size(A,1))
    error('the matrix C must have the same number of coloumns of A.');
end
if (size(D,2)~=size(A,1))
    error('the matrix D must have the same number of coloumns of A.');
end
if (beta<0)
    error('\beta must be a positive real number.');
end
if (gamma<0)
    error('\gamma must be a positive real number.');
end
classq=class(q);
if (classq(1,1)~='f')
    error('q must be a function of time.');
elseif(size(q(0))~=[size(B,2),1])
    error('q must be: q:[0,T]\longrightarrow M(size(B,2),1;\mathbb{R}).');
end
classz=class(z);
if (classz(1,1)~='f')
    error('z must be a function of time.');
elseif(size(z(0))~=[size(A,1),1])
    error('z must be: z:[0,T]\longrightarrow M(size(A,1),1;\mathbb{R}).');
end
if (size(x0,1)~=size(A,1) | size(x0,2)~=1)
    error('size(x0) must be equal to (size(A,1),1).');
end
if (T<=0)
    error('T must be a strictly positive real number.');
end
if (floor(Nt)~=Nt | Nt<=2)
    error('Nt must be a natural number >= 3.');
    %WARNING:
    %If Nt=2, then "ode45" confuses the time parttion "tout"
    %with the time interval where the solution should be computed
    %[t_{initial},t_{final}].
end


% Size of state vector
n = size(A,1);
% Size of control vector 
m = size(B,2);

%STEP 1. We solve Riccati Differential Equation
%related to the control problem with zero targets.

[ tout, E ] = RiccatiDiff( A,B,C, D, beta, gamma,T,Nt );

%STEP 2. Determination of the remainder function h^T.

%We firstly compute \frac{d}{dt}z.
lt=length(tout);
zvect=zeros(lt,n);
for i=1:lt
   zvect(i,:)=z((i/lt)*T);
end

zdervect=zeros(lt,n);
zdervect(1,:)=(zvect(2,:)-zvect(1,:))/(T/lt);
for i=2:lt
    zdervect(i,:)=(zvect(i,:)-zvect(i-1,:))/(T/lt);
end

%We determine \eta(t)=h^T(T-t).
[toutode45, eta] = ode45(@(s, eta) -(transpose(A)*eta+(squeeze(interp1(tout, E, s)))*(B*transpose(B)*eta+(A*z(T-s)+transpose(interp1(tout,zdervect , T-s))-B*q(T-s)))), tout, zeros(n,1));

%we determine h^T,
%by the transformation t--->T-t.
h=zeros(size(eta));
for i=1:lt
       h(i,:)=eta(lt-i+1,:);
end

%STEP 3: We determine the optimal state x
%by solving a closed loop system given by Riccati's operator.

options=odeset('RelTol',1e-12,'AbsTol',1e-14);
[toutode45, x] = ode45(@(s, x) -A*x+B*(-transpose(B)*(squeeze(interp1(tout, E, T-s)))*(x-z(s))-transpose(B)*(transpose(interp1(tout, h, s)))+q(s)), tout, x0,options);

%STEP 4: We determine optimal control by using the optimal feedback law
%given by Riccati's operator.

uopt=zeros(lt,m);

for i=1:lt
    uopt(i,:)=-transpose(B)*(squeeze(E(lt-i+1,:,:))*(transpose(x(i,:))-z(((i-1)/lt)*T)))-transpose(B)*(transpose(h(i,:)))+q((i/lt)*T);
end

%We plot the optimal pair
%(optimal control, optimal state)
%in case m <= 2 and n>=2.

%We construct the vectors
%ct and st containing the values
%resp. of the state target and the control target
%corresponding to the time instances in the partion tout.

ct=zeros(Nt,m);
st=zeros(Nt,n); 
for i=1:lt
    ct(i,:)=transpose(q(tout(i)));
    st(i,:)=transpose(z(tout(i)));
end                

if (n <= 2 & m <= 2)

switch m
    case 1
        switch n
            case 1
                figure(1)
                clf(1);
                plot(tout, uopt, 'b', 'LineWidth', 2);
                hold on;
                plot(tout,ct(:,1),'k', 'LineWidth', 1,'LineStyle','--');
                title("optimal control for LQ")
                xlabel("t")
                ylabel("u")
                legend("optimal control", "control target",'Location','southeast');
                set(gca,'FontSize',26)
                figure(2)
                clf(2);
                plot(tout, x(:,1), 'r', 'LineWidth', 2);
                hold on;
                plot(tout,st,'k', 'LineWidth', 1,'LineStyle','--');
                title("optimal state for LQ")
                xlabel("t")
                ylabel("x")
                legend("optimal state", "state target",'Location','southeast');
                set(gca,'FontSize',26)
            case 2
                figure(1)
                clf(1);
                plot(tout, uopt, 'b', 'LineWidth', 2);
                hold on;
                plot(tout,ct(:,1),'k', 'LineWidth', 1,'LineStyle','--');
                title("optimal control for LQ")
                xlabel("t")
                ylabel("u")
                legend("optimal control", "control target",'Location','southeast');
                set(gca,'FontSize',26)
    
                figure(2)
                clf(2);
                plot(tout, x(:,1), 'r', 'LineWidth', 2);
                hold on;
                plot(tout,st(:,1),'k', 'LineWidth', 1,'LineStyle','--')
                title("1^{st} component of optimal state for LQ")
                xlabel("t")
                ylabel("x_1")
                legend("1^{st} component of optimal state", "1^{st} component of state target",'Location','southeast');
                set(gca,'FontSize',26)
    
                figure(3)
                clf(3);
                plot(tout, x(:,2), 'r', 'LineWidth', 2);
                hold on;
                plot(tout,st(:,2),'k', 'LineWidth', 1,'LineStyle','--');
                title("2^{nd} component of optimal state for LQ")
                xlabel("t")
                ylabel("x_2")
                legend("2^{nd} component of optimal state", "2^{nd} component of state target",'Location','southeast');
                set(gca,'FontSize',26)
    end

    case 2
        
        switch n
            case 1
                figure(1)
                clf(1);
                plot(tout, uopt(:,1), 'b', 'LineWidth', 2);
                hold on;
                plot(tout,ct(:,1),'k', 'LineWidth', 1,'LineStyle','--');
                title("1^{st} component of optimal control for LQ")
                xlabel("t")
                ylabel("u_1")
                legend("1^{st} component of optimal control", "1^{st} component of control target",'Location','southeast');
                set(gca,'FontSize',26)
   
                figure(2)
                clf(2);
                plot(tout, uopt(:,2), 'b', 'LineWidth', 2);
                hold on;
                plot(tout,ct(:,2),'k', 'LineWidth', 1,'LineStyle','--');
                title("2^{nd} of optimal control for LQ")
                xlabel("t")
                ylabel("u_2")
                legend("2^{nd} component of optimal control", "2^{nd} component of control target",'Location','southeast');
                set(gca,'FontSize',26)
    
                figure(3)
                clf(3);
                plot(tout, x(:,1), 'r', 'LineWidth', 2);
                hold on;
                plot(tout,st(:,1),'k', 'LineWidth', 1,'LineStyle','--');
                title("optimal state for LQ")
                xlabel("t")
                ylabel("x")
                plot(tout,st,'k', 'LineWidth', 1,'LineStyle','--');
                legend("optimal state", "state target",'Location','southeast');
                set(gca,'FontSize',26)
            case 2
                figure(1)
                clf(1);
                plot(tout, uopt(:,1), 'b', 'LineWidth', 2);
                hold on;
                plot(tout,ct(:,1),'k', 'LineWidth', 1,'LineStyle','--');
                title("1^{st} component of optimal control for LQ")
                xlabel("t")
                ylabel("u_1")
                legend("1^{st} component of optimal control", "1^{st} component of control target",'Location','southeast');
                set(gca,'FontSize',26)
   
                figure(2)
                clf(2);
                plot(tout, uopt(:,2), 'b', 'LineWidth', 2);
                hold on;
                plot(tout,ct(:,2),'k', 'LineWidth', 1,'LineStyle','--');
                title("2^{nd} component of optimal control for LQ")
                xlabel("t")
                ylabel("u_2")
                legend("2^{nd} component of optimal control", "2^{nd} component of control target",'Location','southeast');
                set(gca,'FontSize',26)
                
                figure(3)
                clf(3);
                plot(tout, x(:,1), 'r', 'LineWidth', 2);
                hold on;
                plot(tout,st(:,1),'k', 'LineWidth', 1,'LineStyle','--')
                title("1^{st} component of optimal state for LQ")
                xlabel("t")
                ylabel("x_1")
                legend("1^{st} component of optimal state", "1^{st} component of state target",'Location','southeast');
                set(gca,'FontSize',26)
    
                figure(4)
                clf(4);
                plot(tout, x(:,2), 'r', 'LineWidth', 2);
                hold on;
                plot(tout,st(:,2),'k', 'LineWidth', 1,'LineStyle','--');
                title("2^{nd} component of optimal state for LQ")
                xlabel("t")
                ylabel("x_2")
                legend("2^{nd} component of optimal state", "2^{nd} component of state target",'Location','southeast');
                set(gca,'FontSize',26)
        end

    end

end




end



