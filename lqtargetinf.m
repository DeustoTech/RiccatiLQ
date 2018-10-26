function [ uopt, x] = lqtargetinf( A,B,C,beta, q, z, x0,T,Nt )
%We look for the optimal pair (uopt, xopt) for the optimal control problem:
%min
%J(u)=\frac{1}{2}[\int_0^{\infty}|u(t)-q(t)|^2dt+\beta\int_0^{+\infty}|C(x(t)-z(t))|^2dt,
%where
%x_t+Ax=Bu,    t\in (0,+\infty)
%x(0)=x_0
%and
%%z_t+Az=Bq,    t\in (0,+\infty).
%The solution is plotted in the time interval [0,T].
%The partition of the time interval in the discretization is:
%tout=linspace(0,T,Nt).
%We employ Riccati's theory. We follow:
%"Long time versus steady state optimal control",
%by professor Alessio Porretta and professor Enrique Zuazua,
%Lemma 2.6 page 4249.
%WARNING(1):
%both the control target and the state target are function of time,
%beta\geq 0,
%x0\in\mathbb{R}^{size(A,1)}
%T>0
%Nt\in\mathbb{N}.

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
if (beta<0)
    error('\beta must be a positive real number.');
end
classq=class(q);
if (classq(1,1)~='f')
    error('q must be a function of time.');
elseif(size(q(0))~=[size(B,2),1])
    error('q must be: q:[0,T]\longrightarrow M(size(B,2),1;\mathbb{R}).');
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

%STEP 1. We solve Algebraic Riccati Equation
%related to the control problem with zero targets.

[X,L,K] = care(-A,B,transpose(C)*C);

%STEP 2: We determine the optimal state x
%by solving a closed loop system given by Riccati's operator.

options=odeset('RelTol',1e-12,'AbsTol',1e-14);

%We start determining tildex = optimal state - state target.

tout=linspace(0,T,Nt);
lt=length(tout);
[toutode45, tildex] = ode45(@(s, tildex) -A*tildex+B*(-K*tildex), tout, x0-z(0),options);

%We determine the optimal state, by using the formula:
%x=tildex+z

x=zeros(lt,n);
for i=1:lt
x(i,:)=tildex(i)+transpose(z(i*(T/Nt)));
end

%STEP 3: We determine optimal control by using the optimal feedback law
%given by Riccati's operator.

uopt=zeros(lt,m);

for i=1:lt
    uopt(i,:)=-transpose(K*(transpose(tildex(i,:))))+transpose(q((i/lt)*T));
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