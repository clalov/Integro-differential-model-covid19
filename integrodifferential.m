function [YY,D,T]=integrodifferential()

% This code aims to simulate the intero-differential model that follows:
%
% dS/dt         = - S(t)*f(t)/(S(t)+I(t)+C(t)) \int phi(x)*I(t,x) dx
% I(t,0)        = - S(t)*f(t)/(S(t)+I(t)+C(t)) \int phi(x)I(t,x)
% dI/dx + dI/dt = -(gammaI(x)+eta(x)+q(x))*I(t,x)
% I(t)          = \int I(t,x) dx
% Y(t,0)        = \int eta(x)*I(t,x) dx
% dY/dx + dY/dt = -h(x)*Y(t,x)
% Y(t)          = \int Y(t,x) dx
% H(t,0)        = \int h(x)*Y(t,x) dx
% dH/dx + dH/dt = -(gammaY(x) + d(x))*H(t,x)
% H(t)          = \int H(t,x) dx
% Q(t,0)        = \int q(x)*I(t,x) dx
% dQ/dx + dQ/dt = - gammaQ(x)*Q(t,x)
% Q(t)          = \int Q(t,x) dx
% C(t,0)        = \int gammaI(x)*I(t,x) + gammaQ(x)*Q(t,x) + gammaY(x)*H(t,x) dx
% C(t,x)        = C(t-x,0)
% C(t)          = \int C(t,x) dx
% D(t,0)        = \int d(x)*H(t,x) dx
% D(t,x)        = D(t-x,0)
% D(t)          = \int D(t,x) dx
%
% The model introduces the AGE OF INFECTION that stands for the time
% elapsed since the infection occoured. 
%
% We use the day as the basic unit for time and age:
%                   -> time l \in {0,1,2...,tf}
%                   -> age a \in {1,2,...,af}
%                                                  
% OUTPUT: YY = [I(t); Y(t); H(t); Q(t); C(t)] 
%               -> infected, symptomatic, hospitalized, quarantined,
%                  recovered at time t, for t in T
%         D(t)  -> deceased at time t, for t in T
%         T     -> time interval

% Initial values
% (population size N = 1000) 
S0 = 10.000; % S(0)
I0 = 1.000;  % I(0,1) 
Y0 = 0; % Y(0,1) 
C0 = 0;  % C(0,1)
H0 = 0;   % H(0,1)
Q0 = 0;  % Q(0,1)
D0 = 0;  % D(0,1)


tf=30;        % final time 
T = [0:1:tf]; % interval time
af=10;         

% Compartments variables initialization
S = zeros(1,tf+1);  % suceptibles
I = zeros(tf+1,af); % asymptomatic infected/infectious 
Y = I;              % symptomatics
C = S;              % recovered 
H = I;              % hospitalized or treated at home 
Q = I;              % quarantined  
D = S;              % deceased 
S(1) = S0;
I(1,1) = I0;
Y(1,1) = Y0;
C(1) = C0;
H(1,1) = H0;
Q(1,1) =  Q0;
D(1) = D0;

% Output matrix initialization
YY = zeros(tf+1,5);
       
% Vectors of parameters
% (random values within realistic intervals are taken;
% we assume them to be stepwise functions with costant values)
%
% recovery rates
gammaI = 0.1 + (0.3-0.1).*rand(1,af);       % 0.1429*ones(1,af); 
gammaQ = 0.05 + (0.1 - 0.05).*rand(1,af);   % 0.0714*ones(1,af); 
gammaY = 0.05 + (0.1 - 0.05).*rand(1,af);   % 0.0714*ones(1,af);
% rate of becoming symptomatic
eta =  0.02 + (0.06 - 0.02).*rand(1,af);    % 0.0333*ones(1,af);
% hospitalization rate of symptomatic people
h =  0.2 + (0.5 - 0.2).*rand(1,af);         % 0.2500*ones(1,af);
% death rate of hospitalized people
d =  0.1 + (0.2 - 0.1).*rand(1,af);         % 0.0333*ones(1,af);
% isolation rate of infected people
q = 0.07 + (0.1 - 0.07).*rand(1,af);        % 0.07732*ones(1,af);

% survival transition rates initialization 
sigmaI = zeros(1,af);
sigmaY = sigmaI; sigmaH = sigmaI; sigmaQ = sigmaI; phi = sigmaI;
% pphi = zeros(1,af);

for a=1:af 
    % survival transition rates
    sigmaI(a) = exp(-(gammaI(a) + eta(a) + q(a)));
    sigmaY(a) = exp(-h(a));
    sigmaH(a) = exp(-(gammaY(a) + d(a)));
    sigmaQ(a) = exp(-gammaQ(a));
    % rate of secondary transmissions per single infectious case
    n=7.5; beta=1/2;
    phi(a) = gammainc(n,beta*a)-gammainc(n, beta*(a-1))./gamma(af);    
    f=@(t) 1;
    % pphi(a)=@(t) phi(a)*f(t); 
end

% Matrices of survival transition rates
SigmaI = [diag(sigmaI(2:end)) zeros(1,af-1)']; % size (af-1)xaf
SigmaY = [diag(sigmaY(2:end)) zeros(1,af-1)'];
SigmaH = [diag(sigmaH(2:end)) zeros(1,af-1)'];
SigmaQ = [diag(sigmaQ(2:end)) zeros(1,af-1)'];

X = zeros(tf+1,4*af+1);
X(1,:) = [I(1,:), Y(1,:), H(1,:), Q(1,:), C(1)];

for l=1:tf+1
    x = X(l,:);
    U = [ones(1,af), zeros(1,3*af), 1]; 
    Phi = [phi, zeros(1,3*af+1)];
    Z = [zeros(1,2*af), d, zeros(1,af+1)];
    
    ff=f(l-1);
    
    % blocks of the matrix M
    A = [(sigmaI(1)*S(l)*ff./(S(l) + U*x'))*phi; SigmaI]; % (af)x(af)
    B = [sigmaY(1)*eta; zeros(af-1,af)];  % (af+1)x(af+1)
    DD = [zeros(1,af); SigmaY];
    E = [sigmaH(1)*h ;zeros(af-1,af)];
    F = [zeros(1,af); SigmaH];
    G = [sigmaQ(1)*q; zeros(af-1,af)];
    L = [zeros(1, af); SigmaQ];
    
    M = [A,         zeros(af),   zeros(af), zeros(af), zeros(af,1);
         B,         DD,          zeros(af), zeros(af), zeros(af,1);
         zeros(af), E,           F,         zeros(af), zeros(af,1);
         G,         zeros(af),   zeros(af), L,         zeros(af,1);
         gammaI,    zeros(1,af), gammaY,    gammaQ,    1];
    
    S(l+1) = S(l) - (S(l)*ff)./(S(l)+ U*x')*Phi*x';
    X(l+1,:)=M*x';
    D(l+1) = D(l) + Z*x';
    
    N = [ones(1,af),  zeros(1,af), zeros(1,af), zeros(1,af), 0;
         zeros(1,af), ones(1,af),  zeros(1,af), zeros(1,af), 0;
         zeros(1,af), zeros(1,af), ones(1,af),  zeros(1,af), 0;
         zeros(1,af), zeros(1,af), zeros(1,af), ones(1,af),  0;
         zeros(1,af), zeros(1,af), zeros(1,af), zeros(1,af), 1];
    YY(l,:) = N*x';
end




