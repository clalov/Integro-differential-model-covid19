% MAIN integro-differential_2

T = 10;     % final time
k = 5;      % final age
hh = 1;     % time step 


% Vectors of parameters
% (random values within realistic intervals are taken;
% we assume them to be stepwise functions with costant values)
%
% recovery rates
gammaI = 0.1 + (0.3-0.1).*rand(1,k);       % 0.1429*ones(1,af); 
gammaQ = 0.05 + (0.1 - 0.05).*rand(1,k);   % 0.0714*ones(1,af); 
gammaY = 0.05 + (0.1 - 0.05).*rand(1,k);   % 0.0714*ones(1,af);
% rate of becoming symptomatic
eta =  0.02 + (0.06 - 0.02).*rand(1,k);    % 0.0333*ones(1,af);
% hospitalization rate of symptomatic people
h =  0.2 + (0.5 - 0.2).*rand(1,k);         % 0.2500*ones(1,af);
% death rate of hospitalized people
d =  0.1 + (0.2 - 0.1).*rand(1,k);         % 0.0333*ones(1,af);
% isolation rate of infected people
q = 0.07 + (0.1 - 0.07).*rand(1,k);        % 0.07732*ones(1,af);

% --------------------------------------------------------------
% Susceptible
i = randi(k,1,20);
c = randi(k,1,20);
phi = randi(k,1,20);
rhs = susceptibles(i,c,phi); % MODIFICA 
[S,tt] = RK_esplicito(@rhs,0,T,S0,hh,@CooperVerner8);

% Infectious 
c1 = gammaI + h + q;
c2 = @(t) 46; % MODIFICA
[A,b] = integrodifferential_2(c1,c2);
rhs = @(t,y) A*y(t) + b(t);
[I,t1] = RK_esplicito(@rhs,0,T,I0,hh,@CooperVerner8);
% [t1,I]=ode45(@rhs,[0 T],I0);
% I =@(t) sum(I);

% Symptomatic
c1 = h;
c2 =@(t) 433; % MODIFICA
[A,b] = integrodifferential_2(c1,c2);
rhs = @(t,y) A*y(t) + b(t);
[Y,t2] = RK_esplicito(@rhs,0,T,Y0,hh,@CooperVerner8);
% [t2,Y]=ode45(@rhs,[0 T],Y0);
% Y =@(t) sum(Y);

% Hospedalized
c1 = gammaY + d;
c2 =@(t) 65; % MODIFICA
[A,b] = integrodifferential_2(c1,c2);
rhs = @(t,y) A*y(t) + b(t);
[H,t3] = RK_esplicito(@rhs,0,T,H0,hh,@CooperVerner8);
% [t3,H]=ode45(@rhs,[0 T],H0);
% H =@(t) sum(H);

% Quarantined
c1 = gammaQ;
c2 =@(t) 214; % MODIFICA
[A,b] = integrodifferential_2(c1,c2);
rhs = @(t,y) A*y(t) + b(t);
[Q,t4] = RK_esplicito(@rhs,0,T,Q0,hh,@CooperVerner8);
% [t4,Q]=ode45(@rhs,[0 T],Q0);
% Q =@(t) sum(Q);

figure;
plot(t1,I,t2,Y,t3,H,t4,Q);
legend('I','Y','H','Q');
