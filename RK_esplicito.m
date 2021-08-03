% Runge-Kutta Esplicito

function [yy,tt,neval,N] = RK_esplicito(f,t0,tf,y0,h,tableau)

% Inizializzazione delle variabili 
t=t0;           % tempo iniziale
y=y0;           % valore iniziale
neval=0;        % numero di valutazioni di f

d=length(y0);   % dimensione del problema

[c,A,b,p]=feval(tableau); % tableau 
s=length(b);    % numero di tappe del metodo

N = ceil((tf-t0)/h); % nro. di sotto-intervalli di tempo in cui dividere [t0,tf]
tt = [t0:h:tf]; % vettore dei tempi

if length(tt)==N
    tt = [tt tf];
end

yy = zeros(d,N+1); % matrice della soluzione approssimata
yy(:,1)=y0;        % dato iniziale = prima colonna di yy


for j=1:N
    
    % Aggiusto il passo di tempo per l'ultima esecuzione del ciclo
    if j==N
        h = tf - tt(j);
    end
    
    k=zeros(s,d); % k: matrice che contiene k_1,...,k_s
    z=f(t,y);
    neval=neval +1;
    k(1,:)=z;     % k_1
    
    % k_i per i=2,...,s
    for i=2:s
        taux=t + c(i) * h;
        yaux=y + h * (A(i,:)*k); 
        z=f(taux,yaux);
        k(i,:)=z';
        neval=neval +1;
    end
    
    % Passo del metodo
    y=y + h*(k'*b)';  
    t=t+h;
    yy(:,j+1) = y;
end

% disp('---------------------')
% disp(['Runge-Kutta esplicito di ordine ', num2str(p)])
% disp(['Numero di valutazioni di f: ',num2str(neval)]);
end