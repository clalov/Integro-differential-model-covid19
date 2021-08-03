% Runge-Kutta Implicito

function [yy,tt,y,nit,neval_fun,neval_jac,nfail] = implicit_rk(f,Jac,t0,tf,y0,h,tableau,tol,nitmax)

[c,A,b] = feval(tableau);

d = length(y0); % dimensione del problema 

s = length(c); %s: step del metodo rk
N = fix((tf-t0)/h); %numero di sotto-intervalli in cui viene suddiviso [t0,tf]
%tt: griglia dei tempi
%tt = (t0:h:(t0+N*h));
tt = [t0:h:tf];
if length(tt)==N
    tt=[tt tf];
end
%yy: traiettoria
yy=zeros(d,N+1);
yy(:,1)=y0;
%nit: vettore che contiene il numero di iterazioni del metodo di Newton
%ad ogni passo
nit = zeros(1,N);
%nfail: numero di fallimenti del metodo di Newton
nfail = 0;
%neval_fun: numero di valutazione di f
neval_fun = 0;
%neval_jac: numero di valutazioni dello jacobiano
neval_jac = 0;

bigA = kron(A,eye(d));
bigbTrasp_per_bigInvA = kron(b', eye(d)) * kron(inv(A), eye(d));
F = zeros(d*s,1);
DF = zeros(d*s);

for n=1:N
    %aggiorno h del passo finale
    if n==N
        h=tf-tt(n);
    end
    %k: iterazioni del metodo di Newton
    k = 0;
    %Z: iterante iniziale del metodo di Newton
    Z = kron(ones(s,1),zeros(d,1));
    %Costruzione di F(Z)
    for i=1:s
        F(((i-1)*d+1) : (i*d)) = f(tt(n)+h*c(i), yy(:,n)+Z((i-1)*d+1 : (i*d)));
        neval_fun = neval_fun + 1;
    end
    %Definisco la funzione ausiliaria G di cui vado a cercare lo zero
    G = Z - h*(bigA * F);
    
    %Applico il metodo di Newton
    %---------------------------------------------------
    for j=1:s
        DF(:,((j-1)*d+1) : (j*d)) = kron(A(:,j),Jac(tt(n)+c(j)*h,yy(:,n)+Z(((j-1)*d+1) : (j*d))));
        neval_jac = neval_jac + 1;
    end
    
    %Calcolo l'incremento
    incr= -(eye(s*d)-h*DF)\G;
    
    %Calcolo la nuova iterata
    Z = Z + incr;
    %Incremento il contatore
    k = k + 1;
    
    while(norm(incr)>tol && k<nitmax)
        
        %Valuto la funzione G nella nuova iterata
        for i=1:s
            F((i-1)*d+1 : (i*d)) = f(tt(n)+h*c(i), yy(:,n)+Z((i-1)*d+1 : (i*d)));
            neval_fun = neval_fun + 1;
        end
        
        G = Z - h*(bigA * F);
        
        for j=1:s
            DF(:,((j-1)*d+1) : (j*d)) = kron(A(:,j),Jac(tt(n)+c(j)*h,yy(:,n)+Z(((j-1)*d+1) : (j*d))));
            neval_jac = neval_jac + 1;
        end
        
        %Aggiorno l'incremento
        incr= -(eye(s*d)-h*DF)\G;
        %Calcolo la nuova iterata
        Z = Z + incr;
        %Incremento il contatore
        k = k + 1;
    end
    %-----------------------------------------------------
    %Salvo il numero di iterate del metodo di newton al passo n
    nit(n) = k;
    %Calcolo l'approssimazione della traiettoria
    yy(:,n+1) = yy(:,n) + bigbTrasp_per_bigInvA * Z;
    
    %Aggiornamento del numero di fallimenti di Newton
    if k>=nitmax
        nfail = nfail + 1;
        disp(['Il metodo di Newton ha raggiunto il numero massimo di iterazioni in t = ',num2str(tt(n))]);
    end
end

y = yy(:,N+1);
end