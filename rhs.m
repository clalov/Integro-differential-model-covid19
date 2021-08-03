function f = rhs(t,I)

k = length(I); 

% Dati
q = 0.07;
eta = 0.02;
gamma_I = 0.1;
alpha = -(gamma_I + eta + q);
beta = 1;
phi = 3;

J = zeros(k,k+1);  
for a = 1:k  
    alpha_a = alpha;
    h_a = 1/(a - (a-1));
    
    J(a,a) = h_a;
    J(a,a+1) = alpha_a + h_a; 
end

% supponendo che la sommatoria al denominatore della BC parta da a=1
I0 = (beta/gamma_I + sum(I(2:end)))*sum(phi*I(2:end));

i = [I0;I];

f = J*i;