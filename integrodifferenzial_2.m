function [rhs] = integrodifferenzial_2(k,alpha)

I = sym('I%d',[k+1 1]); % symbolic vector [I1, ..., Ik+1] that corresponds to [I0, ..., Ik] in the notes
IC = 5;               % initial condition

J = zeros(k,k+1);   

for a=1:k  
    alpha_a = alpha(a);
    h_a = 1/(a - (a-1));
    
    J(a,a) = h_a;
    J(a,a+1) = alpha_a + h_a; 
end

rhs = J*I;

%[tt, i] =ode45(rhs, [0 5], 5);