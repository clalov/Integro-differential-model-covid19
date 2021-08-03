function [A,b] = integrodifferential_2(c1,c2)

k = length(c1);
A = zeros(k,k);   

A(1,1) = -(c1(1)+1);

for a = 2:k  
    A(a,a) = -(c1(a) + 1);
    A(a,a-1) = 1; 
end

b = @(t) [c2(t); zeros(k-1,1)];

