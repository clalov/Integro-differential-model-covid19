function rhs = susceptible(I,C,phi)

c = sum(phi.*I);
rhs =@(t,S) (S(t)*c)./(I(t)+ S(t)+ C(t));

end

