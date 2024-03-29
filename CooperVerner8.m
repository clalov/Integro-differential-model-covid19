function [c,A,b,p]=CooperVerner8

% embedded RK pair Cooper-Verner order 8, from Butcher,  p. 180

r=sqrt(21);
c=[0 1/2 1/2 (7+r)/14 (7+r)/14 1/2 (7-r)/14 (7-r)/14 1/2 (7+r)/14 1]';
A=[zeros(1,11);
   1/2 zeros(1,10);
   1/4 1/4 zeros(1,9);
   1/7 -(7+3*r)/98 (21+5*r)/49 zeros(1,8);
   (11+r)/84 0 (18+4*r)/63 (21-r)/252 zeros(1,7);
   (5+r)/48 0 (9+r)/36 (-231+14*r)/360 (63-7*r)/80 zeros(1,6);
   (10-r)/42 0 (-432+92*r)/315 (633-145*r)/90 (-504+115*r)/70 (63-13*r)/35 zeros(1,5);
   1/14 0 0 0 (14-3*r)/126 (13-3*r)/63 1/9 zeros(1,4);
   1/32 0 0 0 (91-21*r)/576 11/72 -(385+75*r)/1152 (63+13*r)/128  zeros(1,3);
   1/14 0 0 0 1/9 -(733+147*r)/2205 (515+111*r)/504 -(51+11*r)/56 (132+28*r)/245 zeros(1,2);
   0 0 0 0 (-42+7*r)/18 (-18+28*r)/45 -(273+53*r)/72 (301+53*r)/72 (28-28*r)/45 (49-7*r)/18 0];
b=[1/20 0 0 0 0 0 0 49/180 16/45 49/180 1/20]';
p=8;