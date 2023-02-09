syms x
A=-2; B=2; C=2; eta=1/2; N=1; Q=2; 

a = -eta/2*(B+Q*x)^2;
b = A*x+N*x^2;
c = 1/2*C^2*x^2;


dEx1 = (-b - sqrt(b^2-4*a*c))/(2*a);
dEx2 = (-b + sqrt(b^2-4*a*c))/(2*a);


Ex1 = int(dEx1,x);
Ex2 = int(dEx2,x);

% hold on
% fplot(Ex1,[0,6])
% fplot(Ex2,[-6,0])