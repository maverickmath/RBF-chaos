function [y,lia]=sinemap(lambda,N);
% y=sinemap(lambda,N)
% simulates N points of the sine-map 
%    y(i)=lambda*sin(y(i-1));
% initial condition used is pi/3
% suggestions
% for cubic nonlinearity: lambda=1.2*pi;
% for quadractic nonlinearity: lambda=pi;

y(1)=pi/3;
sum = 0;

for i=2:N
   y(i)=lambda*sin(y(i-1));
   J = lambda*cos(y(i));
   sum = sum + (log10(abs(eig(J)))/log10(2));
   %sum = sum + (log(abs(eig(J))));
   lia(i) = sum/i;
end;

%plot(y(1:499),y(2:500),'.');
