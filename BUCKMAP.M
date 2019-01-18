%function y=buckmap(y0,kapa,N);
function [y,lia]=buckmap(y0,kapa,N); %Possibilidade: determina expoente de Lyapunov

% y=buckmap(y0,kapa,N);
%
% simulates N points of the buck-map with derived in:
% Chaos from a buck switching regulator operating in 
% discontinuous mode, by C.K. Tse, Int. J. of Circuit
% Theory and Applications, vol. 22, 263-278, 1994.
%
% y0 is the initial condition. Suggested value y0=24.
% kapa (feedback factor) can be used as a bifurcation parameter
% suggested range: 0.1 < kapa < 0.26

% LAA, BH, 4/4/96

% Definition of constants
E=33; % (V) - Power supply
X=25; % (V) - Voltage reference value
tcr=0.12; % tcr=t/(C*(R+rc))
rtl=20;  % rtl=R*T/L

% Equation 16
% alpha=1-t/(C*(R+rc))+T^2/(2*C^2*(R+rc)^2);
alpha=1-tcr+tcr^2/2;

% Equation 17
% beta=R*T^2/(2*L*C*(R+rc));
beta=rtl*tcr/2;

% Equation 22
D=sqrt((1-alpha)*X^2/(beta*E*(E-X)));

% initial condition
y(1)=y0;

sum=0;

for i=2:N
%
  % Equation 20
  d=D-kapa*(y(i-1)-X);
%
  % Equation 21
  if d<0
   h=0;
  elseif d>1
   h=1;
  else
   h=d;
  end;

  % Equation 19
  y(i)=alpha*y(i-1)+(h^2*beta*E*(E-y(i-1)))/y(i-1);
  
  %Possibilidade: evolução do expoente de Lyapunov [ainda não funcionando]
  if y(i-1) > (D + kapa*X)/kapa,
     J = alpha;
  end;
  if y(i-1) < (-1 + D + kapa*X)/kapa,
     J = alpha + (beta*E*(2*y(i)-E))/y(i)^2;
  end;   
  if (y(i-1) <= (D + kapa*X)/kapa) | (y(i) >= (-1 + D + kapa*X)/kapa),
     J = alpha + ((beta*E*d)/y(i)^2)*((-2*kapa*(E-y(i))-d)*y(i-1)-d*(E-y(i)));
  end;
  
  %sum = sum + (log10(abs(eig(J)))/log10(2));
  sum = sum + (log(abs(eig(J))));
  lia(i) = sum/i;
   
end;

