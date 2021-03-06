%Monta diagrama de bifurca��o de modelos RBF para o mapa senoidal, variando 
%os par�metros de dispers�o das gaussianas [par�metro de bifurca��o]

close all;
clear;
clc;
warning off;

nu = 0;   %Definido
ny = 1;   %Definido
nul = 0;   %Definido
nyl = 1;   %Definido
np = 4;	%N�mero de centros [valor "padr�o": 4]
%dp = 1.4756553;	%Valor padr�o: 1,8 
%[Mas a concord�ncia matem�tica com os pontos fixos ocorre por volta de dp = 1.1,...]
%[... para o qual o comportamento n�o � ca�tico -- como resolver/tratar?]
lin = 1;  %Presen�a ou n�o do polin�mio linear
dist = 1; %Dist�ncia do ponto fixo [valor "padr�o": 1]

pbif = 50;
cg = 1;

[y,lia] = sinemap(1.2*pi,500); %Mapa senoidal com n�o-linearidade c�bica (pf's: 0; +- 2,439)

pf = [-2.4383; 2.4383]; %Ponto fixo a impor
c = [pf - dist; pf + dist]; %Centros sim�tricos em rela��o ao ponto fixo
cn = sort(c);

for dp = 2.0:.001:2.2, %201 pontos
   
   [Psi] = montaP([],y', nu, ny, nul, nyl, lin, cn, 500, np, dp); 

	Psi_n = Psi(:,2:end);

	%Estimativa de M�nimos Quadrados
	X0 = ((inv((Psi_n')*Psi_n))*(Psi_n'))*y(max([ny nyl])+1:end)';

	%Garante ponto fixo nulo
	A = [1 1 0 0 0; 0 0 1 1 0; 0 0 0 0 1];
	b1 = pf(1)*(1 - X0(5))/exp(-(dist^2)/(dp^2));
	b2 = pf(2)*(1 - X0(5))/exp(-(dist^2)/(dp^2));
	B = [b1; b2; X0(5)];
	X1 = mqermod(Psi_n, X0, A, B); %Novo vetor de par�metros

	const = 0;
	u1n = 0;
	u2n = 0;
	y1n = 1;
	y2n = 0;
	y3n = 0;
	sim = simulacao_errn(cn,length(cn),0,1,1000,[],[pi/3],X1,lin,dp,const,u1n,u2n,y1n,y2n,y3n);
   bifu(:,cg) = sim(end-pbif:end)';
   aux1(:,cg) = dp*ones(pbif+1,1);
	cg = cg + 1;
   
end;

figure;
hold on;
for aux2 = 1:size(bifu,2),
   yk = plot(aux1(:,aux2),bifu(:,aux2),'k.');
   set(yk,'markersize',4);
end
hold off;

axis([1.99 2.21 -4 4]);
set(gca,'FontSize',11,'LineWidth',1);
xlabel('Dispers�o das gaussianas');
ylabel('y(k)');