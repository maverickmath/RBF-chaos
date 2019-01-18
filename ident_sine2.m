%Terceira parte: garantir a localização dos pontos fixos simétricos -- 4 centros
%Variando a dispersão

close all;
clear;
clc;
warning off;

nu = 0;   %Definido
ny = 1;   %Definido
nul = 0;   %Definido
nyl = 1;   %Definido
np = 4;	%Número de centros [valor "padrão": 4]
dp = .5;	%Valor padrão: 1,8 
%[Mas a concordância matemática com os pontos fixos ocorre por volta de dp = 1.1,...]
%[... para o qual o comportamento não é caótico -- como resolver/tratar?]
lin = 1;  %Presença ou não do polinômio linear
dist = 1; %Distância do ponto fixo [valor "padrão": 1]

[y,lia] = sinemap(1.2*pi,10000); %Mapa senoidal com não-linearidade cúbica (pf's: 0; +- 2,439)

pf = [-2.4383; 2.4383]; %Ponto fixo a impor
c = [pf - dist; pf + dist]; %Centros simétricos em relação ao ponto fixo
cn = sort(c);

%Adicionar ruído [futuramente]

[Psi] = montaP([],y(1:500)', nu, ny, nul, nyl, lin, cn, 500, np, dp); 

Psi_n = Psi(:,2:end);

%Estimativa de Mínimos Quadrados
X0 = ((inv((Psi_n')*Psi_n))*(Psi_n'))*y(max([ny nyl])+1:500)';

%Garante ponto fixo nulo
A = [1 1 0 0 0; 0 0 1 1 0; 0 0 0 0 1];
b1 = pf(1)*(1 - X0(5))/exp(-(dist^2)/(dp^2));
b2 = pf(2)*(1 - X0(5))/exp(-(dist^2)/(dp^2));
B = [b1; b2; X0(5)];
X1 = mqermod(Psi_n, X0, A, B); %Novo vetor de parâmetros

%Determina os pontos fixos dos modelos
prox = [-2.3 0 2.3];
for jj = 1:3,
   pfs(jj) = fzero('modelo_seno1',prox(jj),[],X0,cn,dp);
   pfs_e(jj) = X0(5) - (2*X0(1)/dp^2)*(pfs(jj)-cn(1))*exp(-(pfs(jj)-cn(1))^2/dp^2)-(2*X0(2)/dp^2)*(pfs(jj)-cn(2))*exp(-(pfs(jj)-cn(2))^2/dp^2) - (2*X0(3)/dp^2)*(pfs(jj)-cn(3))*exp(-(pfs(jj)-cn(3))^2/dp^2) - (2*X0(4)/dp^2)*(pfs(jj)-cn(4))*exp(-(pfs(jj)-cn(4))^2/dp^2);
   pfm(jj) = fzero('modelo_seno1',prox(jj),[],X1,cn,dp);
   pfm_e(jj) = X1(5) - (2*X1(1)/dp^2)*(pfm(jj)-cn(1))*exp(-(pfm(jj)-cn(1))^2/dp^2)-(2*X1(2)/dp^2)*(pfm(jj)-cn(2))*exp(-(pfm(jj)-cn(2))^2/dp^2) - (2*X1(3)/dp^2)*(pfm(jj)-cn(3))*exp(-(pfm(jj)-cn(3))^2/dp^2) - (2*X1(4)/dp^2)*(pfm(jj)-cn(4))*exp(-(pfm(jj)-cn(4))^2/dp^2);
end;

%Calcula os erros na localização dos pontos fixos
pfo = [-2.4383 0 2.4383];
erropfs = [pfo - pfs].^2;
erropfm = [pfo - pfm].^2;
ssepfs = sum(erropfs);
ssepfs_pf = erropfs(1) + erropfs(3); 
ssepfm = sum(erropfm);
ssepfm_pf = erropfm(1) + erropfm(3); 


%Calcula contribuição extra ["violação das restrições"?]
ex = exp(-(pfo(1)-cn(3))^2/dp^2) + exp(-(pfo(1)-cn(4))^2/dp^2);
in = 2*exp(-(pfo(1)-cn(1))^2/dp^2);
ex_p = (ex*100)/in;

%Simulação

const = 0;
u1n = 0;
u2n = 0;
y1n = 1;
y2n = 0;
y3n = 0;
sim = simulacao_errn(cn,length(cn),0,1,10000,[],[pi/3],X0,lin,dp,const,u1n,u2n,y1n,y2n,y3n);
reta = -4:.1:4;
yk = plot(y(1:499),y(2:500),'.'); 
hold on;
yl = plot(sim(1:499),sim(2:500),'ro'); 
hold on;
plot(reta,reta,'b');
set(gca,'FontSize',11,'LineWidth',1);
set(gcf,'Name','Modelo apenas com centros simétricos','NumberTitle','off');
xlabel('y(k)');
ylabel('y(k+1)');
set(yk,'markersize',5);
set(yl,'markersize',5);
axis([-4 4 -4 4]);

errods = [y(1:500) - sim(1:500)].^2;
sseds = sum(errods);

sim1 = simulacao_errn(cn,length(cn),0,1,10000,[],[pi/3],X1,lin,dp,const,u1n,u2n,y1n,y2n,y3n);
reta = -4:.1:4;
figure;
yk = plot(y(1:499),y(2:500),'.'); 
hold on;
yl = plot(sim1(1:499),sim1(2:500),'ro'); 
hold on;
plot(reta,reta,'b');
set(gca,'FontSize',11,'LineWidth',1);
set(gcf,'Name','Modelo com localização dos pontos fixos simétricos matematicamente garantidos','NumberTitle','off');
xlabel('y(k)');
ylabel('y(k+1)');
set(yk,'markersize',5);
set(yl,'markersize',5);
axis([-4 4 -4 4]);

errodm = [y(1:500) - sim1(1:500)].^2;
ssedm = sum(errodm);

%Evolução do maior expoente de Lyapunov

%Apenas com centros simétricos em relação aos pontos fixos não-triviais
sum = 0;
e=log10(2);
for ii = 2:10000
   J = X0(5) - (2*X0(1)/dp^2)*(sim(ii)-cn(1))*exp(-(sim(ii)-cn(1))^2/dp^2)-(2*X0(2)/dp^2)*(sim(ii)-cn(2))*exp(-(sim(ii)-cn(2))^2/dp^2) - (2*X0(3)/dp^2)*(sim(ii)-cn(3))*exp(-(sim(ii)-cn(3))^2/dp^2) - (2*X0(4)/dp^2)*(sim(ii)-cn(4))*exp(-(sim(ii)-cn(4))^2/dp^2);
   sum = sum + log10(abs(eig(J)))/e;
   %sum = sum + log(abs(eig(J)));
 	lia1(ii) = sum/ii;
end;

%Com imposição matemática da localização dos pontos fixos não-triviais
sum = 0;
e=log10(2);
for ii = 2:10000
   J = X1(5) - (2*X1(1)/dp^2)*(sim1(ii)-cn(1))*exp(-(sim1(ii)-cn(1))^2/dp^2)-(2*X1(2)/dp^2)*(sim1(ii)-cn(2))*exp(-(sim1(ii)-cn(2))^2/dp^2) - (2*X1(3)/dp^2)*(sim1(ii)-cn(3))*exp(-(sim1(ii)-cn(3))^2/dp^2) - (2*X1(4)/dp^2)*(sim1(ii)-cn(4))*exp(-(sim1(ii)-cn(4))^2/dp^2);
   sum = sum + log10(abs(eig(J)))/e;
   %sum = sum + log(abs(eig(J)));
 	lia2(ii) = sum/ii;
end;

%save parte2_05 lia cn dp X0 X1 pfs pfm ssepfs ssepfm ssepfs_pf ssepfm_pf ex_p sim sseds lia1 sim1 ssedm lia2;

%sim_p2_05 = sim;
%sim1_p2_05 = sim1;

%save parte2_05_fig y reta sim_p2_05 sim1_p2_05;
