%Primeira parte: garantir um único ponto fixo (0) -- 2 centros
%(Simetria dos outros 2 também é garantida)

close all;
clear;
clc;
warning off;

nu = 0;   %Definido
ny = 1;   %Definido
nul = 0;   %Definido
nyl = 1;   %Definido
np = 2;	%Número de centros [valor "padrão": 2]
dp = 1.8;	
lin = 1;  %Presença ou não do polinômio linear
dist = 2; %Distância do ponto fixo [valor "padrão": 3]

[y,lia] = sinemap(1.2*pi,10000); %Mapa senoidal com não-linearidade cúbica (pf's: 0; +- 2,4383)

pf = 0; %Ponto fixo a impor
c = [pf - dist; pf + dist]; %Centros simétricos em relação ao ponto fixo

%Adicionar ruído [futuramente]

[Psi] = montaP([],y(1:500)', nu, ny, nul, nyl, lin, c, 500, np, dp); 

Psi_n = Psi(:,2:end);

%Estimativa de Mínimos Quadrados
X0 = ((inv((Psi_n')*Psi_n))*(Psi_n'))*y(max([ny nyl])+1:500)';

%Garante ponto fixo nulo [e simetria dos outros 2]
A = [1 1 0];
B = [0];
X1 = mqermod(Psi_n, X0, A, B); %Novo vetor de parâmetros

%Determina os pontos fixos dos modelos
prox = [-2.3 0 2.3];
for jj = 1:3,
   pfs(jj) = fzero('modelo_seno',prox(jj),[],X0,c,dp);
   pfm(jj) = fzero('modelo_seno',prox(jj),[],X1,c,dp);
end;

%Calcula os erros na localização dos pontos fixos
pfo = [-2.4383 0 2.4383];
erropfs = [pfo - pfs].^2;
erropfm = [pfo - pfm].^2;
ssepfs = sum(erropfs);
ssepfm = sum(erropfm);

%Para garantir o ponto fixo nulo basta fazer com que a soma dos parâmetros das
%RBFs seja nula

%Simulação

const = 0;
u1n = 0;
u2n = 0;
y1n = 1;
y2n = 0;
y3n = 0;
sim=simulacao_errn(c,length(c),0,1,10000,[],[pi/3],X0,lin,dp,const,u1n,u2n,y1n,y2n,y3n);
reta = -4:.1:4;
figure;
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

errods = [y(1:500) - sim(1:500)].^2;
sseds = sum(errods);

sim1 = simulacao_errn(c,length(c),0,1,10000,[],[pi/3],X1,lin,dp,const,u1n,u2n,y1n,y2n,y3n);
reta = -4:.1:4;
figure;
yk = plot(y(1:499),y(2:500),'.'); 
hold on;
yl = plot(sim1(1:499),sim1(2:500),'ro'); 
hold on;
plot(reta,reta,'b');
set(gca,'FontSize',11,'LineWidth',1);
set(gcf,'Name','Modelo com centros matematicamente simétricos [e ponto fixo trivial]','NumberTitle','off');
xlabel('y(k)');
ylabel('y(k+1)');
set(yk,'markersize',5);
set(yl,'markersize',5);

errodm = [y(1:500) - sim1(1:500)].^2;
ssedm = sum(errodm);

%Evolução do maior expoente de Lyapunov

%Sem imposição de simetria nos parâmetros [conseqüentemente, sem garantia do ponto fixo trivial]
sum = 0;
e=log10(2);
for ii = 2:10000
   %J = X0(4) + X0(2)*((-2*(sim(ii-1)-c(1)))/dp^2)*exp(-(sim(ii-1)-c(1))^2/dp^2) + X0(3)*((-2*(sim(ii-1)-c(2)))/dp^2)*exp(-(sim(ii-1)-c(2))^2/dp^2); 
   J = X0(3) - (2*X0(1)/dp^2)*(sim(ii)-c(1))*exp(-(sim(ii)-c(1))^2/dp^2)-(2*X0(2)/dp^2)*(sim(ii)-c(2))*exp(-(sim(ii)-c(2))^2/dp^2);
   sum = sum + log10(abs(eig(J)))/e;
   %sum = sum + log(abs(eig(J)));
 	lia1(ii) = sum/ii;
end;

%Com imposição do ponto fixo trivial [e de simetria dos outros 2]
sum = 0;
e=log10(2);
for ii = 2:10000
   %J = X0(4) + X0(2)*((-2*(sim(ii-1)-c(1)))/dp^2)*exp(-(sim(ii-1)-c(1))^2/dp^2) + X0(3)*((-2*(sim(ii-1)-c(2)))/dp^2)*exp(-(sim(ii-1)-c(2))^2/dp^2); 
   J = X1(3) - (2*X1(1)/dp^2)*(sim1(ii)-c(1))*exp(-(sim1(ii)-c(1))^2/dp^2)-(2*X1(2)/dp^2)*(sim1(ii)-c(2))*exp(-(sim1(ii)-c(2))^2/dp^2);
   sum = sum + log10(abs(eig(J)))/e;
   %sum = sum + log(abs(eig(J)));
 	lia2(ii) = sum/ii;
end;

%Agora gerando um modelo com centros selecionados pelo ERR
nrc = 200;
dpm = 2.5; %[Um valor viável é 5]
X = igespacocp(nrc, nu, ny, 500, [], y(1:500)); 
[Psi] = montaP([], y(1:500)', nu, ny, nul, nyl, lin, X, 500, length(X), dpm); 
[A,err,Piv]=myhouse(Psi,np); 

Piv_n = sort(Piv); 
Psi_n = Psi(:,Piv_n);

cm = X(Piv_n - 1);

Psi_a = montaP([],y(1:500)',0,0,0,1,1,[],500,0,0);
Psi_nn = [Psi_n Psi_a(:,2)];

Teta=((inv((Psi_nn')*Psi_nn))*(Psi_nn'))*y(max([ny nyl])+1:500)';

%Determina o ponto fixo do modelo
prox = [-2.3 0 2.3];
for jj = 1:3,
   pferr(jj) = fzero('modelo_seno',prox(jj),[],Teta,cm,dpm);
end;

%Calcula o erro na localização dos pontos fixos
erropferr = [pfo - pferr].^2;
ssepferr = erropferr(1) + erropferr(2) + erropferr(3);

const = 0;
u1n = 0;
u2n = 0;
y1n = 1;
y2n = 0;
y3n = 0;
sim2 = simulacao_errn(cm,length(cm),0,1,10000,[],[pi/3],Teta,lin,dpm,const,u1n,u2n,y1n,y2n,y3n);
reta = -4:.1:4;
figure;
yk = plot(y(1:499),y(2:500),'.'); 
hold on;
yl = plot(sim2(1:499),sim2(2:500),'ro'); 
hold on;
plot(reta,reta,'b');
set(gca,'FontSize',11,'LineWidth',1);
set(gcf,'Name','Modelo com centros escolhidos pelo ERR','NumberTitle','off');
xlabel('y(k)');
ylabel('y(k+1)');
set(yk,'markersize',5);
set(yl,'markersize',5);

erroderr = [y(1:500) - sim2(1:500)].^2;
ssederr = 0;
for ij = 1:length(erroderr),
   ssederr = ssederr + erroderr(ij);
end;

sum = 0;
e=log10(2);
for ii = 2:10000
   %J = X0(4) + X0(2)*((-2*(sim(ii-1)-c(1)))/dp^2)*exp(-(sim(ii-1)-c(1))^2/dp^2) + X0(3)*((-2*(sim(ii-1)-c(2)))/dp^2)*exp(-(sim(ii-1)-c(2))^2/dp^2); 
   J = Teta(3) - (2*Teta(1)/dpm^2)*(sim2(ii)-cm(1))*exp(-(sim2(ii)-cm(1))^2/dpm^2)-(2*Teta(2)/dpm^2)*(sim2(ii)-cm(2))*exp(-(sim2(ii)-cm(2))^2/dpm^2);
   sum = sum + log10(abs(eig(J)))/e;
   %sum = sum + log(abs(eig(J)));
 	lia3(ii) = sum/ii;
end;


%Plota evolução dos máximos expoentes de Lyapunov

figure;
plot([1:10000],lia,'b-',[1:10000],lia1,'k-',[1:10000],lia2,'r-',[1:10000],lia3,'m:');
set(gca,'FontSize',11,'LineWidth',1);
set(gcf,'Name','Evolução do máximo expoente de Lyapunov','NumberTitle','off');
axis([0 10000 0.8 1.4]);

%Salvando todas as variáveis relevantes
%save parte1 lia c dp X0 X1 pfs pfm ssepfs ssepfm sim sseds lia1 sim1 ssedm lia2 cm dpm Teta pferr ssepferr sim2 ssederr lia3;

%sim_p1 = sim;
%sim1_p1 = sim1;
%sim2_p1 = sim2;

%save fig_parte1 y reta sim_p1 sim1_p1 sim2_p1;