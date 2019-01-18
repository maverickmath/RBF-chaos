%Obtém modelo RBF com ponto fixo igual ao do sistema original [mapa de primeiro retorno de um regulador chaveado
%Buck]; compara com modelo cuja única intervenção é tomar centros simétricos ao ponto fixo conhecido

close all;
clear;
clc;
warning off;

nu = 0;   %Definido
ny = 1;   %Definido
nul = 0;   %Definido
nyl = 1;   %Definido
%%nrc = 300;	%Número de centros candidatos
np = 2;	%Número de centros [valor "padrão": 2]
dp = 2.4;	%Desvio-padrão/Largura da função-base [valor "viável": 4 p/ dist = 5]
%[5 p/ dist = 10, dando menos pontos (menores amplitudes) - 6,6.3 garante mais pontos 
%mas inicia a ondulação]
%[2.5 p/ dist = 3, c/ pequena ondulação e amplitudes um pouco maiores - 2.4, quase perfeito]
lin = 1;  %Presença ou não do polinômio linear
dist = 3; %Distância do ponto fixo [valor "padrão": 3]

%load Bucky.txt;
[y,lia]=buckmap(24,0.245,10000); %kapa=0.245 [parâmetro de bifurcação]
pf = 25.0018; %Ponto fixo regulador buck
c = [pf - dist; pf + dist]; %Centros simétricos em relação ao ponto fixo

%Adicionar ruído

[Psi] = montaP([],y(1:500)', nu, ny, nul, nyl, lin, c, 500, np, dp); 
%[Mas o termo de offset não causa nenhum impedimento à aplicação do método]

%Estimativa de Mínimos Quadrados
X0 = ((inv((Psi')*Psi))*(Psi'))*y(max([ny nyl])+1:500)';

%Simulação
const = 1;
u1n = 0;
u2n = 0;
y1n = 1;
y2n = 0;
y3n = 0;
sim=simulacao_errn(c,length(c),0,1,10000,[],[24],X0,lin,dp,const,u1n,u2n,y1n,y2n,y3n);
reta = 23:.1:33;
yk = plot(y(1:499),y(2:500),'.'); 
hold on;
yl = plot(sim(1:499),sim(2:500),'r.'); 
hold on;
plot(reta,reta,'b-');
hold on;
set(gca,'FontSize',11,'LineWidth',1);
set(gcf,'Name','Modelo apenas com centros simétricos','NumberTitle','off');
xlabel('y(k)');
ylabel('y(k+1)');
set(yk,'markersize',5);
set(yl,'markersize',5);
axis([23 33 23 33]);

pfs = fzero('modelo_buck',25,[],X0,c,dp);
%Este seria o ponto fixo do modelo

rg=23.5:0.25:33;
for i=1:length(rg)
	sime(i) = X0(1) + X0(4)*rg(i) + X0(2)*exp(-(rg(i)-c(1))^2/dp^2) + X0(3)*exp(-(rg(i)-c(2))^2/dp^2);
end;

plot(rg,sime,'k-');

%Evolução do maior expoente de Lyapunov 
sum = 0;
e=log10(2);
for ii = 2:10000
   %J = X0(4) + X0(2)*((-2*(sim(ii-1)-c(1)))/dp^2)*exp(-(sim(ii-1)-c(1))^2/dp^2) + X0(3)*((-2*(sim(ii-1)-c(2)))/dp^2)*exp(-(sim(ii-1)-c(2))^2/dp^2); 
   J = X0(4) - (2*X0(2)/dp^2)*(sim(ii)-c(1))*exp(-(sim(ii)-c(1))^2/dp^2)-(2*X0(3)/dp^2)*(sim(ii)-c(2))*exp(-(sim(ii)-c(2))^2/dp^2);
   %sum = sum + log10(abs(eig(J)))/e;
   sum = sum + log(abs(eig(J)));
 	lia1(ii) = sum/ii;
end;


A = [1 0 0 0; 0 1 1 0; 0 0 0 1]; %Matriz de mapeamento das restrições 
%[Considerando a ordem w0, centros, termos lineares]
b = (pf*(1 - X0(4)) - X0(1))/exp(-(dist^2)/(dp^2));
B = [X0(1); b; X0(4)]; %Vetor com os coeficientes

x = mqermod(Psi, X0, A, B); %Novo vetor de parâmetros

const = 1;
u1n = 0;
u2n = 0;
y1n = 1;
y2n = 0;
y3n = 0;
sim1=simulacao_errn(c,length(c),0,1,10000,[],[24],x,lin,dp,const,u1n,u2n,y1n,y2n,y3n);
figure;
reta = 23:.1:33;
yk = plot(y(1:499),y(2:500),'.'); 
hold on;
yl = plot(sim1(1:499),sim1(2:500),'r.'); 
hold on;
plot(reta,reta,'b-');
hold on;
set(gca,'FontSize',11,'LineWidth',1);
set(gcf,'Name','Modelo com ponto fixo garantido','NumberTitle','off');
xlabel('y(k)');
ylabel('y(k+1)');
set(yk,'markersize',5);
set(yl,'markersize',5);
axis([23 33 23 33]);

pfm = x(1) + x(4)*pf + (x(2) + x(3))*exp(-(dist)^2/dp^2);

for i=1:length(rg)
	sime1(i) = x(1) + x(4)*rg(i) + x(2)*exp(-(rg(i)-c(1))^2/dp^2) + x(3)*exp(-(rg(i)-c(2))^2/dp^2);
end;

plot(rg,sime1,'k-');

%Evolução do maior expoente de Lyapunov 
sum = 0;
e=log10(2);
for ii = 2:10000
   %J = x(4) + x(2)*((-2*(sim1(ii-1)-c(1)))/dp^2)*exp(-((sim1(ii-1)-c(1))^2/dp^2));
   %J = J + x(3)*((-2*(sim1(ii-1)-c(2)))/dp^2)*exp(-((sim1(ii-1)-c(2))^2/dp^2)); 
   J = x(4) - (2*x(2)/dp^2)*(sim1(ii)-c(1))*exp(-(sim1(ii)-c(1))^2/dp^2)-(2*x(3)/dp^2)*(sim1(ii)-c(2))*exp(-(sim1(ii)-c(2))^2/dp^2);
   %J=(2*2.6204*sim1(ii-1)-99.875)/sim1(ii-1)-(2.6204*sim1(ii-1)^2-99.875*sim1(ii-1)+1417.1)/sim1(ii-1)^2;
	%J=J-46.429*exp(22-sim1(ii-1));
   %sum = sum + log10(abs(eig(J)))/e;
   sum = sum + log(abs(eig(J)));
   lia2(ii) = sum/ii;
end;

%Agora gerando um modelo com centros selecionados pelo ERR
nrc = 200;
dpm = 5; %[Um valor viável é ?]
X = igespacocp(nrc, nu, ny, 500, [], y(1:500)); 
[Psi] = montaP([], y(1:500)', nu, ny, nul, nyl, lin, X, 500, length(X), dpm); 
[A,err,Piv]=myhouse(Psi,3); 

Piv_n = sort(Piv); 
Psi_n = Psi(:,Piv_n);

cm = X(Piv_n(2:end) - 1);

Psi_a = montaP([],y(1:500)',0,0,0,1,1,[],500,0,0);
Psi_nn = [Psi_n Psi_a(:,2)];

Teta=((inv((Psi_nn')*Psi_nn))*(Psi_nn'))*y(max([ny nyl])+1:500)';

%Determina o ponto fixo do modelo
pferr = fzero('modelo_buck',25,[],Teta,cm,dpm);

const = 1;
u1n = 0;
u2n = 0;
y1n = 1;
y2n = 0;
y3n = 0;
sim2 = simulacao_errn(cm,length(cm),0,1,10000,[],[24],Teta,lin,dpm,const,u1n,u2n,y1n,y2n,y3n);
figure;
yk = plot(y(1:499),y(2:500),'.'); 
hold on;
yl = plot(sim2(1:499),sim2(2:500),'r.'); 
hold on;
plot(reta,reta,'b');
hold on;
set(gca,'FontSize',11,'LineWidth',1);
set(gcf,'Name','Modelo com centros escolhidos pelo ERR','NumberTitle','off');
xlabel('y(k)');
ylabel('y(k+1)');
set(yk,'markersize',5);
set(yl,'markersize',5);
axis([23 33 23 33]);

for i=1:length(rg)
	sime2(i) = Teta(1) + Teta(4)*rg(i) + Teta(2)*exp(-(rg(i)-cm(1))^2/dpm^2) + Teta(3)*exp(-(rg(i)-cm(2))^2/dpm^2);
end;

plot(rg,sime2,'k-');

sum = 0;
e=log10(2);
for ii = 2:10000
   %J = X0(4) + X0(2)*((-2*(sim(ii-1)-c(1)))/dp^2)*exp(-(sim(ii-1)-c(1))^2/dp^2) + X0(3)*((-2*(sim(ii-1)-c(2)))/dp^2)*exp(-(sim(ii-1)-c(2))^2/dp^2); 
   J = Teta(4) - (2*Teta(2)/dpm^2)*(sim2(ii)-cm(1))*exp(-(sim2(ii)-cm(1))^2/dpm^2)-(2*Teta(3)/dpm^2)*(sim2(ii)-cm(2))*exp(-(sim2(ii)-cm(2))^2/dpm^2);
   %sum = sum + log10(abs(eig(J)))/e;
   sum = sum + log(abs(eig(J)));
 	lia3(ii) = sum/ii;
end;

figure;
plot([1:10000],lia,'b-',[1:10000],lia1,'k-',[1:10000],lia2,'r-',[1:10000],lia3,'m:');
set(gca,'FontSize',11,'LineWidth',1);
set(gcf,'Name','Evolução do máximo expoente de Lyapunov','NumberTitle','off');

%Salvando todas as variáveis relevantes
%save parte1_buck lia c dp X0 x pfs pfm sim sime lia1 sim1 sime1 lia2 cm dpm Teta pferr sim2 sime2 lia3;

%sim_p1 = sim;
%sime_p1 = sime;
%sim1_p1 = sim1;
%sime1_p1 = sime1;
%sim2_p1 = sim2;
%sime2_p1 = sime2;

%save fig_parte1_buck y reta sim_p1 sime_p1 sim1_p1 sime1_p1 sim2_p1 sime2_p1;


