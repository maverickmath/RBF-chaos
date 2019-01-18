function y = modelo_buck(x,teta,c,dp);

y = x - teta(1) - teta(4)*x - teta(2)*exp(-(x-c(1))^2/dp^2)- teta(3)*exp(-(x-c(2))^2/dp^2);