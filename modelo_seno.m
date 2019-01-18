function y = modelo_seno(x,teta,c,dp);

y = x - teta(3)*x - teta(1)*exp(-(x-c(1))^2/dp^2)- teta(2)*exp(-(x-c(2))^2/dp^2);