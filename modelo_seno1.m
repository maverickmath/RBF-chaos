function y = modelo_seno1(x,teta,c,dp);

y = x - teta(5)*x - teta(1)*exp(-(x-c(1))^2/dp^2) - teta(2)*exp(-(x-c(2))^2/dp^2) - ...
teta(3)*exp(-(x-c(3))^2/dp^2) - teta(4)*exp(-(x-c(4))^2/dp^2);