clear all
clc

load evolution.dat
x=evolution;
semilogy(x(:,1),x(:,2),'-k',x(:,1),x(:,3),'-b',x(:,1),x(:,4),'-r')
print -dpdfcrop energy.pdf
plot(x(:,1),x(:,5),'-k',x(:,1),x(:,6),'-b',x(:,1),x(:,7),'-r')
print -dpdfcrop rate.pdf

load ini_phys.dat
x=ini_phys;
plot(x(:,2),x(:,1),'-k',x(:,3),x(:,1),'--k',x(:,4),x(:,1),'-b',x(:,5),x(:,1),'--b',x(:,6),x(:,1),'-r',x(:,7),x(:,1),'--r')
ylim([-0.5 0.5])
print -dpdfcrop ini_phys.pdf

load ini_tor.dat
load ini_pol.dat
load ini_tem.dat
x=ini_tor;
y=ini_pol;
z=ini_tem;
semilogy(x(:,1),x(:,4),'-k',y(:,1),y(:,4),'-b',z(:,1),z(:,4),'-r')
print -dpdfcrop ini_spec.pdf

load fin_phys.dat
x=fin_phys;
plot(x(:,2),x(:,1),'-k',x(:,3),x(:,1),'--k',x(:,4),x(:,1),'-b',x(:,5),x(:,1),'--b',x(:,6),x(:,1),'-r',x(:,7),x(:,1),'--r')
ylim([-0.5 0.5])
print -dpdfcrop fin_phys.pdf

load fin_tor.dat
load fin_pol.dat
load fin_tem.dat
x=fin_tor;
y=fin_pol;
z=fin_tem;
semilogy(x(:,1),x(:,4),'-k',y(:,1),y(:,4),'-b',z(:,1),z(:,4),'-r')
print -dpdfcrop fin_spec.pdf