clear all;
clc;

load 14/fort.2
plot(fort(:,2),fort(:,1),'-k',fort(:,3),fort(:,1),'--k',fort(:,4),fort(:,1),'-b',fort(:,5),fort(:,1),'--b',fort(:,6),fort(:,1),'-r',fort(:,7),fort(:,1),'--r')
ylim([-0.5 0.5])
print -dpdfcrop 14/final.pdf