clear all;
clc;

load init_resonant.dat
z=init_resonant(:,1);
x=init_resonant(:,3);
y=init_resonant(:,4);
plot(x,z,'--k',y,z,'-b')
ylim([-0.5 0.5])
print -dpdfcrop init_resonant.pdf