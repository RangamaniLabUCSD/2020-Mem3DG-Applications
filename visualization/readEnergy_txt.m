clc;
clear;
close all;

load energy.txt;
figure
hold on
plot(energy(:,1),energy(:,2:4))
legend("bending",'surface','pressure');
figure 
hold on 
plot(energy(:,1),sum(energy(:,2:4)')',energy(:,1),energy(:,5),energy(:,1),energy(:,7));
legend("potential",'kinetic','total');