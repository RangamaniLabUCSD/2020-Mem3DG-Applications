clc;
clear;
close all;

max = 50;
bendenergy  = ncread('traj.nc','bendenergy');
surfenergy  = ncread('traj.nc','surfenergy');
pressenergy  = ncread('traj.nc','pressenergy');
kineenergy  = ncread('traj.nc','kineenergy');
chemenergy  = ncread('traj.nc','chemenergy');
totalenergy  = ncread('traj.nc','totalenergy');
bendenergy = bendenergy(1:max);
surfenergy = surfenergy(1:max);
pressenergy = pressenergy(1:max);
kineenergy = kineenergy(1:max);
chemenergy = chemenergy(1:max);
totalenergy = totalenergy(1:max);

time = ncread("traj.nc",'time');
time = time(1:max);

figure
subplot(2,1,1);
hold on
plot(time,[bendenergy,surfenergy,pressenergy,chemenergy]);
legend('bending','surf','press','chem');
xlabel("time (s)")
ylabel("energy (10^{-15} J)")

title("Energy landscape with conformal regularization");

subplot(2,1,2);
hold on
plot(time,[bendenergy + surfenergy + pressenergy + chemenergy, kineenergy, totalenergy]);
legend('potential','kinetic', 'total');
xlabel("time (s)")
ylabel("energy (10^{-15} J)")

saveas(gcf,"conformal_regularization.png")
