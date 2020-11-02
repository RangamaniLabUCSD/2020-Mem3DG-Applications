clc;
clear;
close all;

bendenergy  = ncread('traj.nc','bendenergy');
surfenergy  = ncread('traj.nc','surfenergy');
pressenergy  = ncread('traj.nc','pressenergy');
kineenergy  = ncread('traj.nc','kineenergy');
chemenergy  = ncread('traj.nc','chemenergy');
totalenergy  = ncread('traj.nc','totalenergy');
bendenergy = bendenergy(1:150);
surfenergy = surfenergy(1:150);
pressenergy = pressenergy(1:150);
kineenergy = kineenergy(1:150);
chemenergy = chemenergy(1:150);
totalenergy = totalenergy(1:150);

time = ncread("traj.nc",'time');
time = time(1:150);

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
