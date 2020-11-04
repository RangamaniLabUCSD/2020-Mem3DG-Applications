clc;
clear;
close all;

bendenergy  = ncread('traj.nc','bendenergy');
surfenergy  = ncread('traj.nc','surfenergy');
pressenergy  = ncread('traj.nc','pressenergy');
kineenergy  = ncread('traj.nc','kineenergy');
chemenergy  = ncread('traj.nc','chemenergy');
totalenergy  = ncread('traj.nc','totalenergy');

time = ncread("traj.nc",'time');

figure
subplot(2,1,1);
hold on
plot(time,[bendenergy,surfenergy,pressenergy,chemenergy]);
legend('bending','surf','press','chem');
xlabel("time (s)")
ylabel("energy (10^{-15} J)")

title("Energy landscape with no local area localization");

subplot(2,1,2);
hold on
plot(time,[bendenergy + surfenergy + pressenergy + chemenergy, kineenergy, totalenergy]);
legend('potential','kinetic', 'total');
xlabel("time (s)")
ylabel("energy (10^{-15} J)")
