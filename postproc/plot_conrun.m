% plot the probe result from conrun

clear all; close all;
[t,eta] = textread('../probe.dat','%f %f\n');
plot(t,eta,'r'); hold on;

[tp1,zp1] = textread('exp_data/briggs_front_probe.dat');
plot(tp1,zp1,'k+');
