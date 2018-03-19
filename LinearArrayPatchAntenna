%% Linear Array Calculations

close all
clear all

%Models the radiation pattern for an N by 1 matrix of micropatch antenna
%elements

%Variables

V = 1; %Volts
Rn = 1;
N = 6;
f = 2.4*10^9; %Frequency (Hz)
c = 3*10^8; %Speed of Light

lam = c/f;
k = 2*pi/lam;
d = 1/4*lam;

angle =90/360*2*pi;
pshift = k*d*sin(angle);

%% Patch Antenna 3D

er = 10.7; %Duroid 6010

L = c/(f*2*sqrt(er));
W = L;

phi = linspace(-pi/2, pi/2);
theta = linspace(0,2*pi);

[theta, phi] = meshgrid(theta,phi);

Etheta = sin(0.5*k*W*sin(theta).*sin(phi))./(0.5*k*W.*sin(theta).*sin(phi)).*cos(0.5*k*W*sin(theta).*sin(phi)).*cos(phi);
Ephi = -sin(0.5*k*W*sin(theta).*sin(phi))./(0.5*k*W.*sin(theta).*sin(phi)).*cos(0.5*k*W*sin(theta).*sin(phi)).*cos(theta).*sin(phi);
ET = sqrt(Etheta.^2+Ephi.^2);
AF = cos((k*d*cos(theta)+pshift)/2);
AF1 = sin(N*pi*d/lam*(sin(theta)-sin(angle)))./(N*sin(pi*d/lam*(sin(theta)-sin(angle))));

EArray = AF1.*ET;

[x1,y1,z1] = sph2cart(theta,phi,ET);
[x2,y2,z2] = sph2cart(theta,phi,AF1);
[x3,y3,z3] = sph2cart(theta,phi,EArray);


figure (1)
surf(x1,y1,z1)
title('Radiation Pattern for Micropatch Antenna')

figure (2)
surf(x2,y2,z2)
title('Array Factor for 6x1 Linear Antenna Array')

figure (3)
surf(x3,y3,z3)
title('Radiation Patern for Micropatch Antenna in 6x1 Linear Array')
