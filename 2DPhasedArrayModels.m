%% 2D Micropatch Array Field Calculations

clear all
close all

N = 16;
M = 8;

f = 2.45*10^9;
c = 3*10^8;
lam = c/f;
dz = 1/2*lam;
dy = 1/2*lam;
d = dz;

k = 2*pi/lam;

er = 10.7; %Duroid 6010

L = c/(f*2*sqrt(er));
W = c/(2*f*sqrt((er+1)/2));

angle = 90/360*2*pi;
angle2 = -10/360*2*pi;

pshift = k*dz*sin(angle);

phi = linspace(-pi/2, pi/2);
theta = linspace(0,2*pi);
[theta, phi] = meshgrid(theta,phi);

%Micropatch Radiation Pattern
Etheta = sin(0.5*k*W*sin(theta).*sin(phi))./(0.5*k*W.*sin(theta).*sin(phi)).*cos(0.5*k*L*sin(theta).*sin(phi)).*cos(phi);
Ephi = -sin(0.5*k*W*sin(theta).*sin(phi))./(0.5*k*W.*sin(theta).*sin(phi)).*cos(0.5*k*L*sin(theta).*sin(phi)).*cos(theta).*sin(phi);
ET = sqrt(Etheta.^2+Ephi.^2);

%1D Array
AF1 = sin(N*pi*d/lam*sin(theta))./(N*sin(pi*d/lam*sin(theta)));

%1D Array with Phase Shift
sint = sin(theta)-sin(angle);
cost = cos(theta)-cos(angle);
sinp = sin(phi)-sin(angle2);
cosp = cos(phi)-cos(angle2);

AF1P = sin(N*pi*d/lam*sint)./(N*sin(pi*d/lam*sint));

%2D Array
cosax = sin(theta).*cos(phi);
cosay = sin(theta).*sin(phi);
AF = (sin(N*k*dz*0.5.*cosax).*sin(M*k*dy*0.5.*cosay))./(N*sin(k*dz*0.5*cosax)*M.*sin(k*dy*0.5*cosay));

cosaxp = sin(theta).*cos(phi)-sin(angle).*cos(angle2);
cosayp = sin(theta).*sin(phi)-sin(angle).*sin(angle2);

AFP = (sin(N*k*dz*0.5.*cosaxp).*sin(M*k*dy*0.5.*cosayp))./(N*sin(k*dz*0.5*cosaxp)*M.*sin(k*dy*0.5*cosayp));

ElinArray = AF1.*ET;
ElinArrayP = AF1P.*ET;
EArray = AF.*ET;
EArrayP = AFP.*ET;


%% Gain
psi = k*d*sin(phi)*cos(theta);
psip = k*d*(sin(phi)-sin(angle))*(cos(theta)-cos(angle2));

gtheta = sin(N*psip/2)./(N*sin(psip/2)).*sin(N*d*k/2)*cos(theta);

[x,y,z] = sph2cart(theta,phi,ET);
[x1,y1,z1] = sph2cart(theta,phi,ElinArray);
[x1P,y1P,z1P] = sph2cart(theta,phi,ElinArrayP);
[x2,y2,z2] = sph2cart(theta,phi,EArray);
[x2P,y2P,z2P] = sph2cart(theta,phi,EArrayP);

r = x^2+y^2+z^2;

%% Directivity

M1 = max(max(ET.^2));
M2 = max(max(ElinArray.^2));
M2P = max(max(ElinArrayP.^2));
M3 = max(max(EArray.^2));
M3P = max(max(EArrayP.^2));

A1 = mean(ET.^2);
A2 = mean(ElinArray.^2);
A2P = mean(ElinArrayP.^2);

EArrayM = EArray(1:length(EArray)-1,2:length(EArray));
EArrayMP = EArray(1:length(EArray),2:length(EArray));

A3 = mean(mean(EArrayM.^2));
A3P = mean(mean(EArrayMP.^2));


D1 = M1/mean(A1(2:length(A1)));
d1db = 10*log10(D1);
D2 = M2/mean(A2(2:length(A2)));
d2db = 10*log10(D2);
D2P = M2/mean(A2P(2:length(A2P)));
d2Pdb = 10*log10(D2P);
D3 = M3/A3;
d3db = 10*log10(D3);
D3P = M3P/A3P;
d3Pdb = 10*log10(D3P);

%% dB
theta1 = 0:pi/32:2*pi;
phi1 = -10;

Etheta1 = sin(0.5*k*W*sin(theta1).*sin(phi1))./(0.5*k*W.*sin(theta1).*sin(phi1)).*cos(0.5*k*L*sin(theta1).*sin(phi1)).*cos(phi1);
Ephi1 = -sin(0.5*k*W*sin(theta1).*sin(phi1))./(0.5*k*W.*sin(theta1).*sin(phi1)).*cos(0.5*k*L*sin(theta1).*sin(phi1)).*cos(theta1).*sin(phi1);
ET1 = sqrt(Etheta1.^2+Ephi1.^2);

cosaxp = sin(theta1).*cos(phi1)-sin(angle).*cos(angle2);
cosayp = sin(theta1).*sin(phi1)-sin(angle).*sin(angle2);

AFP1 = (sin(N*k*dz*0.5.*cosaxp).*sin(M*k*dy*0.5.*cosayp))./(N*sin(k*dz*0.5*cosaxp)*M.*sin(k*dy*0.5*cosayp));

EdB = ET1.*AFP1;

polar(theta1,10*log10(EdB))
ylabel('dB')
%% Plots

figure (1)
surf(x,y,z.*d1db)
title('Radiation Pattern for Single Patch Antenna')


%figure (2)
%surf(x1,y1,z1.*d2db)
%title('Radiation Pattern for 16x1 Patch Antenna Array')

figure (3)
surf(x1P,y1P,z1P.*d2db)
title('Radiation Pattern for 16x1 Patch Antenna Array Targeted to 70 Degrees')


figure (6)
surf(x2P,y2P,z2P.*d3db)
title('Radiation Pattern for 16x4 Patch Antenna Array Targeted to 90 Degrees and 5 degrees')


%% AF Plots
[xaf,yaf,zaf] = sph2cart(theta,phi,AF1P);
figure (7)
surf(x1P,y1P,z1P.*d2db)
title('Array Factor for 16x1 Array Targeted to 90 degrees')

[xaf2,yaf2,zaf2] = sph2cart(theta,phi,AFP);
figure (7)
surf(xaf2,yaf2,zaf2.*d2db)
title('Array Factor for 16x4 Array Targeted to 90 degrees and -10 degrees')
