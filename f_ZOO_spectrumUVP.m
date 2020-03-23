% function to compute spectrum from ZOOscan data
% input:
% - smin = minimum size
% - smax= maximum size
% - y = array of size measure
% - k = exponant of the logarithmic binning
% - opt = 1 for abundance spectrum
% - opt = 2 for volume spectra -> in this case array y must be in ESD

% made by Pieter Vandromme on june 2008

function [X,Y,Z,X1]=f_ZOO_spectrumUVP(sss,vol,smin,smax,kk,opt,conver)

x(1)=smin;t=1;
while x(t,1)<smax
    t=t+1;
    x(t,1)=x(1)*kk^(t-1);
end
x1=diff(x);
x=log(x);
% x2=diff(x);
X=(x(1:end-1,1)+x(2:end,1))./2; % choix de la taille nominale de chaque classe

if opt==1
    y1(1:size(sss,1),1)=1; %abundance
else
    y1=vol; %biovolume
end
n=size(x);
for t=1:n-1
    aa=find(((sss>=x(t,1)).*(sss<x(t+1,1)))==1);
    Y(t,1)=sum(y1(aa,1).*conver(aa,1));
end
Z=Y;%spectre non normalisé
Y=Y./x1; % normalization
X1=x1;