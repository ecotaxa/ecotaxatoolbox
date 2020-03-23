function [flux,classes,DSE,tailles,varargout]=nb2fluxECOpart(X,mini,maxi,pas,A,B,varargin)
% Estimate carbon flux from particle size distribution (ECOpart)
%
%           [flux,classes,DSE,sizeclasses]=nb2fluxECOpart(X,mini,maxi,pas,A,B)
% Input:
%   X : Concentration of particles (Nb/l per size classes)
%   mini : Minimum size for particles (1 micron) .. in mm
%   maxi : Maximum size for particles (26 cm) .. in mm
%   pas : Geometric progression for the size classes (2^1/3 pour ECOpart)
%   A : Coefficient from Guidi et al., 2008 (default 12.5 for Carbon)
%   B : Exponent from Guidi et al., 2008    (defaut 3.81 for Carbon)
% 
% Output:
%   flux : Matrix of fluxes en mg.m^(-2).d^(-1) : Tot Flux (classes 25 to 32)
%   classes : Size classes in mm
%   DSE : Mean of the size classes in mm
%   sizeclasses : Size of the size bins in mm


[a,b]=size(X);
flux=zeros(a,b);

[nl,nc]=size(X);
if nargin==1            % for the 45 size classes (ECOpart)
    mini=1*10^-3;       % in mm
    maxi=30;            % in mm
    pas=2^(1/3);
end

if nargin<5
    A=12.5;
    B=3.81;
end

% Direct estimates of classes in mm
[classe,tailles,DSE]=pasvar(mini,maxi,pas);
classes=[classe(1:end-1)' classe(2:end)'];

% Calcul du flux
flu=A*DSE.^B;        % Formule pour le calcul du flu basee sur les minimisation sur les pieges a sediment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fin de la Troisieme partie%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:nl
    flux(i,:)=X(i,:).*flu;                     % Matrice des flux (en mgC.m^(-2).j^(-1))
end