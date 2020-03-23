function [ves,new,med]=pasvar(minor,maxor,pas)
%           [ves,new,med]=pasvar(minor,maxor,pas)
%
% Avec :
% Construction d'un vecteur à pas variable pour les histogrammes (ves)
% Construction d'un vecteur mesurant la taille des intervals entre les valeurs extremes (new)
% Construction d'un vecteur avec les valeurs médianes des pas (med)


ves=minor;
n=1000;
for i=2:n
    ves(i)=pas*ves(i-1);
    if ves(i)>maxor
        break
    end
end

ves; 
h=length(ves);
sous=[ves(:,2:end)];
new=sous-ves(1:h-1);

for i=1:h-1
    med(i)=(ves(i)+ves(i+1))/2;
end