%% ZOOSCAN : process spectra and abundances for all groups
% this extract the following measure on the database:
% all particles abundance and biovolume (in spherical equivalent)
%   biovolume are in mm3/m3
%   spectra y unit is m3
%   to plot the spectra log transform the YY values
% made by Pieter Vandromme on july 2008
% modified and commented by jean-baptiste romagnan on january 2010
% adapted by Marc Picheral January 2013
% adapted for ecotaxa use by Fabien Lombard February 2016
% adapted for other instruments FL 2017
% adapted to calculate size quartiles + mean/std for each taxonomic group
% (FL January 2019)
function [base_spectres SStot] = process_abundances_biovolumes_NBSSspectra_IFCB(subset,base_spectres,smin,smax,k,uu,zoo_groups,i,fracnb);



%global ID aaa pred
%base_spectres = [];
test=strcmp(fracnb,'tot');

SStot = [];
warning off
% ----------- name of the groups (identifications of objects) -----------
zoo_groups = zoo_groups;
%zoo_groups = idnames;
m=length(zoo_groups);
disp([num2str(m),' groups']);
for j=1:m
    f=find(zoo_groups{j,:}=='-');
    zoo_groups{j,1}(f)='_';
end

% ------------ WAIT BAR -------------------------
%h=waitbar(0,'Computing SPECTRA ...');%create and display a bar that show progess of the analysis

%% inutile?
% % ------- Boucle sur les échantillons de la base des spectres ----------------
% n = length(base_spectres);
% % ---------- Fractions ------------------------------
% NBR = 4;
% fraction_names = {'da' 'db' 'dc' 'tot'};
% compteur_frac=0;
%% ---------- dimensions des vecteurs classes de taille ---------
x(1)=smin;
t=1;
while x(t,1)<smax
    t=t+1;
    x(t,1)=x(1)*k^(t-1);
end
x1=diff(x);
x=log(x);
X=(x(1:end-1,1)+x(2:end,1))./2; % choix de la classe milieu
X1=x1; % taille des classes
nb_int = length(X);
%%

method_list = {'Plain_Area_BV_spectra' 'Ellipsoid_BV_spectra' 'Summed_plain_Area_BV_spectra'};


for meth = 1:3
    clear mj mn sss vol pred ar aarea esd R3 area_int  aareaexc fferet % clear variables that will be reused later
    %% exctracting data 'area' 'major' 'minor' 'perimferet' 'area_exc'
    mj = subset.major; % (mm)
    mn = subset.minor; % (mm)               eval(['mn(:,1) = base_zooscan(base_spectres(t).',char(fraction_names(nbr)),'.base_zooscan_index).ids.Minor;']);% (mm²)
    %aareaexc = subset.area_exc; % (mm2)               eval(['aareaexc(:,1) = base_zooscan(base_spectres(t).',char(fraction_names(nbr)),'.base_zooscan_index).ids.Area_exc;']);   % (mm²)
    aarea = subset.area; % (mm2)               eval(['aarea(:,1) = base_zooscan(base_spectres(t).',char(fraction_names(nbr)),'.base_zooscan_index).ids.Area;']);% (mm²)
    %fferet = subset.perimferet; % (mm)               eval(['fferet(:,1)
    summedbiovolume=subset.summedbiovolume; %(mm3)
    %% calculating volumes
    if meth == 3;
        % ---------- calculation of summed biovolume from IFCB -----------
       
        sss=summedbiovolume;    % calculation of biovolume (mm3/m3)
    elseif meth == 2;
         % ---------- calculation of area of best fitting ellipse -----------
        ar=pi.*(mj./2).*(mn./2);
        sss=(4/3)*pi.*(mj./2).*(mn./2).*(mn./2);    % calculation of biovolume (mm3/m3)
    else
                area_int=aarea./pi;
        esd=2*(sqrt(area_int));
        R3=(esd./2).^3;
        sss=(4/3)*pi.*R3;
    end
    vol=double(sss);         % vol=biovolume
    sss=double(log(sss));    % log of biovolume
    pred=subset.object_annotation_hierarchy;
    
    for j=1:m % m = # de gpes de zoo
        id(:,j)=strcmp(zoo_groups{j,1},pred); % id matrix for initially identified groups
        id=double(id);
    end
    ID = id; % matrice des indices d'identification de chaque objets l = size(pred,1) et c = # de gpes: on a un seul 1 par ligne!!!
    
    conver = subset.conver; %per cubic meter
    if test==1
        eval(['base_spectres(i).tot.conver=conver;'])
    else
        eval(['base_spectres(i).d' num2str(fracnb) '.conver=conver;'])
    end
    for j=1:m
        % -------- conversion factor to have all in cubic meter ---
        % eval(['conver = base_spectres(t).',char(fraction_names(nbr)),'.Fracnb/(base_zooscan(base_spectres(t).',char(fraction_names(nbr)),'.base_zooscan_index).Vol * base_zooscan(base_spectres(t).',char(fraction_names(nbr)),'.base_zooscan_index).Townb);']);
        
        
        if test==1
            Ab(j)=sum(ID(:,j)).*conver; % abundance per fraction rapportée au volume
            eval(['base_spectres(i).tot.Ab(j,1)=Ab(j);'])              % abundance per fraction rapportée au volume (#/m3)
            Abtot(j)=sum(ID(:,j));    % abondance totale sur toute la profondeur echantillonnée (= # d'organismes EN TOUT, pas rapporté au volume filtré)
            eval(['base_spectres(i).tot.Abtot(j,1)=Abtot(j);'])    % # d'organismes EN TOUT pour la fraction scannée, pas rapporté au volume filtré
            
        else
            Ab(j)=sum(ID(:,j)).*conver; % abundance per fraction rapportée au volume
            eval(['base_spectres(i).d' num2str(fracnb) '.Ab(j,1)=Ab(j);'])              % abundance per fraction rapportée au volume (#/m3)
            Abtot(j)=sum(ID(:,j));    % abondance totale sur toute la profondeur echantillonnée (= # d'organismes EN TOUT, pas rapporté au volume filtré)
            eval(['base_spectres(i).d' num2str(fracnb) '.Abtot(j,1)=Abtot(j);'])    % # d'organismes EN TOUT pour la fraction scannée, pas rapporté au volume filtré
        end
        
        
        if meth == 1
            sizes(:,j)=ID(:,j).*esd; % abundance per fraction rapportée au volume
            
            %eval(['base_spectres(i).d' num2str(fracnb) '.sizes(:,j)=sizes(:,j);'])
        end
        
        if meth == 3
            if test==1
                Bv(j)=sum(ID(:,j).*vol).*conver ;% biovolume per fraction rapportée au volume
                eval(['base_spectres(i).tot.Bv(j,1)=Bv(j);'])               % abundance per fraction rapportée au volume (#/m3)
                Bvtot(j)=sum(ID(:,j).*vol);    % biovolume total sur toute la profondeur echantillonnée (= # d'organismes EN TOUT, pas rapporté au volume filtré)
                eval(['base_spectres(i).tot.Bvtot(j,1)=Bvtot(j);'])    % # d'organismes EN TOUT pour la fraction scannée, pas rapporté au volume filtré
            else
                Bv(j)=sum(ID(:,j).*vol).*conver ;% biovolume per fraction rapportée au volume
                eval(['base_spectres(i).d' num2str(fracnb) '.Bv(j,1)=Bv(j);'])               % abundance per fraction rapportée au volume (#/m3)
                Bvtot(j)=sum(ID(:,j).*vol);    % biovolume total sur toute la profondeur echantillonnée (= # d'organismes EN TOUT, pas rapporté au volume filtré)
                eval(['base_spectres(i).d' num2str(fracnb) '.Bvtot(j,1)=Bvtot(j);'])    % # d'organismes EN TOUT pour la fraction scannée, pas rapporté au volume filtré
            end
        end
        
        aaa=find(ID(:,j)==1);
        f=find(zoo_groups{j,1}=='-');
        zoo_groups{j,1}(f)='_';
        % ----------- spectra per fraction -------------
        % ----------- Presence d'organismes de cette categorie pour cette fraction -------------
        if isempty(aaa) == 0;
            
            % ---------- Abundances ---------------------
%             eval(['[X,Yab_zoo_groups_',num2str(j),',Z,X1]=f_ZOO_spectrum(sss(aaa,1),vol(aaa,1),smin,smax,k,1);']);% clear aaa
%             eval(['Yab_zoo_groups_',num2str(j),'=Yab_zoo_groups_',num2str(j),'.*conver;']) %'Ygroup'= size distribution for a group for a fraction
%             
            % ---------- Biovolume ---------------------
            eval(['[X,Ybv_zoo_groups_',num2str(j),',Z,X1]=f_ZOO_spectrum(sss(aaa,1),vol(aaa,1),smin,smax,k,2);']);% clear aaa
            eval(['Ybv_zoo_groups_',num2str(j),'=Ybv_zoo_groups_',num2str(j),'.*conver;']) %'Ygroup'= size distribution for a group for a fraction
            
            % ----------- Pas d'organismes de cette categorie pour cette fraction ---------
        elseif isempty(aaa) == 1;
            %eval(['Yab_zoo_groups_',num2str(j),'=zeros(size(X,1),1);']);
            eval(['Ybv_zoo_groups_',num2str(j),'=zeros(size(X,1),1);']);
        end
        if test==1
            % ----------- On charge les données pour les fractions dans la base ---------
            %eval(['base_spectres(i).tot.Yab(:,j) = Yab_zoo_groups_',num2str(j),';']);
            eval(['base_spectres(i).tot.Ybv_',char(method_list(meth)),'(:,j) = Ybv_zoo_groups_',num2str(j),';']);
            
        else
            % ----------- On charge les données pour les fractions dans la base ---------
            %eval(['base_spectres(i).d' num2str(fracnb) '.Yab(:,j) = Yab_zoo_groups_',num2str(j),';']);
            eval(['base_spectres(i).d' num2str(fracnb) '.Ybv_',char(method_list(meth)),'(:,j) = Ybv_zoo_groups_',num2str(j),';']);
        end
        

    end

        I=sizes==0;
    sizes(I)=NaN;
    quartiles=quantile(sizes,[0.05 0.25 0.5 0.75 0.95]);
    moyenne=nanmean(sizes);
    ecartype=nanstd(sizes);
    if size(quartiles,1)==1
        quartilesmoyenne=[quartiles';moyenne;ecartype];
    else
    quartilesmoyenne=[quartiles;moyenne;ecartype];
    end
    if test==1
        eval(['base_spectres(i).tot.ESDquartilesmoyenne=quartilesmoyenne;'])
    else
        eval(['base_spectres(i).d' num2str(fracnb) '.ESDquartilesmoyenne=quartilesmoyenne;']) 
    end
end

if test==1
    eval(['base_spectres(i).tot.X = X;']);
    eval(['base_spectres(i).tot.X1 = X1;']);
    eval(['base_spectres(i).tot.Zoo_groups = zoo_groups;']);
else
    eval(['base_spectres(i).d' num2str(fracnb) '.X = X;']);
    eval(['base_spectres(i).d' num2str(fracnb) '.X1 = X1;']);
    eval(['base_spectres(i).d' num2str(fracnb) '.Zoo_groups = zoo_groups;']);
end