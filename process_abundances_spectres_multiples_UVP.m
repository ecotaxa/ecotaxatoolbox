%% ZOOSCAN : process spectra and abundances for all groups
% this extract the following measure on the database:
% all particles abundance and biovolume (in spherical equivalent)
%   biovolume are in mm3/m3
%   spectra y unit is m3
%   to plot the spectra log transform the YY values
% made by Pieter Vandromme on july 2008
% modified and commented by jean-baptiste romagnan on january 2010
% adapted by Marc Picheral January 2013
% adapted by Fabien Lombard 2016-2018

function [base_spectres SStot] = process_abundances_spectres_multiples_UVP(subset,base_spectres,smin,smax,k,uu,zoo_groups,i,fracnb);
%base_spectres=base_Zooscan;

%global ID aaa pred
%base_spectres = [];
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
% % ------- Boucle sur les �chantillons de la base des spectres ----------------
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

method_list = {'Plain_Area_BV_spectra' 'Riddled_Area_BV_spectra' 'Ellipsoid_BV_spectra'};


for meth = 1:3
    clear mj mn sss vol pred ar aarea esd R3 area_int  aareaexc fferet % clear variables that will be reused later
    %% exctracting data 'area' 'major' 'minor' 'perimferet' 'area_exc'
    mj = subset.major; % (mm)
    mn = subset.minor; % (mm)               eval(['mn(:,1) = base_zooscan(base_spectres(t).',char(fraction_names(nbr)),'.base_zooscan_index).ids.Minor;']);% (mm�)
    aareaexc = subset.area_exc; % (mm2)               eval(['aareaexc(:,1) = base_zooscan(base_spectres(t).',char(fraction_names(nbr)),'.base_zooscan_index).ids.Area_exc;']);   % (mm�)
    aarea = subset.area; % (mm2)               eval(['aarea(:,1) = base_zooscan(base_spectres(t).',char(fraction_names(nbr)),'.base_zooscan_index).ids.Area;']);% (mm�)
    fferet = subset.perimferet; % (mm)               eval(['fferet(:,1)
    %% calculating volumes
    if meth == 3;
        % ---------- calculation of area of best fitting ellipse -----------
        ar=pi.*(mj./2).*(mn./2);
        sss=(4/3)*pi.*(mj./2).*(mn./2).*(mn./2);    % calculation of biovolume (mm3/m3)
    elseif meth == 2;
        area_int=aareaexc./pi;
        esd=2*(sqrt(area_int));
        R3=(esd./2).^3;
        sss=(4/3)*pi.*R3;
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
    eval(['base_spectres(i).tot.conver=conver;'])
    for j=1:m
        % -------- conversion factor to have all in cubic meter ---
        % eval(['conver = base_spectres(t).',char(fraction_names(nbr)),'.Fracnb/(base_zooscan(base_spectres(t).',char(fraction_names(nbr)),'.base_zooscan_index).Vol * base_zooscan(base_spectres(t).',char(fraction_names(nbr)),'.base_zooscan_index).Townb);']);
        
        
        Ab(j)=sum(ID(:,j).*conver); % abundance per fraction rapport�e au volume
        eval(['base_spectres(i).tot.Ab(j,1)=Ab(j);'])              % abundance per fraction rapport�e au volume (#/m3)
        Abtot(j)=sum(ID(:,j));    % abondance totale sur toute la profondeur echantillonn�e (= # d'organismes EN TOUT, pas rapport� au volume filtr�)
        eval(['base_spectres(i).tot.Abtot(j,1)=Abtot(j);'])    % # d'organismes EN TOUT pour la fraction scann�e, pas rapport� au volume filtr�
        
        if meth == 3
            Bv(j)=sum(ID(:,j).*vol.*conver) ;% biovolume per fraction rapport�e au volume
            eval(['base_spectres(i).tot.Bv(j,1)=Bv(j);'])               % abundance per fraction rapport�e au volume (#/m3)
            Bvtot(j)=sum(ID(:,j).*vol);    % biovolume total sur toute la profondeur echantillonn�e (= # d'organismes EN TOUT, pas rapport� au volume filtr�)
            eval(['base_spectres(i).tot.Bvtot(j,1)=Bvtot(j);'])    % # d'organismes EN TOUT pour la fraction scann�e, pas rapport� au volume filtr�
        end
        
        aaa=find(ID(:,j)==1);
        f=find(zoo_groups{j,1}=='-');
        zoo_groups{j,1}(f)='_';
        % ----------- spectra per fraction -------------
        % ----------- Presence d'organismes de cette categorie pour cette fraction -------------
        if isempty(aaa) == 0;
            
            % ---------- Abundances ---------------------
            eval(['[X,Yab_zoo_groups_',num2str(j),',Z,X1]=f_ZOO_spectrumUVP(sss(aaa,1),vol(aaa,1),smin,smax,k,1,conver(aaa));']);% clear aaa
            eval(['Yab_zoo_groups_',num2str(j),'=Yab_zoo_groups_',num2str(j),';']) %'Ygroup'= size distribution for a group for a fraction
            
            % ---------- Biovolume ---------------------
            eval(['[X,Ybv_zoo_groups_',num2str(j),',Z,X1]=f_ZOO_spectrumUVP(sss(aaa,1),vol(aaa,1),smin,smax,k,2,conver(aaa));']);% clear aaa
            eval(['Ybv_zoo_groups_',num2str(j),'=Ybv_zoo_groups_',num2str(j),';']) %'Ygroup'= size distribution for a group for a fraction
            
            % ----------- Pas d'organismes de cette categorie pour cette fraction ---------
        elseif isempty(aaa) == 1;
            eval(['Yab_zoo_groups_',num2str(j),'=zeros(size(X,1),1);']);
            eval(['Ybv_zoo_groups_',num2str(j),'=zeros(size(X,1),1);']);
        end
        % ----------- On charge les donn�es pour les fractions dans la base ---------
        eval(['base_spectres(i).tot.Yab(:,j) = Yab_zoo_groups_',num2str(j),';']);
        eval(['base_spectres(i).tot.Ybv_',char(method_list(meth)),'(:,j) = Ybv_zoo_groups_',num2str(j),';']);
        
    end
end


eval(['base_spectres(i).tot.X = X;']);
eval(['base_spectres(i).tot.X1 = X1;']);
eval(['base_spectres(i).tot.Zoo_groups = zoo_groups;']);
    
    
return
for t=1:n; % boucle sur les �chantillons 
    
    waitbar(t / n)
    disp([num2str(t),' / ',num2str(n),'  :  ',char(base_spectres(t).SampleId)]);
    %strcmp(base_zooscan(2).ids.tag_nb_objs,'ok')
    
    % --------------- Boucle sur les 3 methodes de calcul du BV -------
    for meth = 1:3
        % ---------- Boucle sur les fractions pour 1 �chantillon ------
        for nbr=1:NBR
            clear mj mn sss vol pred ar aarea esd R3 area_int  aareaexc fferet % clear variables that will be reused later
            % ---------------------- Test existence de donn�es pour la fraction ---------------
            eval(['testpred = isempty(base_spectres(t).',char(fraction_names(nbr)),');']);
            if testpred == 0;
                compteur_frac=compteur_frac+1;
                % ----------- Fraction contains data -------------------------
                eval(['mj(:,1) = base_zooscan(base_spectres(t).',char(fraction_names(nbr)),'.base_zooscan_index).ids.Major;']); % (mm)
                eval(['mn(:,1) = base_zooscan(base_spectres(t).',char(fraction_names(nbr)),'.base_zooscan_index).ids.Minor;']);% (mm�)
                eval(['aareaexc(:,1) = base_zooscan(base_spectres(t).',char(fraction_names(nbr)),'.base_zooscan_index).ids.Area_exc;']);   % (mm�)
                eval(['aarea(:,1) = base_zooscan(base_spectres(t).',char(fraction_names(nbr)),'.base_zooscan_index).ids.Area;']);% (mm�)
                eval(['fferet(:,1) = base_zooscan(base_spectres(t).',char(fraction_names(nbr)),'.base_zooscan_index).ids.Feret;']);% (mm�)
                if meth == 3;
                    % ---------- calculation of area of best fitting ellipse -----------
                    ar=pi.*(mj./2).*(mn./2);
                    sss=(4/3)*pi.*(mj./2).*(mn./2).*(mn./2);    % calculation of biovolume (mm3/m3)
                elseif meth == 2;
                    area_int=aareaexc./pi;
                    esd=2*(sqrt(area_int));
                    R3=(esd./2).^3;
                    sss=(4/3)*pi.*R3;
                else
                    area_int=aarea./pi;
                    esd=2*(sqrt(area_int));
                    R3=(esd./2).^3;
                    sss=(4/3)*pi.*R3;
                end
                vol=double(sss);         % vol=biovolume
                sss=double(log(sss));    % log of biovolume
                % ----------------- definition of 'pred' variable corresponding to the identifications
                eval(['pred = base_zooscan(base_spectres(t).',char(fraction_names(nbr)),'.base_zooscan_index).ids.Pred;']);
                eval(['base_spectres(t).',char(fraction_names(nbr)),'.Raw_Ids=pred;']);
                eval(['base_spectres(t).',char(fraction_names(nbr)),'.Raw_Major=double(mj);']);
                eval(['base_spectres(t).',char(fraction_names(nbr)),'.Raw_Minor=double(mn);']);
                eval(['base_spectres(t).',char(fraction_names(nbr)),'.Raw_Area_exc=double(aareaexc);']);
                eval(['base_spectres(t).',char(fraction_names(nbr)),'.Raw_Area=double(aarea);']);
                eval(['base_spectres(t).',char(fraction_names(nbr)),'.Raw_Feret=double(fferet);']);
                
                % -------------- Boucle sur les groupes --------------------
                
                for j=1:m % m = # de gpes de zoo
                    id(:,j)=strcmp(zoo_groups{j,1},pred); % id matrix for initially identified groups
                    id=double(id);
                end
                
                clear j
                ID = id; % matrice des indices d'identification de chaque objets l = size(pred,1) et c = # de gpes: on a un seul 1 par ligne!!!
                if isequal(size(vol,1),size(ID,1))==1
                    % -------------- Boucle sur les groupes --------------------
                    for j=1:m
                        % -------- conversion factor to have all in cubic meter ---
                        % eval(['conver = base_spectres(t).',char(fraction_names(nbr)),'.Fracnb/(base_zooscan(base_spectres(t).',char(fraction_names(nbr)),'.base_zooscan_index).Vol * base_zooscan(base_spectres(t).',char(fraction_names(nbr)),'.base_zooscan_index).Townb);']);
                        eval(['conver = base_spectres(t).',char(fraction_names(nbr)),'.Fracnb/(base_zooscan(base_spectres(t).',char(fraction_names(nbr)),'.base_zooscan_index).Vol);']);
                        
                        eval(['base_spectres(t).',char(fraction_names(nbr)),'.conver=conver;']);
                        Ab(nbr,j)=sum(ID(:,j)).*conver; % abundance per fraction rapport�e au volume
                        eval(['base_spectres(t).',char(fraction_names(nbr)),'.Ab(j,1)=Ab(nbr,j);']);               % abundance per fraction rapport�e au volume (#/m3)
                        Abtot(nbr,j)=sum(ID(:,j));    % abondance totale sur toute la profondeur echantillonn�e (= # d'organismes EN TOUT, pas rapport� au volume filtr�)
                        eval(['base_spectres(t).',char(fraction_names(nbr)),'.Abtot(j,1)=Abtot(nbr,j);']);    % # d'organismes EN TOUT pour la fraction scann�e, pas rapport� au volume filtr�
                        
                        if meth == 3
                            Bv(nbr,j)=sum(ID(:,j).*vol).*conver ;% biovolume per fraction rapport�e au volume
                            eval(['base_spectres(t).',char(fraction_names(nbr)),'.Bv(j,1)=Bv(nbr,j);']);               % abundance per fraction rapport�e au volume (#/m3)
                            Bvtot(nbr,j)=sum(ID(:,j).*vol);    % biovolume total sur toute la profondeur echantillonn�e (= # d'organismes EN TOUT, pas rapport� au volume filtr�)
                            eval(['base_spectres(t).',char(fraction_names(nbr)),'.Bvtot(j,1)=Bvtot(nbr,j);']);    % # d'organismes EN TOUT pour la fraction scann�e, pas rapport� au volume filtr�
                        end
                        
                        aaa=find(ID(:,j)==1);
                        f=find(zoo_groups{j,1}=='-');
                        zoo_groups{j,1}(f)='_';
                        % ----------- spectra per fraction -------------
                        % ----------- Presence d'organismes de cette categorie pour cette fraction -------------
                        if isempty(aaa) == 0;
                            
                            % ---------- Abundances ---------------------
                            eval(['[X,Yab_zoo_groups_',num2str(j),'(:,nbr),Z,X1]=f_ZOO_spectrum(sss(aaa,1),vol(aaa,1),smin,smax,k,1);']);% clear aaa
                            eval(['Yab_zoo_groups_',num2str(j),'(:,nbr)=Yab_zoo_groups_',num2str(j),'(:,nbr).*conver;']) %'Ygroup'= size distribution for a group for a fraction

                            % ---------- Biovolume ---------------------
                            eval(['[X,Ybv_zoo_groups_',num2str(j),'(:,nbr),Z,X1]=f_ZOO_spectrum(sss(aaa,1),vol(aaa,1),smin,smax,k,2);']);% clear aaa
                            eval(['Ybv_zoo_groups_',num2str(j),'(:,nbr)=Ybv_zoo_groups_',num2str(j),'(:,nbr).*conver;']) %'Ygroup'= size distribution for a group for a fraction
                            
                            % ----------- Pas d'organismes de cette categorie pour cette fraction ---------
                        elseif isempty(aaa) == 1;
                            eval(['Yab_zoo_groups_',num2str(j),'(:,nbr)=zeros(size(X,1),1);']);
                            eval(['Ybv_zoo_groups_',num2str(j),'(:,nbr)=zeros(size(X,1),1);']);
                        end
                        % ----------- On charge les donn�es pour les fractions dans la base ---------
                        eval(['base_spectres(t).',char(fraction_names(nbr)),'.Yab(:,j) = Yab_zoo_groups_',num2str(j),'(:,nbr);']);
                        eval(['base_spectres(t).',char(fraction_names(nbr)),'.Ybv_',char(method_list(meth)),'(:,j) = Ybv_zoo_groups_',num2str(j),'(:,nbr);']);
                        
                    end
                    
                    clear id aaa pred f% clean variables that will be reused later
   
                else
                    for j=1:m
                        eval(['conver = base_spectres(t).',char(fraction_names(nbr)),'.Fracnb/(base_zooscan(base_spectres(t).',char(fraction_names(nbr)),'.base_zooscan_index).Vol);']);
                        eval(['base_spectres(t).',char(fraction_names(nbr)),'.conver=conver;']);
                        
                        % abondances
                        Ab(nbr,j)=nan;
                        Abtot(nbr,j)=nan;
                        
                        eval(['base_spectres(t).',char(fraction_names(nbr)),'.Ab(j,1)=Ab(nbr,j);']);               % abundance per fraction rapport�e au volume (#/m3)
                        eval(['base_spectres(t).',char(fraction_names(nbr)),'.Abtot(j,1)=Abtot(nbr,j);']);    % # d'or
                        
                        % biovolumes
                        Bv(nbr,j)=nan;
                        Bvtot(nbr,j)=nan;
                        
                        eval(['base_spectres(t).',char(fraction_names(nbr)),'.Bv(j,1)=Bv(nbr,j);']);               % abundance per fraction rapport�e au volume (#/m3)
                        eval(['base_spectres(t).',char(fraction_names(nbr)),'.Bvtot(j,1)=Bvtot(nbr,j);']);    % # d'organismes EN TOUT pour la fraction scann�e, pas rapport� au volume filtr�
                        
                        % spectres
                        eval(['Yab_zoo_groups_',num2str(j),'(:,nbr)=nan(size(X,1),1);']);
                        eval(['Ybv_zoo_groups_',num2str(j),'(:,nbr)=nan(size(X,1),1);']);
                        
                        eval(['base_spectres(t).',char(fraction_names(nbr)),'.Yab(:,j) = Yab_zoo_groups_',num2str(j),'(:,nbr);']);
                        eval(['base_spectres(t).',char(fraction_names(nbr)),'.Ybv_',char(method_list(meth)),'(:,j) = Ybv_zoo_groups_',num2str(j),'(:,nbr);']);
                    end
                end
            end
        end
        
        Ab1=nan_sum(Ab,1)';
        Abtot1=nan_sum(Abtot,1)';
        if meth == 3
            Bv1=nan_sum(Bv,1)';
            Bvtot1=nan_sum(Bvtot,1);
        end
        % ----------- On travaille maintenant a rassembler les fractions ---------
        for j=1:m
            
            eval(['Yab_zoo_groups_',num2str(j),'= nan_sum(Yab_zoo_groups_',num2str(j),',2);']) % sum of the several fractions spectra for a sample
            eval(['Yab_' ,char(method_list(meth)),'(:,j)=Yab_zoo_groups_',num2str(j),';']);
            eval(['Ybv_zoo_groups_',num2str(j),'= nan_sum(Ybv_zoo_groups_',num2str(j),',2);']) % sum of the several fractions spectra for a sample
            eval(['Ybv_' ,char(method_list(meth)),'(:,j)=Ybv_zoo_groups_',num2str(j),';']);
        end
        
        %clear id
        clear Ab Abtot Bv Bvtot
    end
    
    base_spectres(t).combined_frac.smin = smin;
    base_spectres(t).combined_frac.smax = smax;
    base_spectres(t).combined_frac.k = k;
    base_spectres(t).combined_frac.uu = uu;
    base_spectres(t).combined_frac.X = X;
    base_spectres(t).combined_frac.X1 = X1;
    base_spectres(t).combined_frac.Zoo_groups = zoo_groups;
    base_spectres(t).combined_frac.Ab = Ab1;   % abundance of groups for a station, all fractions combined
    base_spectres(t).combined_frac.Abtot = Abtot1;%
    base_spectres(t).combined_frac.Bv = Bv1;   % biovolume of groups for a station, all fractions combined
    base_spectres(t).combined_frac.Bvtot = Bvtot1;%
    eval(['base_spectres(t).combined_frac.Abundance_spectra=Yab_',char(method_list(meth)),';']);%
    for j=1:size(method_list,2)
        eval(['base_spectres(t).combined_frac.',char(method_list(j)),' = Ybv_',char(method_list(j)),';']); % Somme des distributions de taille
    end
    eval(['SStot(:,:,t)=base_spectres(t).combined_frac.',char(method_list(3)),';']);
    for i=1:size(method_list)
        a=method_list(i);
        eval(['clear Ybv_',cell2mat(a)])
    end
    for i=1:m
        eval(['clear Ybv',zoo_groups{i,1}])
        eval(['clear Yab',zoo_groups{i,1}])
    end
end

SStot=permute(SStot,[1 3 2]);
for i=1:size(SStot,3)
    Sz=SStot(:,:,i);
    S(i,1)={Sz};
end

SStot=S;

% -------- Fermeture ------------
close(h)
