% fonction de groupage du zooplankton pour calculer abondance et spectres
% des PFT

%variables d'entrée:
% 1 table_groupage : un tableau xls avec 3 colonnes. 
%       1iere colonne = groupes originaux
%       2ieme colonne = nouveaux noms des groupes living vs non_living, 
%       3ieme colonne = nouveaux noms des groupes PFT
% 2 abundances : tableau d'abondances (ind/m3)
%       # lignes = # de groupes de zoo identifiés à la base 
%       # colonnes = # d'échantillons
% 3 biovolumes : tableau de biovolumes (mm3/m3)
%       # lignes = # de groupes de zoo identifiés à la base 
%       # colonnes = # d'échantillons
% 4 size_spectra : SStot en sortie de la routine de creation de la
% base_spectres

function [base_regroup] = f_regroup_all(table_groupage,base)

tab=table_groupage;
%tab_orig=id;

tab_orig=base.Zoo_groups;
id=tab_orig;

g_pft2=tab(:,3);
%g_trophic=tab(:,4);

% not_to_keep = strfind(id,'not_living');
% to_keep = (cellfun(@isempty,not_to_keep)==1);

to_keep=[];
for i=1:size(id,1)
    a=strmatch(id(i),tab(:,1),'exact');
    if isempty(a)==0
        if length(a)>1
           warning(['warning' id(i) 'is in double, please check reference excell file before going further']);
        end
        to_keep=[to_keep;a];
    end
end

tab=tab(to_keep,:);

g_orig=tab(:,1);
g_living=tab(:,2);% idée: recherche des strings "living" et "not-living"
g_pft=tab(:,3);
g_trophic=tab(:,4);


%données de spectres de taille
% Ab
% Yab
% Ybv_Plain_Area_BV_spectra
% Ybv_Riddled_Area_BV_spectra
% Bv
% Ybv_Ellipsoid_BV_spectra

Ab=base.Ab;
%Yab=base.Yab;
Ybv_Plain_Area_BV_spectra=base.Ybv_Plain_Area_BV_spectra;
Ybv_Riddled_Area_BV_spectra=base.Ybv_Riddled_Area_BV_spectra;
Bv=base.Bv;
ss=base.Ybv_Ellipsoid_BV_spectra;


clear SA
%% regroupement total = 'all'

Ab_t=sum(Ab);
%Yab_t=sum(Yab,2);
Ybv_Plain_Area_BV_spectra_t=sum(Ybv_Plain_Area_BV_spectra,2);
Ybv_Riddled_Area_BV_spectra_t=sum(Ybv_Riddled_Area_BV_spectra,2);
Bv_t=sum(Bv);

sst=sum(ss,2);



%% regroupement living vs non living
new_groups=unique(g_living);
for i=1:size(new_groups,1)
    ng=new_groups(i);
    f_ng=strmatch(ng,g_living);

    
    % spectres bv
    ss_r=ss(:,f_ng);
    Ab_r=Ab(f_ng);
    %Yab_r=Yab(:,f_ng);
    Ybv_Plain_Area_BV_spectra_r=Ybv_Plain_Area_BV_spectra(:,f_ng);
    Ybv_Riddled_Area_BV_spectra_r=Ybv_Riddled_Area_BV_spectra(:,f_ng);
    Bv_r=Bv(f_ng);
    
    Ab_n1(i)=sum(Ab_r);
    %Yab_n1(:,i)=sum(Yab_r,2);
    Ybv_Plain_Area_BV_spectra_n1(:,i)=sum(Ybv_Plain_Area_BV_spectra_r,2);
    Ybv_Riddled_Area_BV_spectra_n1(:,i)=sum(Ybv_Riddled_Area_BV_spectra_r,2);
    Bv_n1(i)=sum(Bv_r);
    ssn1(:,i)=sum( ss_r,2);
    
    
    clear SSJ ssj ss_r f_ng ng Bv_r Ybv_Riddled_Area_BV_spectra_r Ybv_Plain_Area_BV_spectra_r  Ab_r
end

%% regroupement pft
new_groups=unique(g_pft2);
for i=1:size(new_groups,1)
    ng=new_groups(i);
    f_ng=strmatch(ng,g_pft);
   
    
    %spectres bv
    ss_r=ss(:,f_ng);
        Ab_r=Ab(f_ng);
    %Yab_r=Yab(:,f_ng);
    Ybv_Plain_Area_BV_spectra_r=Ybv_Plain_Area_BV_spectra(:,f_ng);
    Ybv_Riddled_Area_BV_spectra_r=Ybv_Riddled_Area_BV_spectra(:,f_ng);
    Bv_r=Bv(f_ng);
    
        Ab_n2(i)=sum(Ab_r);
    %Yab_n2(:,i)=sum(Yab_r,2);
    Ybv_Plain_Area_BV_spectra_n2(:,i)=sum(Ybv_Plain_Area_BV_spectra_r,2);
    Ybv_Riddled_Area_BV_spectra_n2(:,i)=sum(Ybv_Riddled_Area_BV_spectra_r,2);
    Bv_n2(i)=sum(Bv_r);
    ssn2(:,i)=sum( ss_r,2) ;
 
    clear SSJ ssj ss_r  f_ng ng Bv_r Ybv_Riddled_Area_BV_spectra_r Ybv_Plain_Area_BV_spectra_r Ab_r
end
    
    %% regroupement trophic
new_groups=unique(cell2mat(g_trophic));
for i=1:size(new_groups,1)
    ng=new_groups(i);
    f_ng=ng==cell2mat(g_trophic);
   
    
    %spectres bv
    ss_r=ss(:,f_ng);
            Ab_r=Ab(f_ng);
    %Yab_r=Yab(:,f_ng);
    Ybv_Plain_Area_BV_spectra_r=Ybv_Plain_Area_BV_spectra(:,f_ng);
    Ybv_Riddled_Area_BV_spectra_r=Ybv_Riddled_Area_BV_spectra(:,f_ng);
    Bv_r=Bv(f_ng);
    
        Ab_n3(i)=sum(Ab_r);
    %Yab_n3(:,i)=sum(Yab_r,2);
    Ybv_Plain_Area_BV_spectra_n3(:,i)=sum(Ybv_Plain_Area_BV_spectra_r,2);
    Ybv_Riddled_Area_BV_spectra_n3(:,i)=sum(Ybv_Riddled_Area_BV_spectra_r,2);
    Bv_n3(i)=sum(Bv_r);
    ssn3(:,i)=sum( ss_r,2);
   
    
    clear SSJ ssj ss_r  f_ng ng Bv_r Ybv_Riddled_Area_BV_spectra_r Ybv_Plain_Area_BV_spectra_r Ab_r
end
%%

%Spec_regroup=[ss sst ssn1 ssn2 ssn3];

Ab_regroup=[Ab' Ab_t Ab_n1 Ab_n2 Ab_n3];
%Yab_regroup=[Yab Yab_t Yab_n1 Yab_n2 Yab_n3];
Ybv_Plain_Area_BV_spectra_regroup=[Ybv_Plain_Area_BV_spectra Ybv_Plain_Area_BV_spectra_t Ybv_Plain_Area_BV_spectra_n1 Ybv_Plain_Area_BV_spectra_n2 Ybv_Plain_Area_BV_spectra_n3];
Ybv_Riddled_Area_BV_spectra_regroup=[Ybv_Riddled_Area_BV_spectra Ybv_Riddled_Area_BV_spectra_t Ybv_Riddled_Area_BV_spectra_n1 Ybv_Riddled_Area_BV_spectra_n2 Ybv_Riddled_Area_BV_spectra_n3];
Bv_regroup=[Bv' Bv_t Bv_n1 Bv_n2 Bv_n3];
Ybv_Ellipsoid_BV_spectra_regroup=[ss sst ssn1 ssn2 ssn3];
base_regroup=[];
%base_regroup.Abtot=base.Abtot; % ne sert a rien
base_regroup.Ab=Ab_regroup;
%base_regroup.Yab=Yab_regroup;
base_regroup.Ybv_Plain_Area_BV_spectra=Ybv_Plain_Area_BV_spectra_regroup;
base_regroup.Ybv_Riddled_Area_BV_spectra=Ybv_Riddled_Area_BV_spectra_regroup;
base_regroup.Bv=Bv_regroup;
%base_regroup.Bvtot=base.Bvtot; % ne sert a rien
base_regroup.Ybv_Ellipsoid_BV_spectra=Ybv_Ellipsoid_BV_spectra_regroup;

base_regroup.Zoo_groups=[g_orig;{'all'};unique(g_living);unique(g_pft2);cellstr(num2str(unique(cell2mat(g_trophic))))];
base_regroup.originalplace=[1 length(g_orig)];
base_regroup.allplace=length(g_orig)+1;
base_regroup.livingplace=length(g_orig)+2;
base_regroup.notlivingplace=length(g_orig)+3;
base_regroup.pftplace=[length(g_orig)+4 length(g_orig)+3+length(unique(g_pft2))];
base_regroup.trophicplace=[length(g_orig)+4+length(unique(g_pft2)) length(g_orig)+3+length(unique(g_pft2))+length(unique(cell2mat(g_trophic)))];
% groups=[g_orig(:,1);{'all'};unique(g_living);unique(g_pft)];
