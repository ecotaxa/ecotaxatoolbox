
%% now working on getting correct abundances/biovolumes and NBSS spectra calculations
%% -------------- setting for calculus of the size spectra ------------------------------
% Area, Minor and Major are in mm in base_zooscan
smin=0.000000000001;     %set the lower limit of the biovolume spectra that will be calculated
smax=10000;      %set the upper limit of the biovolume spectra that will be calculated
k=2^(1/4);      %set logarithmic base used to calculate bins of the size spectra
%according to platt & denman's theory (1978). scaling
%exponent set to 0.25.
uu=1; % change here the first size class for the regression (and put option = 1 for performing on the spectra from uu to end)
%zoo_groups = sort(unique(zoo_groups));
%% starting to process the base
%return
for i=1:m
    %[nfrac, no_use]=size(base_Zooscan(i).FracIds);
    nfrac=1;
    for fracnb=1:nfrac
        eval(['subset=base_Zooscan(i).d' num2str(fracnb) ';'])
        [base_Zooscan SStot] = process_abundances_biovolumes_NBSSspectra_IFCB(subset,base_Zooscan,smin,smax,k,uu,zoo_groups,i,fracnb);
    end
end


%% scanning the different fractions and assemblating them
for i=1:m
    %[nfrac, no_use]=size(base_Zooscan(i).FracIds);
    nfrac=1;
        Ab=[]; % abundance per groups (ind.m-3)
        Bv=[]; % biovolume per groups (mm3.m-3)
        Yab=[];  %% abundance per groups per size class (ind.m-3) only using Ellipsoid_BV_spectra biovolumes
        Ybv_Plain_Area_BV_spectra=[];  % NBSSbiovolume per groups per size class (mm3.mm-3.m-3)
        Ybv_Summed_plain_Area_BV_spectra=[]; % NBSSbiovolume per groups per size class (mm3.mm-3.m-3)
        Ybv_Ellipsoid_BV_spectra=[]; % NBSSbiovolume per groups per size class (mm3.mm-3.m-3)
        for fracnb=1:nfrac
        if fracnb==1
            eval(['Ab=base_Zooscan(i).d' num2str(fracnb) '.Ab;'])
            eval(['Bv=base_Zooscan(i).d' num2str(fracnb) '.Bv;'])
            %eval(['Yab=base_Zooscan(i).d' num2str(fracnb) '.Yab;'])
            eval(['Ybv_Plain_Area_BV_spectra=base_Zooscan(i).d' num2str(fracnb) '.Ybv_Plain_Area_BV_spectra;'])
            eval(['Ybv_Summed_plain_Area_BV_spectra=base_Zooscan(i).d' num2str(fracnb) '.Ybv_Summed_plain_Area_BV_spectra;'])
            eval(['Ybv_Ellipsoid_BV_spectra=base_Zooscan(i).d' num2str(fracnb) '.Ybv_Ellipsoid_BV_spectra;'])
            
            
        else
            eval(['Ab=nansum([Ab base_Zooscan(i).d' num2str(fracnb) '.Ab],2);'])
            eval(['Bv=nansum([Bv base_Zooscan(i).d' num2str(fracnb) '.Bv],2);'])
            %eval(['Yab=nansum(cat(3,Yab, base_Zooscan(i).d' num2str(fracnb) '.Yab),3);'])
            eval(['Ybv_Plain_Area_BV_spectra=nansum(cat(3,Ybv_Plain_Area_BV_spectra, base_Zooscan(i).d' num2str(fracnb) '.Ybv_Plain_Area_BV_spectra),3);'])
            eval(['Ybv_Summed_plain_Area_BV_spectra=nansum(cat(3,Ybv_Summed_plain_Area_BV_spectra, base_Zooscan(i).d' num2str(fracnb) '.Ybv_Summed_plain_Area_BV_spectra),3);'])
            eval(['Ybv_Ellipsoid_BV_spectra=nansum(cat(3,Ybv_Ellipsoid_BV_spectra, base_Zooscan(i).d' num2str(fracnb) '.Ybv_Ellipsoid_BV_spectra),3);'])
            
            
        end
    end
    
    base_Zooscan(i).tot.Ab=Ab;              % abundance per fraction rapportée au volume (#/m3)
    base_Zooscan(i).tot.Bv=Bv;               % abundance per fraction rapportée au volume (#/m3)
    
    %base_Zooscan(i).tot.Yab = Yab;
    base_Zooscan(i).tot.Ybv_Plain_Area_BV_spectra = Ybv_Plain_Area_BV_spectra;
    base_Zooscan(i).tot.Ybv_Summed_plain_Area_BV_spectra = Ybv_Summed_plain_Area_BV_spectra;
    base_Zooscan(i).tot.Ybv_Ellipsoid_BV_spectra = Ybv_Ellipsoid_BV_spectra;
    
    
eval(['base_Zooscan(i).tot.X = base_Zooscan(i).d' num2str(fracnb) '.X;'])
    eval(['base_Zooscan(i).tot.X1 = base_Zooscan(i).d' num2str(fracnb) '.X1;'])
    eval(['base_Zooscan(i).tot.Zoo_groups = base_Zooscan(i).d' num2str(fracnb) '.Zoo_groups;'])
    eval(['X = base_Zooscan(i).d' num2str(fracnb) '.X;'])
    sizelist=exp(X);%in mm3
    ESD=2*(sizelist*3/(4*pi)).^(1/3); % in mm
    ESD=ESD*1000; %in Âµm
    base_Zooscan(i).tot.ESDvector=ESD;
end


%% plus producing regrouped groups

%% loading default mapping
table_groupage=readtable('zooregroup_IFCB.xlsx','ReadVariableNames',false);  %all copoda as herbivorous
table_groupage=table2cell(table_groupage);

%% if needing to add new functional/trophic finction, mofify and add your desired within the excell file
% currently trophic groups includes
% -1= do not feed
% 1 phototrophs
% 1.5 mixotrophs
% 2 grazers
% 2.5 omnivorous
% 3 predators
% 3.5 unknown trophic group
% (note a 0.5 place is available for bacteria living from dissolved
% matter and potentially a 0 place is possible for viruses... but the
% "placement" is still subject to debate)

%% now working on regrouping taxa per functional/trophic groups
%
Zoo_groups=table_groupage(:,1);
[n,p]=size(zoo_groups);

%return
%% finding if temporary groups are used
istemporary=0;
k = strfind(zoo_groups,'temporary_t0');
test=sum(cell2mat(k));
test2=cellfun(@isempty,k);
test2=test2==0;

if  test>0
    answer = questdlg('Your files includes one or several temporary "t00X" categories. Do you have any "functional/trophic" mapping existing for those', ...
        'temporary categories mapping', ...
        'Yes please load them','No please create them','No please ignore them (not recommended)','Yes please load them');
    
    switch answer
        case 'No please ignore them (not recommended)'
            zoo_groups(test2)=[];
            [n,p]=size(zoo_groups);
            %% updating to remove temporary groups
            for i=1:m
                base_Zooscan(i).tot.Zoo_groups(test2)=[];
                base_Zooscan(i).tot.Ab(test2)=[];              % abundance per fraction rapportée au volume (#/m3)
                base_Zooscan(i).tot.Bv(test2)=[];               % abundance per fraction rapportée au volume (#/m3)
                
                %base_Zooscan(i).tot.Yab(:,test2)=[];
                base_Zooscan(i).tot.Ybv_Plain_Area_BV_spectra (:,test2)=[];
                base_Zooscan(i).tot.Ybv_Summed_plain_Area_BV_spectra(:,test2)=[];
                base_Zooscan(i).tot.Ybv_Ellipsoid_BV_spectra(:,test2)=[];
            end
        case 'Yes please load them'
            [file,path] = uigetfile('*.xlsx')
            addontemp=readtable([path file],'ReadVariableNames',false);  %all copoda as herbivorous
            addontemp=table2cell(addontemp);
            Zoo_groups=[Zoo_groups; addontemp(:,1)];
            istemporary=1;
        case 'No please create them'
            
            groups=table_groupage(:,2:end);
            [n,p]=size(groups);
            
            group1=unique(groups(:,1));
            group2=unique(groups(:,2));
            group3=unique(cellstr(num2str(cell2mat(groups(:,3)))));
            
            %return
            I=cellfun(@isempty,k);
            I=I==0;
            new_taxa=zoo_groups(I);
            [p]=length(new_taxa);
            newfunctional=[];
            
            for i=1:p
                
                settings = settingsdlg('Description', ['A new temporary taxonomic group have been found ' char(new_taxa(i))],...
                    'title' , 'New taxa functional mapping',...
                    'Alive' , group1 ,...
                    'functional group' , group2 , ...
                    'trophic group' , group3 ,...
                    'WindowWidth' , 800)
                newfunctional=[newfunctional settings];
            end
            newfunctional=struct2table(newfunctional)
            newfunctional(:,4)=[];
            %% updating the xls reference list
            %
            addontemp=[new_taxa table2cell(newfunctional)];
            tosave=array2table(addontemp);
            [file,path] = uiputfile('temporarymapping_instrument_net_location.xlsx');
            filename = fullfile(path,file);
            writetable(tosave,filename,'WriteVariableNames',0);
            Zoo_groups=[Zoo_groups; addontemp(:,1)];           
            istemporary=1;

    end

end
%% checking if no "new" groups are present
%
[n,p]=size(zoo_groups);
new_taxa={};p=0;

for i=1:n
    %J=strcmp(char(zoo_groups(i)),Zoo_groups(1:I-1,:));
    J=strcmp(char(zoo_groups(i)),Zoo_groups);
    if sum(J)==0
        
        p=p+1;
        new_taxa(p,1)=zoo_groups(i);
    end
end

clear Zoo_group


%% proposing a mapping for the new groups
%

groups=table_groupage(:,2:end);
[n,p]=size(groups);

group1=unique(groups(:,1));
group2=unique(groups(:,2));
group3=unique(cellstr(num2str(cell2mat(groups(:,3)))));

%return
[p]=length(new_taxa);
newfunctional=[];

for i=1:p
    
    settings = settingsdlg('Description', ['A new taxonomic group have been found ' char(new_taxa(i))],...
        'title' , 'New taxa functional mapping',...
        'Alive' , group1 ,...
        'functional group' , group2 , ...
        'trophic group' , group3 ,...
        'WindowWidth' , 800)
    newfunctional=[newfunctional settings];
end

if p>0
newfunctional=struct2table(newfunctional)
newfunctional(:,4)=[];
%% updating the xls reference list
%
addon=[new_taxa table2cell(newfunctional)];

table_groupage=[table_groupage; addon];


tosave=array2table(table_groupage);
cd(directoryoftoolbox);
writetable(tosave,'zooregroup_IFCB.xlsx','WriteVariableNames',0);
cd(folder);
end

if istemporary==1
    table_groupage=[table_groupage; addontemp];
end


%%

for i=1:m
    
    [base_regroup] = f_regroup_allIFCB(table_groupage,base_Zooscan(i).tot);
    base_Zooscan(i).regroupped=base_regroup;
end


%% preparing resume files on abundance / BV per taxa per sample

Ab_resume=[];
Bv_resume=[];
samplelist=[];

for i=1:m
    
    Ab_resume=[Ab_resume;base_Zooscan(i).regroupped.Ab];
    Bv_resume=[Bv_resume;base_Zooscan(i).regroupped.Bv];
    samplelist=[samplelist;base_Zooscan(i).SampleID];
    Zoo_groups=base_Zooscan(i).regroupped.Zoo_groups;
    base_Zooscan(i).regroupped.Idlistshort=Idlistshort;
end

%% saving the final bases and resume files
instrument=char(list2(indx2));
prompt = {'instrument:','project:'};
title = 'save base under the name:';
dims = [1 35];
definput = {instrument,'pointB_Regent_1995_2019'};
answer = inputdlg(prompt,title,dims,definput)


save(['base_' char(answer(1)) '_' char(answer(2))],'base_Zooscan','-v7.3')
%save base_spectre_zooscan_regent_point_B base_Zooscan
%save base_spectre_flowcam_168b20 base_spectres

%Abtable=table(Ab_resume,'VariableNames',Zoo_groups,'RowNames',samplelist);    % do not work because the name of taxa are TOO LONG
Abtable=array2table(Ab_resume','VariableNames',samplelist,'rownames',Zoo_groups);
Bvtable=array2table(Bv_resume','VariableNames',samplelist,'rownames',Zoo_groups);
writetable(Abtable,'Abundance_resume.csv','WriteRowNames',true)
writetable(Bvtable,'Biovolume_resume.csv','WriteRowNames',true)


