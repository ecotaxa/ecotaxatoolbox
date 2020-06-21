%% Construct a matlab base structured in size for a given depth strata (UVP)
% by Fabien Lombard 2016-2018

% clear all
% close all

%% offset_depth : 
%currently ecotaxa does not implement the offset dept of 1.2 m between the 
%depth sensor and the imaged zone
%this will be corrected in future versions of ecotaxa (and thus depth
%offset would have to be fixed to 0m

offset_depth=1.2;
%%

%for i = 1:width(t), if iscell(t.(i)), t.(i) = cell2mat(t.(i)); end, end

prompt = {'Enter max depth'};%
dlg_title = 'Input';
num_lines = 1;
defaultans = {'200'};
maxdepth = inputdlg(prompt,dlg_title,num_lines,defaultans);
maxdepth=str2num(maxdepth{1});


A=dir('*tsv');
filenames={A.name};
[n,m]=size(filenames);
f = msgbox('select the particle base')
[file,path] = uigetfile('*.mat')

load(file)  % load tha particles data base in which volume per depth intervals are storred
samplebase={base(:).profile};
%samplebase=samplebase';
%% caution this assumes that the same UVP was used for an entire project (not anymore= corrected )
expo=[base(:).exp];	
aa=[base(:).aa];
%pixelsize=aa.^(expo); %
pixelsize=[base(:).PixelSize]; %in mm
%pixelsize=pixelsize.^0.5; % in mm2


% sampleheaders=readtable('allsampleheaders.xlsx');
% scanheaders=readtable('allscanheaders.xls');
base_Zooscan=[];
Idlist=[];
%return
h = waitbar(0,'Please wait...');

%return

%filenames([425 466])=[];

[n,m]=size(filenames);
%return

for i=1:m
    
%     if i==430
%         continue
%     end
%     if i==425
%         continue
%     end
%     if i==466
%         continue
%     end
    
    S=readtable(char(filenames(i)),'Filetype','text','ReadVariableNames',1);
    sample=unique(S.sample_id);
    %for i = 1:width(S), if iscell(t.(i)), t.(i) = cell2mat(t.(i)); end, end
    
    %% cleaning anotation hierarchy
    [n,no_use]=size(S);
    
    for j=1:n
        f=find(S.object_annotation_hierarchy{j,:}=='-');
        S.object_annotation_hierarchy{j,1}(f)='_';
        
        f=find(S.object_annotation_hierarchy{j,:}=='>');
        S.object_annotation_hierarchy{j,1}(f)='_';
    end
    
    
    %% cleaning for messy (text) entrance of files
    
    if iscell(S.object_depth_max)==1
    S.object_depth_max=cellfun(@str2num,S.object_depth_max);
    end

    if iscell(S.object_major)==1
    S.object_major=cellfun(@str2num,S.object_major);
    end

    if iscell(S.object_minor)==1
    S.object_minor=cellfun(@str2num,S.object_minor);
    end

    if iscell(S.object_area_exc)==1
    S.object_area_exc=cellfun(@str2num,S.object_area_exc);
    end

    if iscell(S.object_area)==1
    S.object_area=cellfun(@str2num,S.object_area);
    end


    if iscell(S.object_feret)==1
    S.object_feret=cellfun(@str2num,S.object_feret);
    end

    
    
    base_Zooscan(i).SampleID=sample;
    base_Zooscan(i).DN=unique(S.sample_dn);
    %base_Zooscan(i).Idstatus=
    base_Zooscan(i).Ship=unique(S.sample_ship);
    base_Zooscan(i).Scientificprog=unique(S.sample_cruise);
    base_Zooscan(i).StationId=sample;
    temp=cell2mat(sample);
    %base_Zooscan(i).StationIdnum=num2str(temp(6:8));  % to extract number in tara
      base_Zooscan(i).StationIdnum=(temp);
    base_Zooscan(i).Date=unique(S.object_date);
    base_Zooscan(i).time=unique(S.object_time);
    %base_Zooscan(i).Datenum=
    base_Zooscan(i).Latitude=unique(S.object_lat);
    base_Zooscan(i).Longitude=unique(S.object_lon);
    base_Zooscan(i).Depth=str2num(cell2mat(unique(S.sample_bottomdepth)));
    base_Zooscan(i).CTDref=unique(S.sample_ctdrosettefilename);
    base_Zooscan(i).profileid=unique(S.sample_profileid);
%     base_Zooscan(i).Townb=str2num(cell2mat(unique(S.sample_tow_nb)));
%     base_Zooscan(i).Towtype=str2num(cell2mat(unique(S.sample_tow_type)));
%     base_Zooscan(i).Nettype=unique(S.sample_net_type);
%     base_Zooscan(i).Netmesh=str2num(cell2mat(unique(S.sample_net_mesh)));
%     base_Zooscan(i).Netsurf=str2num(cell2mat(unique(S.sample_net_surf)));
%     base_Zooscan(i).Zmax=str2num(cell2mat(unique(S.sample_zmax)));
%     base_Zooscan(i).Zmin=str2num(cell2mat(unique(S.sample_zmin)));
%     base_Zooscan(i).Vol=str2num(cell2mat(unique(S.sample_tot_vol)));
%     base_Zooscan(i).Sample_comments=unique(S.sample_comment);
    
    
    %% getting pixel size in micrometer and converting in mm
    %unique(process_particle_pixel_size__m)
    %pixelsize=(str2num((cell2mat(unique(S.process_pixel))))).^(0.5);  % in mm/pixel
    %base_Zooscan(i).pixelsize=pixelsize;
    %pixelsize=0.174^(0.5);
    %%
            test=strfind(samplebase,temp);
index = (cellfun(@isempty,  test)==0);
I=find(index==1);
if length(I)>1
    I=I(1);
end
     volumes=base(I).histnb.data.SampledVolume_L_;   %basetot(1,I).hisnb(:,3).*basetot(1,I).volimg0;   
     depthUVP=base(I).histnb.data.Depth_m_; %basetot(1,I).hisnb(:,1);
%      if isempty(basetot(I).zoopuvp5)==0;
%     
% pixelsize=basetot(I).zoopuvp5.pixel;
% 
%      end
     base_Zooscan(i).pixelsize=pixelsize(I); 
     
     I=depthUVP<(maxdepth+2.5); %because layer of UVP particles counts are centered: eg 200 = 197.5-202.5
     
     %return

     base_Zooscan(i).tot.vol=(volumes(I)/1000); %in m3
    base_Zooscan(i).tot.conver=1./(volumes(I)/1000);
    base_Zooscan(i).tot.depthstrata=depthUVP; %in m3

base_Zooscan(i).tot.totvol=nansum(volumes(I)/1000); %in m3
base_Zooscan(i).tot.totconver=1/nansum(volumes(I)/1000);
%return
 
    %base_Zooscan(i).tot.Scanfilename=
%     base_Zooscan(i).FracIds=unique(S.acq_id);
%     [nfrac, no_use]=size(base_Zooscan(i).FracIds);

        S.object_depth_max=S.object_depth_max+offset_depth;

        I=S.object_depth_max<(maxdepth+2.5);  %because layer of UVP particles counts are centered: eg 200 = 197.5-202.5
        %base_Zooscan(i).tot.Fracmin=unique(S.acq_min_mesh(I));
        %base_Zooscan(i).tot.Fracsup=unique(S.acq_max_mesh(I));
        %base_Zooscan(i).tot.Fracnb= unique(S.acq_sub_part(I));
        %base_Zooscan(i).tot.Scanned_Objects=
        %base_Zooscan(i).tot.Resolution=unique(S.process_img_resolution(I));
        
        base_Zooscan(i).tot.object_annotation_hierarchy=unique(S.object_annotation_hierarchy(I));
        Idlist=unique([Idlist; S.object_annotation_hierarchy(I)]);
        base_Zooscan(i).tot.object_annotation_hierarchy=S.object_annotation_hierarchy(I);
        base_Zooscan(i).tot.depth=S.object_depth_max(I);%,S.object_depth_max
        %base_Zooscan(i).tot.object_annotation_hierarchy=S.object_annotation_hierarchy;
        base_Zooscan(i).tot.major=S.object_major(I)*base_Zooscan(i).pixelsize; %object_perimmajor
        base_Zooscan(i).tot.minor=S.object_minor(I)*base_Zooscan(i).pixelsize;
        base_Zooscan(i).tot.area_exc=S.object_area_exc(I)*(base_Zooscan(i).pixelsize^2);
        base_Zooscan(i).tot.area=S.object_area(I)*(base_Zooscan(i).pixelsize^2);    %object__area
        base_Zooscan(i).tot.perimferet=S.object_feret(I)*base_Zooscan(i).pixelsize;   % object_perimferet   %object_feretareaexc
        %base_Zooscan(i).tot.conver=  str2num(cell2mat(base_Zooscan(i).d.Fracnb))./base_Zooscan(i).Vol; %per cubic meter
        
        

        

    
    
waitbar(i/m)
end
close(h)

%save base_temporary base_Zooscan
 %load base_temporary
%return
%% now working on spectra
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
[n,m]=size(base_Zooscan);
zoo_groups=Idlist;
%return
for i=1:m
    %[nfrac, no_use]=size(base_Zooscan(i).FracIds);
        subset=base_Zooscan(i).tot;
        fracnb=1;
        conver1=subset.conver;
        depth=subset.depth;
        depthstrata=subset.depthstrata;
        
        [ndepth,no_use]=size(conver1);
        depthstrata=depthstrata(1:ndepth);
        [nobject,no_use]=size(depth);
        conver2=NaN*zeros(nobject,1);
        %% caution only works with 5m depth intervals
        toto=0.5+depthstrata/5;
        toto2=ceil(depth/5);
        I=find(toto2>ndepth);% for the few cases where depth observations sligthly overpass max depth (by max 2.5m...)
        toto2(I)=ndepth;
        conver2=conver1(toto2);
        base_Zooscan(i).tot.conver=conver2;
        base_Zooscan(i).tot.converorig=conver1;
        subset.conver=conver2;
        %%
        
        [base_Zooscan SStot] = process_abundances_spectres_multiples_UVP(subset,base_Zooscan,smin,smax,k,uu,zoo_groups,i,fracnb);
        
    
    
    
end


 %save base_temporary base_Zooscan
  %load base_temporary

%return

%% now correcting from the fact that the data are in ind.m3 for each strata observed.... and cumulated
% thus abundances/biovolumes needs to be divided by the number of depth
% strata accumulated
for i=1:m
    [ndepth,no_use]=size(base_Zooscan(i).tot.converorig);
    base_Zooscan(i).tot.Ab=base_Zooscan(i).tot.Ab/ndepth;
    base_Zooscan(i).tot.Yab=base_Zooscan(i).tot.Yab/ndepth;
    base_Zooscan(i).tot.Ybv_Plain_Area_BV_spectra=base_Zooscan(i).tot.Ybv_Plain_Area_BV_spectra/ndepth;
    base_Zooscan(i).tot.Ybv_Riddled_Area_BV_spectra=base_Zooscan(i).tot.Ybv_Riddled_Area_BV_spectra/ndepth;
    base_Zooscan(i).tot.Bv=base_Zooscan(i).tot.Bv/ndepth;
    base_Zooscan(i).tot.Ybv_Ellipsoid_BV_spectra=base_Zooscan(i).tot.Ybv_Ellipsoid_BV_spectra/ndepth;

end


    %% plus producing regrouped groups
    table_groupage=readtable('zooregroup_zooscan.xlsx','ReadVariableNames',false);  %all copoda as omnivorous
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
k = strfind(zoo_groups,'temporary_');
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
                
                base_Zooscan(i).tot.Yab(:,test2)=[];
                base_Zooscan(i).tot.Ybv_Plain_Area_BV_spectra (:,test2)=[];
                base_Zooscan(i).tot.Ybv_Riddled_Area_BV_spectra(:,test2)=[];
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
writetable(tosave,'zooregroup_zooscan.xlsx','WriteVariableNames',0);
cd(folder);
end
if istemporary==1
    table_groupage=[table_groupage; addontemp];
end



%% producing the regrouped groups
    
    for i=1:m
        
        [base_regroup] = f_regroup_all(table_groupage,base_Zooscan(i).tot);
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
    end
    
    
    
%% saving the final bases and resume files

instrument=char(list2(indx2));

prompt = {'instrument:','project:'};
title = 'save base under the name:';
dims = [1 35];
definput = {instrument,'pointB_Regent_1995_2019'};
answer = inputdlg(prompt,title,dims,definput)


save(['base_instrument_' char(answer(1)) '_' char(answer(2))],'base_Zooscan','-v7.3')
%save base_spectre_zooscan_regent_point_B base_Zooscan
%save base_spectre_flowcam_168b20 base_spectres

%Abtable=table(Ab_resume,'VariableNames',Zoo_groups,'RowNames',samplelist);    % do not work because the name of taxa are TOO LONG
Abtable=array2table(Ab_resume','VariableNames',samplelist,'rownames',Zoo_groups);
Bvtable=array2table(Bv_resume','VariableNames',samplelist,'rownames',Zoo_groups);
writetable(Abtable,'Abundance_resume.csv','WriteRowNames',true)
writetable(Bvtable,'Biovolume_resume.csv','WriteRowNames',true)
    
  
