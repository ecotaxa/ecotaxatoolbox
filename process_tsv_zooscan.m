function [base_Zooscan,m,zoo_groups,Idlistshort]=process_tsv_zooscan



A=dir('*.tsv');
filenames={A.name};
[n,m]=size(filenames);

% sampleheaders=readtable('allsampleheaders.xlsx');
% scanheaders=readtable('allscanheaders.xls');
base_Zooscan=[];
Idlist=[];
Idlistshort=[];
h = waitbar(0,'Please wait...');

for i=1:m
    
    S=readtable(char(filenames(i)),'Filetype','text','ReadVariableNames',1);
    sample=unique(S.sample_id);
    
    
    %% cleaning anotation hierarchy
    [n,no_use]=size(S);
    
    for j=1:n
        f=find(S.object_annotation_hierarchy{j,:}=='-');
        S.object_annotation_hierarchy{j,1}(f)='_';
        
        f=find(S.object_annotation_hierarchy{j,:}=='>');
        S.object_annotation_hierarchy{j,1}(f)='_';
    end
    
    % I=strfind(sampleheaders.sampleid,sample);
    % index = (cellfun(@isempty,  I) ==0);
    % sample_header = sampleheaders(index,:);
    %
    % I=strfind(scanheaders.sampleid,sample);
    % index = (cellfun(@isempty,  I) ==0);
    % scan_header = scanheaders(index,:);
    
    %% if µ sign is used
%     f=find(char(sample)=='µ');
%     sample=char(sample);
%     sample(f:f+1)='mu';
%     sample=cellstr(sample);
    %%
    
    base_Zooscan(i).SampleID=sample;
    base_Zooscan(i).Scanop=unique(S.sample_scan_operator);
    %base_Zooscan(i).Idstatus=
    base_Zooscan(i).Ship=unique(S.sample_ship);
    base_Zooscan(i).Scientificprog=unique(S.sample_program);
    base_Zooscan(i).StationId=str2num(cell2mat(unique(S.sample_stationid)));
    base_Zooscan(i).Date=unique(S.object_date);
    base_Zooscan(i).time=unique(S.object_time);
    %base_Zooscan(i).Datenum=
%         temp=unique(S.object_lat); % to be use when several fractions     entered with round numbers
%     base_Zooscan(i).Latitude=temp(1);
    base_Zooscan(i).Latitude=unique(S.object_lat);
    
    base_Zooscan(i).Longitude=unique(S.object_lon);
    base_Zooscan(i).Depth=str2num(cell2mat(unique(S.sample_bottomdepth)));
    base_Zooscan(i).CTDref=unique(S.sample_ctdrosettefilename);
    base_Zooscan(i).Otherref=unique(S.sample_other_ref);
     base_Zooscan(i).barcode=unique(S.sample_barcode);
    base_Zooscan(i).Townb=str2num(cell2mat(unique(S.sample_tow_nb)));
    base_Zooscan(i).Towtype=str2num(cell2mat(unique(S.sample_tow_type)));
    base_Zooscan(i).Nettype=unique(S.sample_net_type);
    base_Zooscan(i).Netmesh=str2num(cell2mat(unique(S.sample_net_mesh)));
    base_Zooscan(i).Netsurf=str2num(cell2mat(unique(S.sample_net_surf)));
    base_Zooscan(i).Zmax=str2num(cell2mat(unique(S.sample_zmax)));
    base_Zooscan(i).Zmin=str2num(cell2mat(unique(S.sample_zmin)));
%     temp=unique(S.sample_tot_vol);  % to be use when several fractions     entered with round numbers
%      base_Zooscan(i).Vol=str2num(cell2mat(temp(1)));
    base_Zooscan(i).Vol=str2num(cell2mat(unique(S.sample_tot_vol)));
    base_Zooscan(i).Sample_comments=unique(S.sample_comment);
    
    
    %% getting pixel size in micrometer and converting in mm
    %unique(process_particle_pixel_size__m)
    pixelsize=str2num((cell2mat(unique(S.process_particle_pixel_size_mm))));  % in mm/pixel
    
    if isempty(pixelsize)==1
        pixelsize=str2num((cell2mat(unique(S.process_particle_pixel_size__m))));  % in mm/pixel
    end
    %%
    
    %base_Zooscan(i).tot.Scanfilename=
    base_Zooscan(i).FracIds=unique(S.acq_id);
    [nfrac, no_use]=size(base_Zooscan(i).FracIds);
    for fracnb=1:nfrac

    
    
        I=strcmp(S.acq_id,base_Zooscan(i).FracIds(fracnb));
        eval(['base_Zooscan(i).d' num2str(fracnb) '.Fracmin=unique(S.acq_min_mesh(I));']);
        eval(['base_Zooscan(i).d' num2str(fracnb) '.Fracsup=unique(S.acq_max_mesh(I));'])
        eval(['base_Zooscan(i).d' num2str(fracnb) '.Fracnb= unique(S.acq_sub_part(I));'])
        %base_Zooscan(i).tot.Scanned_Objects=
        eval(['base_Zooscan(i).d' num2str(fracnb) '.Resolution=unique(S.process_img_resolution(I));'])
        
        eval(['base_Zooscan(i).d' num2str(fracnb) '.object_annotation_hierarchy=unique(S.object_annotation_hierarchy(I));'])
        
        [Idlist, ia, ic]=unique([Idlist; S.object_annotation_hierarchy(I)]);
        
        Idlistshort=[Idlistshort; S.object_annotation_category(I)];
        Idlistshort=Idlistshort(ia);
        %base_Zooscan(i).tot.Comments=
        eval(['base_Zooscan(i).d' num2str(fracnb) '.object_annotation_hierarchy=S.object_annotation_hierarchy(I);'])
                    I2=strfind(S.object_annotation_status(I),'validated');
    index = (cellfun(@isempty,  I2) ==0);
    
    eval(['base_Zooscan(i).d' num2str(fracnb) '.percentValidated=100*sum(index)/length(index);'])
        eval(['base_Zooscan(i).d' num2str(fracnb) '.major=S.object_major(I)*pixelsize;']) %object_perimmajor
        eval(['base_Zooscan(i).d' num2str(fracnb) '.minor=S.object_minor(I)*pixelsize;'])
        eval(['base_Zooscan(i).d' num2str(fracnb) '.area_exc=S.object_area_exc(I)*(pixelsize^2);'])
        eval(['base_Zooscan(i).d' num2str(fracnb) '.area=S.object_area(I)*(pixelsize^2);'])    %object__area
        eval(['base_Zooscan(i).d' num2str(fracnb) '.ESD=2*(((S.object_area(I)*(pixelsize^2))/pi).^0.5);'])    %object__area
        eval(['base_Zooscan(i).d' num2str(fracnb) '.perimferet=S.object_feret(I)*pixelsize;'])   % object_perimferet   %object_feretareaexc
        eval(['base_Zooscan(i).d' num2str(fracnb) '.conver=  str2num(cell2mat(base_Zooscan(i).d' num2str(fracnb) '.Fracnb))./base_Zooscan(i).Vol;']) %per cubic meter
        
    end
    
    
waitbar(i/m)
end
close(h)

[n,m]=size(base_Zooscan);
zoo_groups=Idlist;