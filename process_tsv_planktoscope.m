function [base_Zooscan,m,zoo_groups,Idlistshort]=process_tsv_planktoscope




A=dir('*.tsv');
filenames={A.name};
[n,m]=size(filenames);





%filenames(31:33)=[]; % sample from easter island are in double

% sampleheaders=readtable('allsampleheaders.xlsx');
% scanheaders=readtable('allscanheaders.xls');


base_Zooscan=[];
Idlist=[];
Idlistshort=[];


h = waitbar(0,'Please wait...');

for i=1:m
    i
    %     if i==60
    %         continue
    %     end
    
    
    
    
    S=readtable(char(filenames(i)),'Filetype','text','ReadVariableNames',1);
    
    %sample=unique(S.sample_id);
    sample=unique(S.object_date);
    
    %return
    %% cleaning empty cells
    
    [no_use,n]=size(S);
    
    
    varnames=S.Properties.VariableNames;
    for j=16:n
        eval(['var_j=S.' char(varnames(j)) ';']);
        
        if iscellstr(var_j)
            index = (cellfun(@isempty,  var_j) ==1);
            if sum(index)>0
                var_j(index)={'NaN'};
                var_j=(cellfun(@str2num,  var_j));
                %var_j=cell2mat(var_j);
                eval(['S.' char(varnames(j)) '=var_j;']);
                
            end
            
        end
        
    end
    
    
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
    
    %% if samples contains µ
    % f=find(char(sample)=='µ');
    % sample=char(sample);
    % sample(f:f+1)='mu';
    % sample=cellstr(sample);
    
    %% getting relevant info
    
    base_Zooscan(i).SampleID=sample;
    %      i
    %      sample
    %
    %base_Zooscan(i).Scanop=unique(S.sample_scan_operator);
    %base_Zooscan(i).Idstatus=
    base_Zooscan(i).Ship=unique(S.sample_ship);
    %base_Zooscan(i).Scientificprog=unique(S.sample_program);
    %base_Zooscan(i).StationId=str2num(cell2mat(unique(S.sample_stationid)));
    base_Zooscan(i).Date=unique(S.object_date);
    base_Zooscan(i).time=unique(S.object_time);
    %base_Zooscan(i).Datenum=
    base_Zooscan(i).Latitude=unique(S.object_lat);
    base_Zooscan(i).Longitude=unique(S.object_lon);
    %base_Zooscan(i).Depth=str2num(cell2mat(unique(S.sample_bottomdepth)));
    %base_Zooscan(i).CTDref=unique(S.sample_ctdrosettefilename);
    %base_Zooscan(i).Otherref=unique(S.sample_other_ref);
    %base_Zooscan(i).Townb=str2num(cell2mat(unique(S.sample_tow_nb)));
    %base_Zooscan(i).Towtype=str2num(cell2mat(unique(S.sample_tow_type)));
    %base_Zooscan(i).Nettype=unique(S.sample_net_type);
    %base_Zooscan(i).Netmesh=str2num(cell2mat(unique(S.sample_net_mesh)));
    %base_Zooscan(i).Netsurf=str2num(cell2mat(unique(S.sample_net_surf)));
    base_Zooscan(i).Zmax=unique(S.object_depth_min);
    base_Zooscan(i).Zmin=unique(S.object_depth_min);
    
    base_Zooscan(i).Volconc=1/(10^6); %1L en m3
    base_Zooscan(i).Volpump=1;
    
    
    
    
    
    
    %%
    %     test=unique(S.sample_comment_or_volume);
    %     if length(test)==1
    %         base_Zooscan(i).comment_or_volume=str2num(cell2mat(unique(S.sample_comment_or_volume)))/1000000;
    %     else
    %         base_Zooscan(i).comment_or_volume=str2num(cell2mat(unique(S.sample_comment)))/1000000;
    %     end
    
    %   if isempty(base_Zooscan(i).comment_or_volume)==1
    base_Zooscan(i).comment_or_volume=1;
    %  end
    
    temp=unique(S.acq_fluid_volume_imaged);
    %if length(test)==1
    f=strfind(temp,'_');   %supressing the _ml
    temp=char(temp);
    temp2=temp(1:cell2mat(f)-1);
    
    base_Zooscan(i).fluidimaged=str2num(temp2)/1000000;
    %     else
    %         char(filenames(i))
    %     %% correcting for inomogeneous notations of metadata in TPac part 1 (to supress when corrected)
    %      base_Zooscan(i).fluidimaged=0.1863*base_Zooscan(i).Volpump;
    %     end
    
    
    
    
    
    %% getting pixel size in micrometer and converting in mm
    %unique(process_particle_pixel_size__m)
    pixelsize=str2num((cell2mat(unique(S.process_pixel))))/1000;  % in mm/pixel
    base_Zooscan(i).pixelsize=pixelsize;
    %%
    
    
    base_Zooscan(i).FracIds=unique(S.acq_id);
    [nfrac, no_use]=size(base_Zooscan(i).FracIds);
    for fracnb=1:nfrac
        I=strcmp(S.acq_id,base_Zooscan(i).FracIds(fracnb));
        %return
        %         eval(['base_Zooscan(i).d' num2str(fracnb) '.Fracmin=unique(S.acq_min_esd);']);
        %         eval(['base_Zooscan(i).d' num2str(fracnb) '.Fracsup=unique(S.acq_max_esd);']);
        %base_Zooscan(i).tot.Fracnb= unique(S.acq_sub_part);
        %base_Zooscan(i).tot.Scanned_Objects=
        %base_Zooscan(i).tot.Resolution=unique(S.process_img_resolution);
        
        eval(['base_Zooscan(i).d' num2str(fracnb) '.object_annotation_hierarchy=unique(S.object_annotation_hierarchy(I));']);
        
        %%
        eval(['base_Zooscan(i).d' num2str(fracnb) '.anotStatus=S.object_annotation_status(I);']);
        [Idlist, ia, ic]=unique([Idlist; S.object_annotation_hierarchy(I)]);
        Idlistshort=[Idlistshort; S.object_annotation_category(I)];
        Idlistshort=Idlistshort(ia);
        %base_Zooscan(i).tot.Comments=
        
        I2=strfind(S.object_annotation_status(I),'validated');
        index = (cellfun(@isempty,  I2) ==0);
        
        eval(['base_Zooscan(i).d' num2str(fracnb) '.percentValidated=100*sum(index)/length(index);']);
        
        eval(['base_Zooscan(i).d' num2str(fracnb) '.object_annotation_hierarchy=S.object_annotation_hierarchy(I);']);
        %cellfun(@str2num,S.object_major)*pixelsize;
        if iscell(S.object_major)==0
            eval(['base_Zooscan(i).d' num2str(fracnb) '.major=S.object_major(I)*pixelsize;']); %object_perimmajor
            eval(['base_Zooscan(i).d' num2str(fracnb) '.minor=S.object_minor(I)*pixelsize;']);
            eval(['base_Zooscan(i).d' num2str(fracnb) '.area_exc=S.object_area_exc(I)*(pixelsize^2);']);
            eval(['base_Zooscan(i).d' num2str(fracnb) '.area=S.object_area(I)*(pixelsize^2);']);    %object__area
            eval(['base_Zooscan(i).d' num2str(fracnb) '.area_origin=S.object_area(I);']);    %object__area
            eval(['base_Zooscan(i).d' num2str(fracnb) '.ESD=2*(((S.object_area(I))*(pixelsize^2)/pi).^0.5);']);
            %eval(['base_Zooscan(i).d' num2str(fracnb) '.perimferet=S.object_feret(I)*pixelsize;']);   % object_perimferet   %object_feretareaexc
        else
            eval(['base_Zooscan(i).d' num2str(fracnb) '.major=cellfun(@str2num,S.object_major(I))*pixelsize;']); %object_perimmajor
            eval(['base_Zooscan(i).d' num2str(fracnb) '.minor=cellfun(@str2num,S.object_minor(I))*pixelsize;']);
            eval(['base_Zooscan(i).d' num2str(fracnb) '.area_exc=cellfun(@str2num,S.object_area_exc(I))*(pixelsize^2);']);
            eval(['base_Zooscan(i).d' num2str(fracnb) '.ESD=2*((cellfun(@str2num,S.object_area(I))*(pixelsize^2)/pi).^0.5);']);
            eval(['base_Zooscan(i).d' num2str(fracnb) '.area=cellfun(@str2num,S.object_area(I))*(pixelsize^2);']);    %object__area
            eval(['base_Zooscan(i).d' num2str(fracnb) '.area_origin=cellfun(@str2num,S.object_area(I));']);    %object__area
            %eval(['base_Zooscan(i).d' num2str(fracnb) '.perimferet=cellfun(@str2num,S.object_feret(I))*pixelsize;']);   % object_perimferet   %object_feretareaexc
        end
        
        
        eval(['base_Zooscan(i).d' num2str(fracnb) '.conver=  base_Zooscan(i).Volconc./(base_Zooscan(i).fluidimaged.*base_Zooscan(i).comment_or_volume);']); %ind per cubic meter
    end
    
    
    
    
    waitbar(i/m)
end
close(h)





[n,m]=size(base_Zooscan);
zoo_groups=Idlist;