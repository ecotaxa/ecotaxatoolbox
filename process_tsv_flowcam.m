function [base_Zooscan,m,zoo_groups,Idlistshort]=process_tsv_flowcam

answer = questdlg('Flowcam could be done using one concentration/dilution step of the collector (entered as volpump) is it your case? (caution, for the moment it convert back any volpump>10 to a concentration factor of 1)', ...
    'concentration step', ...
    'Yes','No, no concentration have been done','No, no concentration have been done');

answer2 = questdlg('Flowcam images could be separated in two groups by size with a different subsampling (needs to be in two different folder), is it the case?', ...
    'concentration step', ...
    'Yes','No','No');

switch answer2
    case 'Yes'
        f = msgbox('Please select the folder containing the *.tsv files of the first fraction');
        uiwait(f)
        folder1=uigetdir;
        
        prompt = {'Fractionation operated? (in %)'};
        title = 'fraction1';
        dims = [1 35];
        definput = {'100'};
        answerfrac1 = inputdlg(prompt,title,dims,definput);
        answerfrac1=str2num(cell2mat(answerfrac1))/100;
        
        f = msgbox('Please select the folder containing the *.tsv files of the second fraction');
        uiwait(f)
        folder2=uigetdir;
        prompt = {'Fractionation operated? (in %)'};
        title = 'fraction1';
        dims = [1 35];
        definput = {'30'};
        answerfrac2 = inputdlg(prompt,title,dims,definput);
        answerfrac2=str2num(cell2mat(answerfrac2))/100;
        
        cd(folder1)
        A=dir('*.tsv');
        filenames={A.name};
        [n,m]=size(filenames);
        
        %% to check that the same number of tsv apre present in both folders
    case 'No'
        A=dir('*.tsv');
        filenames={A.name};
        [n,m]=size(filenames);
        
end




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
    switch answer2
        case 'Yes'
            cd(folder1)
            S=readtable(char(filenames(i)),'Filetype','text','ReadVariableNames',1);
            S.acq_id(:)={'d1'};
            cd(folder2)
            S2=readtable(char(filenames(i)),'Filetype','text','ReadVariableNames',1);
            S2.acq_id(:)={'d2'};
            S=[S;S2];
           
        case 'No'
            
            
            S=readtable(char(filenames(i)),'Filetype','text','ReadVariableNames',1);
    end
    sample=unique(S.sample_id);
    
    %return
    %% cleaning empty cells
    
    [no_use,n]=size(S);
    
    
    varnames=S.Properties.VariableNames;
    for j=21:n
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
    base_Zooscan(i).raw_image_total=str2num(cell2mat(unique(S.acq_raw_image_total)));
     base_Zooscan(i).process_nb_images=str2num(cell2mat(unique(S.process_nb_images)));
    
    base_Zooscan(i).Volconc=str2num(cell2mat(unique(S.sample_volconc)))/1000000;
    base_Zooscan(i).Volpump=str2num(cell2mat(unique(S.sample_volpump)));
    
    
    switch answer
        case 'Yes'
            %% correcting the first leg which volpum is afterward the factor of concentration/dillution (passing everything in concentration factor)
            if base_Zooscan(i).Volpump>10
                base_Zooscan(i).Volpump=1;
            end
        case 'No, no concentration have been done'
    end
    
    
    %%
    test=unique(S.sample_comment_or_volume);
    if length(test)==1
        base_Zooscan(i).comment_or_volume=str2num(cell2mat(unique(S.sample_comment_or_volume)))/1000000;
    else
        base_Zooscan(i).comment_or_volume=str2num(cell2mat(unique(S.sample_comment)))/1000000;
    end
    
    if isempty(base_Zooscan(i).comment_or_volume)==1
        base_Zooscan(i).comment_or_volume=1;
    end
    
    temp=unique(S.acq_fluid_volume_imaged);
    if length(test)==1
    f=strfind(temp,'_');   %supressing the _ml
    temp=char(temp);
    temp2=temp(1:cell2mat(f)-1);
    
    base_Zooscan(i).fluidimaged=str2num(temp2)/1000000;
    else
        char(filenames(i))
    %% correcting for inomogeneous notations of metadata in TPac part 1 (to supress when corrected)
     base_Zooscan(i).fluidimaged=0.1863*base_Zooscan(i).Volpump;
    end

    base_Zooscan(i).Volume_Imaged_Processed= (base_Zooscan(i).fluidimaged/ base_Zooscan(i).raw_image_total)* min(base_Zooscan(i).raw_image_total,base_Zooscan(i).process_nb_images);
    
    
    
    
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
        eval(['base_Zooscan(i).d' num2str(fracnb) '.Fracmin=unique(S.acq_min_esd);']);
        eval(['base_Zooscan(i).d' num2str(fracnb) '.Fracsup=unique(S.acq_max_esd);']);
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
            eval(['base_Zooscan(i).d' num2str(fracnb) '.perimferet=S.object_feret(I)*pixelsize;']);   % object_perimferet   %object_feretareaexc
        else
            eval(['base_Zooscan(i).d' num2str(fracnb) '.major=cellfun(@str2num,S.object_major(I))*pixelsize;']); %object_perimmajor
            eval(['base_Zooscan(i).d' num2str(fracnb) '.minor=cellfun(@str2num,S.object_minor(I))*pixelsize;']);
            eval(['base_Zooscan(i).d' num2str(fracnb) '.area_exc=cellfun(@str2num,S.object_area_exc(I))*(pixelsize^2);']);
            eval(['base_Zooscan(i).d' num2str(fracnb) '.ESD=2*((cellfun(@str2num,S.object_area(I))*(pixelsize^2)/pi).^0.5);']);
            eval(['base_Zooscan(i).d' num2str(fracnb) '.area=cellfun(@str2num,S.object_area(I))*(pixelsize^2);']);    %object__area
            eval(['base_Zooscan(i).d' num2str(fracnb) '.area_origin=cellfun(@str2num,S.object_area(I));']);    %object__area
            eval(['base_Zooscan(i).d' num2str(fracnb) '.perimferet=cellfun(@str2num,S.object_feret(I))*pixelsize;']);   % object_perimferet   %object_feretareaexc
        end
        
        switch answer2
            case 'Yes'
                switch answer
                    case 'Yes'
                        %% (temporary correction) correcting the first leg which volpum is afterward the factor of concentration/dillution (passing everything in concentration factor)
                        if base_Zooscan(i).Volpump>10
                            base_Zooscan(i).Volpump=1;
                        end
                        
                        if fracnb==1
                            falsemotoda=answerfrac1;
                        elseif fracnb==2
                            falsemotoda=answerfrac2;
                        end
                        
                        %% to adjust depending on the project (not all have been marked the same way)
                        
                        eval(['base_Zooscan(i).d' num2str(fracnb) '.conver=  base_Zooscan(i).Volpump.*base_Zooscan(i).Volconc./(base_Zooscan(i).Volume_Imaged_Processed.*base_Zooscan(i).comment_or_volume.*falsemotoda);']); %ind per cubic meter
                    case 'No, no concentration have been done'
                        eval(['base_Zooscan(i).d' num2str(fracnb) '.conver=  base_Zooscan(i).Volconc./(base_Zooscan(i).Volume_Imaged_Processed.*base_Zooscan(i).comment_or_volume.*falsemotoda);']); %ind per cubic meter
                end
            case 'No'
                switch answer
                    case 'Yes'
                        %% (temporary correction) correcting the first leg which volpum is afterward the factor of concentration/dillution (passing everything in concentration factor)
                        if base_Zooscan(i).Volpump>10
                            base_Zooscan(i).Volpump=1;
                        end
                        %% to adjust depending on the project (not all have been marked the same way)
                        eval(['base_Zooscan(i).d' num2str(fracnb) '.conver=  base_Zooscan(i).Volpump.*base_Zooscan(i).Volconc./(base_Zooscan(i).Volume_Imaged_Processed.*base_Zooscan(i).comment_or_volume);']); %ind per cubic meter
                    case 'No, no concentration have been done'
                        eval(['base_Zooscan(i).d' num2str(fracnb) '.conver=  base_Zooscan(i).Volconc./(base_Zooscan(i).Volume_Imaged_Processed.*base_Zooscan(i).comment_or_volume);']); %ind per cubic meter
                end
        end
        
    end
    
    waitbar(i/m)
end
close(h)


answer = questdlg('is this Tara Polar circle data', ...
    'concentration step', ...
    'Yes','No','No');

switch answer
    case 'Yes'
        %% correcting from logsheets sample 158 in tara having missing metadata
        
        base_Zooscan(22).comment_or_volume=50*10000000/1500;
        base_Zooscan(23).comment_or_volume=50*6000000/1500;
        base_Zooscan(24).comment_or_volume=50*10000000/1500;
        base_Zooscan(22).d1.conver=  base_Zooscan(22).Volconc./(base_Zooscan(22).Volpump.*base_Zooscan(22).comment_or_volume); %ind per cubic meter
        base_Zooscan(23).d1.conver=  base_Zooscan(23).Volconc./(base_Zooscan(23).Volpump.*base_Zooscan(23).comment_or_volume); %ind per cubic meter
        base_Zooscan(24).d1.conver=  base_Zooscan(24).Volconc./(base_Zooscan(24).Volpump.*base_Zooscan(24).comment_or_volume); %ind per cubic meter
        
end


[n,m]=size(base_Zooscan);
zoo_groups=Idlist;