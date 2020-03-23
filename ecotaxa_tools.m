%% Ecotaxa tool (F. Lombard 2019 lombard@obs-vlfr.fr)
% ecotaxa tool helps you to process raw results originating from quantitative
% imaging devices and originating from Ecotaxa (https://ecotaxa.obs-vlfr.fr/)
% ecotaxa tool includes the initial process of raw *.tsv files available from
% the export function in ecotaxa (Warning, please use the option "one file/sample")
%
% In the future ecotaxa tool will also includes several functions allowing
% the raw visialisation and raw analysis of your data. Those analysis are only
% of indicational value and don't dispense to further analyse some more
% specific features by building/refining the analysis yourself
%
% Zooscan analysis is supporting several scans per samples (and assembling
% them accordingly to the fractionation)
%
% the resulting structured base includes the following variables as final
% results
%     Ab  abundance per groups (ind.m-3)
%     Bv biovolume per groups (mm3.m-3)
%     Ybv_Plain_Area_BV_spectra NBSSbiovolume per groups per size class (mm3.mm-3.m-3)
%     Ybv_Riddled_Area_BV_spectra NBSSbiovolume per groups per size class (mm3.mm-3.m-3)
%     Ybv_Ellipsoid_BV_spectra NBSSbiovolume per groups per size class (mm3.mm-3.m-3)
%     X  is the Middle of each biovolume size class (caution it is in log here) (log(mm3)
%     X1 is the amplitude of each biovolume size class - used for de-normalizing NBSS (mm3)
%     ESD vector is the conversion of X in ESD (µm this time)
%     ESDquartilesmoyenne are the 5 25 50 75 95 % quartiles of ESD (+ mean and std) of each groups

%% Fabien Lombard February 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

close all
clear all


directoryoftoolbox='/Users/fabienlombard/Documents/travail/projets/ecotaxatoolbox';


list = {'process bases (classic-quantitative)','process bases for images features analysis (qualitative) - not yet available','analyse base - not yet available'};

% to be added in future: gross analysis of images properties

[indx,tf] = listdlg('ListString',list,'SelectionMode','single','ListSize',[300,150]);


%% classic process
if indx==1 
    f = msgbox('Please select the folder containing the *.tsv files or allready processed bases');
uiwait(f)
folder=uigetdir;

% in the future, make communication with the API here (if/when available) to extract tsv directly from ecotaxa

cd(folder)

A=dir('*.tsv');

a=length(A);
    if a==0
    warning('The folder does not contain any *.tsv files from ecotaxa, please extract your data before from the Ecotaxa software')
    return
    end
    
    
    
    list2 = {'Zooscan','Flowcam','IFCB','UVP step 1-from "ecopart-detailed export" (no size)','UVP step 2-add ecotaxa export (size and other features)','Accuri flowcytometry','eHFCM_H5','eHFCM_H20 initial big file','eHFCM_H20 by sample update','Planktoscope trial 2020'};
    
    
    [indx2,tf] = listdlg('PromptString','What instrument was used?','ListString',list2,'SelectionMode','single','ListSize',[300,150]);
    
    
    if indx2==5 
        
        B=dir('*.mat');

b=length(B);
if b==0
    warning('The folder does not contain any processed bases, please process your "detailled particles" bases first and save it to the same folder as your raw files')
    return
end
    end

    
    if indx2==1
        [base_Zooscan,m,zoo_groups,Idlistshort]=process_tsv_zooscan;
        process_spectra_and_save
    elseif indx2==2
        [base_Zooscan,m,zoo_groups,Idlistshort]=process_tsv_flowcam;
        process_spectra_and_save
    elseif indx2==3
        [base_Zooscan,m,zoo_groups,Idlistshort]=process_tsv_IFCB;
        process_spectra_and_saveIFCB
    elseif indx2==4
         ecotaxa_particles_base
    elseif indx2==5
        UVPread_persample
    elseif indx2==6
        [base_Zooscan,m,zoo_groups]=process_accuri;
    elseif indx2==7
        [base_Zooscan,m,zoo_groups,Idlistshort]=process_tsv_H5;
        process_spectra_and_saveH5
            elseif indx2==8
        [base_Zooscan,m,zoo_groups,Idlistshort]=process_tsv_H20;
        process_spectra_and_saveH20
        
        elseif indx2==9
        [base_Zooscan,m,zoo_groups,Idlistshort]=process_tsv_H20_bysample;
        process_spectra_and_saveH20
                elseif indx2==10
                    %return
        [base_Zooscan,m,zoo_groups,Idlistshort]=process_tsv_planktoscope;
        process_spectra_and_save
    end
    
    
end
%%



%% process images (not yet integrated)
if indx==2
    warning('feature not yet integrated, may be prone to bugs')
    % to be done:
    %to test the presence of images
    %to migrate all images to one single folder
    % to select the object properties to analyse / excluse
    %
    return
end


%% analyze bases
if indx==3 
    

B=dir('*.mat');

b=length(B);

if b==0
    warning('The folder does not contain any processed bases, please process your bases first')
    return
end

%% start the app here
%
%%

end




% %% starting from here this should become a function latter (activated
% 
% 
% %% reading the tsv files and starting to create the bases
% if indx==1
%     
%     list = {'Zooscan','Flowcam','IFCB','UVP-from "ecopart-detailed export" (no size)','UVP detailled export +ecotaxa export (size and other features)','Accuri flowcytometry'};
%     
%     
%     [indx2,tf] = listdlg('PromptString','What instrument was used?','ListString',list,'SelectionMode','single','ListSize',[300,150]);
%     
%     
%     if indx2==5 && b==0
%     warning('The folder does not contain any processed bases, please process your "detailled particles" bases first and save it to the same folder as your raw files')
%     return
%     end
% 
%     
%     if indx2==1
%         [base_Zooscan,m,zoo_groups,Idlistshort]=process_tsv_zooscan;
%         process_spectra_and_save
%     elseif indx2==2
%         [base_Zooscan,m,zoo_groups,Idlistshort]=process_tsv_flowcam;
%         process_spectra_and_save
%     elseif indx2==3
%         [base_Zooscan,m,zoo_groups,Idlistshort]=process_tsv_IFCB;
%         process_spectra_and_saveIFCB
%     elseif indx2==4
%          ecotaxa_particles_base
%     elseif indx2==5
%         UVPread_persample
%     elseif indx2==6
%         %[base_Zooscan,m,zoo_groups]=process_accuri;
%     end
%     
% end


