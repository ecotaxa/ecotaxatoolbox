% clear all
% close all

% tara pac S=readtable('export_201_20161213_1207.tsv','Filetype','text','ReadVariableNames',1);

A=dir('*Export_metadata_summary.tsv');
filenames={A.name};
metadata=readtable(char(filenames),'Filetype','text','ReadVariableNames',1);

filenames=metadata.PlanktonFilename;
[m,n]=size(filenames);

% sampleheaders=readtable('allsampleheaders.xlsx');
% scanheaders=readtable('allscanheaders.xls');

base=table2struct(metadata);
Idlist=[];
h = waitbar(0,'Please wait...');

for i=1:m
     test=char(filenames(i));
     switch test 
         case 'no data available'
             disp(test)
             
             otherwise
    
    S=readtable(char(filenames(i)),'Filetype','text','ReadVariableNames',1);
    [titi,toto]=size(S);
    if titi==0
        continue
    end
    
    zoonames=S.Properties.VariableDescriptions;
    zoonames=zoonames';
    [n,no_use]=size(zoonames);
for j=1:n
    f=strfind(zoonames{j,:},'''');
    zoonames{j,:}(f)=' ';
    zoonames{j,:} = strrep(zoonames{j,:},'Original column heading: ','');
end
f=strfind(zoonames,' living [# m-3] '); %living[#]
index = (cellfun(@isempty,  f) ==0);
place=find(index==1);
place1=place(1);

f=strfind(zoonames,' living biovolume [mm3 l-1] '); %living biovolume[ppm]
index = (cellfun(@isempty,  f) ==0);
place=find(index==1);
place2=place(1);

f=strfind(zoonames,' living avgesd [mm] '); %living avgesd[mm3]
index = (cellfun(@isempty,  f) ==0);
place=find(index==1);
place3=place(1);

%extracting the sampled volume
subtable1=S(:,[1:place1-1 place1:place2-1]);
[n,o]=size(subtable1);
vol=table2array(subtable1(:,5));
vol=vol*ones(1,o-5)/1000; %from liters to m3
%extracting zoo counts and normalizing them by sampled volume
subtable1(:,6:end)=array2table(table2array(subtable1(:,6:end))./vol);
base(i).zoo_nb_per_m3.data=subtable1;
base(i).zoo_nb_per_m3.original_labels=zoonames([1:place1-1 place1:place2-1]);
   

subtable2=S(:,[1:place1-1 place2:place3-1]);
base(i).zoo_mm3_per_m3.data=subtable2;
base(i).zoo_mm3_per_m3.original_labels=zoonames([1:place1-1 place2:place3-1]);


subtable3=S(:,[1:place1-1 place3:j]);
base(i).zoo_mean_bv_mm3.data=subtable3;
base(i).zoo_mean_bv_mm3.original_labels=zoonames([1:place1-1 place3:j]);
 waitbar(i/m)
     end
end
close(h)
filenames=metadata.ParticleFilename;
[m,n]=size(filenames);
h = waitbar(0,'Please wait...');
for i=1:m
     
    S=readtable(char(filenames(i)),'Filetype','text','ReadVariableNames',1);
    zoonames=S.Properties.VariableDescriptions;
    zoonames=zoonames';
    [n,no_use]=size(zoonames);
for j=1:n
    f=strfind(zoonames{j,:},'''');
    zoonames{j,:}(f)=' ';
    zoonames{j,:} = strrep(zoonames{j,:},'Original column heading: ','');
end

f=strfind(zoonames,'LPM (');
index = (cellfun(@isempty,  f) ==0);
place=find(index==1);
place1=place(1);

f=strfind(zoonames,'LPM biovolume');
index = (cellfun(@isempty,  f) ==0);
place=find(index==1);
place2=place(1);
place3=place(end)+1;



subtable1=S(:,[1:place1-1 place1:place2-1]);
[n,o]=size(subtable1);


base(i).histnb.data=subtable1;
base(i).histnb.original_labels=zoonames([1:place1-1 place1:place2-1]);
   

subtable2=S(:,[1:place1-1 place2:place3-1]);
base(i).histbv.data=subtable2;
base(i).histbv.original_labels=zoonames([1:place1-1 place2:place3-1]);


subtable3=S(:,[1:place1-1 place3:j]);
base(i).CTD.data=subtable3;
base(i).CTD.original_labels=zoonames([1:place1-1 place3:j]);
 waitbar(i/m)
end
close(h)

uisave('base','base_UVP_location')