%% move project to a new folder

old = '/Volumes/Vision/MRI/Decoding/derivatives/fmriprep/';
new = '/Volumes/Vision/MRI/DecodingPublic/derivatives/fmriprep/';

sub = {'0201','0202','0204','0205','0206','0228','0229','0248','0903'};


for iSub = 1:9
    
    for ses = 1:4
        
        [iSub ses]
        
        nowDir = sprintf('%ssub-%s/ses-0%s/anat/',new,sub{iSub},num2str(ses));
        oldDir = sprintf('%ssub-%s/ses-0%s/anat/',old,sub{iSub},num2str(ses));
        
        if exist(oldDir, 'dir')
            
            if ~exist(nowDir, 'dir')
                mkdir(nowDir)
            end
            
            system(['rsync -av ' oldDir ' ' nowDir])
        end
        
        nowDir = sprintf('%ssub-%s/ses-0%s/func/',new,sub{iSub},num2str(ses));
        oldDir = sprintf('%ssub-%s/ses-0%s/func/',old,sub{iSub},num2str(ses));
        
        if ~exist(nowDir, 'dir')
            mkdir(nowDir)
        end
        
        system(['rsync -av ' oldDir '*fsaverage6* ' nowDir])
        system(['rsync -av ' oldDir '*space-T1w_desc-preproc_bold* ' nowDir])
        
    end
    
end

%% deface


bidsDir = '/Volumes/Vision/MRI/DecodingPublic/';

sub = {'0201','0202','0204','0205','0206','0228','0229','0248','0903'};
twoFolder = {'rawdata','derivatives/fmriprep'};

tic
for iSub = 1:9 %for every subject
    
    for ses = 1:4 %for every session
        
        [iSub ses]
        
        for whichfolder = 1:2 % switching between the rawdata folder and the fmriprep folder as both contains anat
            
            %find the file with face
            input1 = dir(sprintf('%s%s/sub-%s/ses-0%s/anat/*w.nii.gz',bidsDir,twoFolder{whichfolder},sub{iSub},num2str(ses)));
            
            for ii = 1:length(input1) % for every face file
                f1 = [input1(ii).folder '/' input1(ii).name]; % get file name
                system(['/Users/pw1246/opt/anaconda3/bin/pydeface ' f1]) % deface
                f2 = [f1(1:end-7) '_defaced.nii.gz']; % find face file
                % if deface file exists, delete face file
                if exist(f2,'file')
                    delete(f1);
                end
            end
        end
        
        %remove mni files
        mni = dir(sprintf('%sderivatives/fmriprep/sub-%s/ses-0%s/anat/*MNI152NLin2009cAsym*',bidsDir,sub{iSub},num2str(ses)));
        for ii = 1:numel(mni)
            delete([mni(ii).folder '/' mni(ii).name])
        end
        
    end
    
end
toc

%% rename

dir_path = fullfile('/Volumes/Vision/MRI/DecodingPublic/'); 
mydir = dir(fullfile(dir_path, '**/ses-01/anat/*_defaced.nii.gz'));

for i = 1:numel(mydir);
   folders = strsplit(mydir(i).folder, filesep);  %split path into individual folders
   newname = [folders{end-1}, '.log'];
   f1 = fullfile(mydir(i).folder, mydir(i).name);
   f2 = [f1(1:end-15) '.nii.gz'];
   movefile(f1,f2); %rename file  
end

%% fix intended for 
bidsDir = '/Volumes/Vision/MRI/DecodingPublic';
subject = {'0201','0202','0204','0205','0206','0228','0229','0248','0903'};
ses = {'01','02','03','04'};
for iSub = 1:9

    for ii = 1:length(ses)
        
        % Fix fmap json files by adding run information
        jsons = dir(fullfile(bidsDir, ['sub-' subject{iSub}], ['ses-' ses{ii}], 'fmap', '*.json'));
        for ff = 1:length(jsons) % for each json file
            fname = fullfile(jsons(ff).folder, jsons(ff).name);
            str = fileread(fname);
            val = jsondecode(str);
            
            funcfiles = dir(fullfile(bidsDir, ['sub-' subject{iSub}], ['ses-' ses{ii}], 'func', '*.nii.gz'));
            funcfiles = struct2cell(funcfiles);
            funcfiles = strcat(['ses-' ses{ii}],'/func/', funcfiles(1,:)');
            val.IntendedFor = funcfiles;
            val.TaskName = "3dmotion";
            str = jsonencode(val);
            % Make the json output fi       le more human readable
            str = strrep(str, ',"', sprintf(',\n"'));
            str = strrep(str, '[{', sprintf('[\n{\n'));
            str = strrep(str, '}]', sprintf('\n}\n]'));
            
            fid = fopen(fname,'w');
            fwrite(fid,str);
            fclose(fid);
            
        end
        
    end
end
