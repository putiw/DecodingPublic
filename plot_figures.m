clear all;clc;close all; restoredefaultpath;
addpath(genpath([pwd '/helper_functions']));
bidsDir = '/Volumes/Vision/MRI/DecodingPublic';
fsDir = '/Applications/freesurfer/7.2.0';
gitDir = '~/Documents/GitHub';
set_up(bidsDir,gitDir,fsDir)
%% load takfap results

% 3 datasets (200,all,chance), 9 subjects, 2 conditions, 22 rois,
sub = {'0201','0202','0204','0205','0206','0228','0229','0248','0903'};
nRoi = 22;
con = {'ses-0102','ses-0304'};
datasets = {'TAFKAP','TAFKAP_200voxels','TAFKAP_200voxels_chance'};

CNFM =  cell(nRoi,length(sub),2,3);
acc = zeros(nRoi,length(sub),2,3);
se = zeros(nRoi,length(sub),2,3);
for iSub = 1:length(sub)
    for iCon = 1:length(con)
        for iSet = 1:length(datasets)
            f1 = load(sprintf('%s/derivatives/resultMat/sub-%s_%s_%s.mat',bidsDir,sub{iSub},con{iCon},datasets{iSet}));
            for iRoi = 1:nRoi
                CNFM{iRoi,iSub,iCon,iSet} = f1.saveresult{iRoi,1}(1:8,1:8)./100;
                tmp = squeeze(mean(reshape((f1.pres{iRoi}==f1.ests{iRoi}),16,32,[])));
                acc(iRoi,iSub,iCon,iSet) = 100.*mean(tmp);
                se(iRoi,iSub,iCon,iSet) = 100.*std(tmp)/sqrt(32);
            end
        end
    end
end

meanacc = reshape(mean(acc,2),nRoi,length(con),length(datasets));
see = zeros(nRoi,length(con),length(datasets));
for iRoi = 1:nRoi
    for iCon = 1:length(con)
        for iSet = 1:length(datasets)
            see(iRoi,iCon,iSet) = std(acc(iRoi,:,iCon,iSet))/sqrt(length(sub));
        end
    end
end

miss = zeros(nRoi,length(sub),2);
whichSet = 2;
for iRoi = 1:nRoi
    for iSub = 1:length(sub)
        miss(iRoi,iSub,1) = (CNFM{iRoi,iSub,1,whichSet}(1,2) + CNFM{iRoi,iSub,1,whichSet}(5,4) + CNFM{iRoi,iSub,1,whichSet}(5,6) + CNFM{iRoi,iSub,1,whichSet}(1,8))/4;
        miss(iRoi,iSub,2) = (CNFM{iRoi,iSub,1,whichSet}(3,2) + CNFM{iRoi,iSub,1,whichSet}(3,4) + CNFM{iRoi,iSub,1,whichSet}(7,6) + CNFM{iRoi,iSub,1,whichSet}(7,8))/4;
    end
end

cnfm = cell(nRoi,length(con),length(datasets));
for iRoi = 1:nRoi
    for iCon = 1:length(con)
        for iSet = 1:length(datasets)
            sumMat = zeros(8,8);
            for iSub = 1:length(sub)
                sumMat = sumMat + CNFM{iRoi,iSub,iCon,iSet};
            end
            cnfm{iRoi,iCon,iSet} = sumMat./length(sub);
        end
    end
end

%% figure 3 (Confusion matrix) - (A C)
close all
figure;

whichSet = [1 2];
whichRoi = [1 9];

ind = 1;
for iSet = 1:numel(whichSet)
    for iRoi = 1:length(whichRoi);
        subplot(2,2,ind);ind=ind+1;
        title([f1.roi{whichRoi(iRoi)}])
        vals = cnfm{whichRoi(iRoi),1,iSet}.*100;
        vals(9,:) = vals(1,:);
        vals(:,9) = vals(:,1);
        imagesc(vals);
        cb = colorbar(); caxis([0 50]);
        cb.Ruler.TickLabelFormat = '%d%%';
        colorgroup = [200 230 255; 255 255 255; 255 122 122]./255;
        mymap1 = linspace(0, 1, 10)'*colorgroup(2,:)+(1 - linspace(0, 1, 10))'*colorgroup(1,:);
        mymap2 = linspace(0, 1, round(10.*((50-12.5)/12.5)))'*colorgroup(3,:)+(1 - linspace(0, 1, round(10.*((50-12.5)/12.5))))'*colorgroup(2,:);
        colormap([mymap1;mymap2]);
        set(gca,'YDir','normal');
        ylabel('Decoded Direction');
        xlabel('Presented Direction');
        axis square;
        title([f1.roi{whichRoi(iRoi)}])
    end
end

%% figure 3 (bar V1 vs hMT & all vs 200 voxels) - (D)

close all
whichRoi = [1 9];
bar1X = [1.2 3.8 2.2 4.8];
bar1Y = [meanacc(whichRoi,1,1);meanacc(whichRoi,1,2)]';
error1Y = [see(whichRoi,1,1);see(whichRoi,1,1)]';

f3 = figure('Renderer', 'painters', 'Position', [10 10 350 298]);
hold on
bar1 = bar(bar1X,bar1Y,'linewidth',1);
errorbar(bar1X,bar1Y,error1Y,'k.','linewidth',1);
bar1.FaceColor = [1 1 1];
plot([0 6],[12.5 12.5],'k--','linewidth',1);
f3.CurrentAxes.XTick = [1.7 4.3];
f3.CurrentAxes.YTick = [15 30 45];
f3.CurrentAxes.XTickLabel = {'V1','hMT'};
ylabel('Percentage Correct (%)')
title({'Decoding Accuracy'})
l = legend('200 voxels','all voxels');
xlabel('ROI')
ylim([0 50])
xlim([0 6])
ax = gca;
ax.LineWidth = 1;
set(gca,'FontSize',15)
%% figure 4 Miss-classification
close all;clc;
whichRoi = [1 9];
bar1X = [1.2 3.8 2.2 4.8];
bar1Y = [mean(miss(whichRoi,:,2),2);mean(miss(whichRoi,:,1),2)]'.*100;
error1Y = [std(miss(whichRoi,:,2),1,2)/sqrt(length(sub)-1);std(miss(whichRoi,:,1),1,2)/sqrt(length(sub)-1)]'.*100;
f4 = figure('Renderer', 'painters', 'Position', [10 10 350 298]);
hold on
bar(bar1X,bar1Y,'linewidth',1);
errorbar(bar1X,bar1Y,error1Y,'k.','linewidth',1);
plot([0 6],[12.5 12.5],'k--','linewidth',1);
f4.CurrentAxes.XTick = [1.7 4.3];
f4.CurrentAxes.YTick = [15 30 45];
f4.CurrentAxes.XTickLabel = {'V1','hMT'};
ylabel('Decoding accuracy (%)')
title({'Decoding Misclassification'});
xlabel('ROI')

%% figure 5 All ROI bar H/V
close all
f5 = figure('Renderer', 'painters', 'Position', [10 10 950 310]);
hold on
bar1 = bar(1:nRoi,meanacc2(:,1),0.4,'linewidth',1);
bar2 = bar((1:nRoi)+0.4,meanacc2(:,2),0.4,'linewidth',1);
errorbar(1:nRoi,meanacc2(:,1),see2(:,1),'linewidth',1,'Capsize',0);
errorbar((1:nRoi)+0.4,meanacc2(:,2),see2(:,2),'linewidth',1,'Capsize',0);
plot([0 23],[12.5 12.5],'k--','linewidth',1);
legend('Horizontal','Vertical')
f5.CurrentAxes.XTick = (1:nRoi)+0.25;
f5.CurrentAxes.YTick = [15 30 45];
f5.CurrentAxes.XTickLabel = Hor.roi;
ylim([0 50])
h=gca;
ylabel('Decoding Accuracy (%)')
title({'Decoding accuracy averaged across nine subjects'})
xlabel('ROI')

%% Searchlight

% load fsaverage6 .inflated surf
surfL = fullfile([bidsDir '/derivatives/freesurfer/fsaverage6/surf/lh.inflated']);
surfL = fs_read_surf(surfL);
surfR = fullfile([bidsDir '/derivatives/freesurfer/fsaverage6/surf/rh.inflated']);
surfR = fs_read_surf(surfR);
num = size(surfL.coord,2);
% find what is inside for each searchlight
idx = zeros(num,100,2);
for iVer = 1:num
    distance = sqrt(sum((surfL.coord - surfL.coord(:,iVer)).^2));
    [~,idx(iVer,:,1)] = mink(distance,100);
    distance = sqrt(sum((surfR.coord - surfR.coord(:,iVer)).^2));
    [~,idx(iVer,:,2)] = mink(distance,100);
end
surfNeighbour = [idx(:,:,1)' idx(:,:,2)'];

data = zeros(numel(sub),num*2,numel(con));

for iSub = 1:length(sub)
    for iCon = 1:length(con)
        load(sprintf('%s/derivatives/resultMat/sub-%s_%s_searchlight_fsaverage6.mat',bidsDir,sub{iSub},con{iCon}));
        data(iSub,:,iCon) = acc;
    end
end

meanAcc = nanmean(data,1);
cvnlookup('fsaverage6',[],meanAcc(:,:,1)'); % horizonal acc
cvnlookup('fsaverage6',[],meanAcc(:,:,2)'); % vertical acc
meanAcc(:,:,2)'

load('helper_functions/searchlightChance.mat');

pval = zeros(num*2,numel(con));
for iCon = 1:numel(con)
    for kk = 1:num*2
        pval(kk,iCon) = sum(pcChance(:,kk)>meanAcc(:,kk,2)')/1000;
    end
end
cvnlookup('fsaverage6',[],meanAcc(:,:,1)'); % horizonal p val
cvnlookup('fsaverage6',[],meanAcc(:,:,2)'); % vertical p val

[~,p] = ttest(data(:,:,1),data(:,:,2));
neighbour = surfNeighbour(2:7,:);
difference = data(:,:,1)-data(:,:,2);
alpha = 0.04;
idx = 1;
diffCluster = cell(1,805);
for iV = 1:numel(difference)
    if ismember(iV,cat(1,diffCluster{:})) == 0 && difference(iV) >= alpha
        clusterTemp = find_neighbour(neighbour,difference,iV,alpha);
        diffCluster{1,idx} = [iV ; clusterTemp];
        idx = idx + 1;
    end
end
clusterSize =cellfun(@numel,diffCluster)';
trackSize = nan(surfSize*2,1);
for iC = 1:size(diffCluster,2)
    trackSize(diffCluster{1,iC},1) = numel(diffCluster{1,iC});
end

maxSize = 1;
pp =  nan(size(p));
pp(trackSize>maxSize) = p(trackSize>maxSize);
cvnlookup('fsaverage6',[],pp); % h-v
