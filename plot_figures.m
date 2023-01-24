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

clear all;clc;close all;
addpath('helper_functions')
projectDir = '/Users/pw1246/Desktop/MRI/Decoding';
setup_user_decoding('Puti',projectDir);

aa = [];
file1 = '/Users/pw1246/Desktop/dd/l.fsaverage6_horizontal.mgh';
file2 = '/Users/pw1246/Desktop/dd/r.fsaverage6_horizontal.mgh';

funcL = MRIread(file1);
aa{1} = squeeze(funcL.vol)';
funcR = MRIread(file2);
aa{2} = squeeze(funcR.vol)';

file1 = '/Users/pw1246/Desktop/dd/l.fsaverage6_vertical.mgh';
file2 = '/Users/pw1246/Desktop/dd/r.fsaverage6_vertical.mgh';


funcL = MRIread(file1);
aa{3} = squeeze(funcL.vol)';
funcR = MRIread(file2);
aa{4} = squeeze(funcR.vol)';

file1 = '/Users/pw1246/Desktop/dd/l.H_V.mgh';
file2 = '/Users/pw1246/Desktop/dd/r.H_V.mgh';

funcL = MRIread(file1);
aa{5} = squeeze(funcL.vol)';
funcR = MRIread(file2);
aa{6} = squeeze(funcR.vol)';
                    
file1 = '/Users/pw1246/Dropbox (RVL)/pw1246/Home/Decoding/motion/maps/l.sub-0201_0102.mgh';
file2 = '/Users/pw1246/Dropbox (RVL)/pw1246/Home/Decoding/motion/maps/r.sub-0201_0102.mgh';

funcL = MRIread(file1);
aa{7} = squeeze(funcL.vol)';
funcR = MRIread(file2);
aa{8} = squeeze(funcR.vol)';

                    
file1 = '/Users/pw1246/Dropbox (RVL)/pw1246/Home/Decoding/motion/maps/l.sub-0201_0304.mgh';
file2 = '/Users/pw1246/Dropbox (RVL)/pw1246/Home/Decoding/motion/maps/r.sub-0201_0304.mgh';

funcL = MRIread(file1);
aa{9} = squeeze(funcL.vol)';
funcR = MRIread(file2);
aa{10} = squeeze(funcR.vol)';
%%
close all
gitDir = '~/Documents/Github';
fsDir = '/Applications/freesurfer/7.2.0';
bidsDir = projectDir;
setenv('FREESURFER_HOME',fsDir)
setenv('SUBJECTS_DIR',[bidsDir '/derivatives/freesurfer'])

lVal = cvntransfertosubject('fsaverage6','fsaverage', aa{1}, 'lh', 'nearest', 'orig', 'orig');
rVal = cvntransfertosubject('fsaverage6','fsaverage', aa{2}, 'rh', 'nearest', 'orig', 'orig');

%%


vals = [lVal;rVal];
valsalpha = ~isnan(vals);
vals(isnan(vals)) = 0;

vals(vals<=0.15) = 0;

bins = 0.125:0.01:0.22;
cmap0 = cmaplookup(bins,min(bins),max(bins),[],(cmapsign4));

%cvnlookup('fsaverage',13,vals,[min(bins) max(bins)],hot,[],[],[],{'roiname',{'Kastner2015'},'roicolor',{'w'},'drawrpoinames',0,'roiwidth',{3},'fontsize',20,'overlayalpha',valsalpha});

%[rawimg,Lookup,rgbimg] = cvnlookup('fsaverage',13,vals,[min(bins) max(bins)],cmap0,min(bins),[],0,{'roiname',{'Glasser2016'},'roicolor',{'g'},'drawroinames',0,'roiwidth',{0.5},'fontsize',20}); %Kastner2015
[rawimg,Lookup,rgbimg] = cvnlookup('fsaverage',13,vals,[min(bins) max(bins)],cmap0,min(bins),[],0); %Kastner2015

color = [0.5];
[r,c,t] = size(rgbimg);
[i j] = find(all(rgbimg == repmat(color,r,c,3),3));

for ii = 1: length(i)
    rgbimg(i(ii),j(ii),:) = ones(1,3);
end

a = imagesc(rgbimg); axis image tight;

    
set(gcf,'Position',[277 119 1141 898])
axis off
hold on

plot(0,0);
colormap(cmap0);
hcb=colorbar('SouthOutside');
hcb.Ticks = [0 1];
hcb.TickLabels = {'12.5%';'25%'}
hcb.FontSize = 25
hcb.Label.String = 'Searchlight Decoding Accuracy (horizontal motion)'
hcb.TickLength = 0.001;


%print('fsaverage_V','-dpdf','-bestfit','-r0')
%% roi for reviewer #2


roi = {'V1d','TO1','V1v'};
roimask = cell(1,numel(roi));
for i = 1:numel(roi)
[tmpl, ~, ~] = cvnroimask('fsaverage','lh','Kastner2015',roi{i},[],[]);
[tmpr,~,~] = cvnroimask('fsaverage','rh','Kastner2015',roi{i},[],[]);
roimask{i} = [tmpl{:};tmpr{:}];
end
roimas = {};
roimas{1} = roimask{1}|roimask{3};
roimas{2} = roimask{2};



%%

lVal = cvntransfertosubject('fsaverage6','fsaverage', aa{1}, 'lh', 'nearest', 'orig', 'orig');
rVal = cvntransfertosubject('fsaverage6','fsaverage', aa{2}, 'rh', 'nearest', 'orig', 'orig');
valsH = [lVal;rVal];
lVal = cvntransfertosubject('fsaverage6','fsaverage', aa{3}, 'lh', 'nearest', 'orig', 'orig');
rVal = cvntransfertosubject('fsaverage6','fsaverage', aa{4}, 'rh', 'nearest', 'orig', 'orig');
valsV = [lVal;rVal];


lVal = cvntransfertosubject('fsaverage6','fsaverage', data(5,1:40962,1)', 'lh', 'nearest', 'orig', 'orig');
rVal = cvntransfertosubject('fsaverage6','fsaverage', data(5,40963:end,1)', 'rh', 'nearest', 'orig', 'orig');
valsH5 = [lVal;rVal];
lVal = cvntransfertosubject('fsaverage6','fsaverage', data(5,1:40962,2)', 'lh', 'nearest', 'orig', 'orig');
rVal = cvntransfertosubject('fsaverage6','fsaverage', data(5,1:40962,2)', 'rh', 'nearest', 'orig', 'orig');
valsV5 = [lVal;rVal];
lVal = cvntransfertosubject('fsaverage6','fsaverage', data(1,1:40962,1)', 'lh', 'nearest', 'orig', 'orig');
rVal = cvntransfertosubject('fsaverage6','fsaverage', data(1,40963:end,1)', 'rh', 'nearest', 'orig', 'orig');
valsH1 = [lVal;rVal];
lVal = cvntransfertosubject('fsaverage6','fsaverage', data(1,1:40962,2)', 'lh', 'nearest', 'orig', 'orig');
rVal = cvntransfertosubject('fsaverage6','fsaverage', data(1,1:40962,2)', 'rh', 'nearest', 'orig', 'orig');
valsV1 = [lVal;rVal];

%%
close all; 
figure;
hold on
for i = 1:2
mV = mean(valsV(roimas{i}))*100;
mH = mean(valsH(roimas{i}))*100;
sV = std(valsV(roimas{i}))./sqrt(sum(roimas{i}))*100;
sH = std(valsH(roimas{i}))./sqrt(sum(roimas{i}))*100;
bar([i i+0.3],[mH mV])
end
%%
vals = zeros(327684,9,2);
for i = 1:9
    lVal = cvntransfertosubject('fsaverage6','fsaverage', data(i,1:40962,1)', 'lh', 'nearest', 'orig', 'orig');
rVal = cvntransfertosubject('fsaverage6','fsaverage', data(i,40963:end,1)', 'rh', 'nearest', 'orig', 'orig');
vals(:,i,1) = [lVal;rVal];
lVal = cvntransfertosubject('fsaverage6','fsaverage', data(i,1:40962,2)', 'lh', 'nearest', 'orig', 'orig');
rVal = cvntransfertosubject('fsaverage6','fsaverage', data(i,1:40962,2)', 'rh', 'nearest', 'orig', 'orig');
vals(:,i,2) = [lVal;rVal];
end
%%
close all; 
figure;
subplot(2,1,1)
hold on
histogram(reshape(vals(roimas{1},:,1),numel(vals(roimas{1},:,1)),1),'FaceColor',[1 0 0]);
histogram(reshape(vals(roimas{1},:,2),numel(vals(roimas{1},:,2)),1),'FaceColor',[0 0 1]);
subplot(2,1,2)
hold on
histogram(reshape(vals(roimas{2},:,1),numel(vals(roimas{2},:,1)),1),'FaceColor',[1 0 0]);
histogram(reshape(vals(roimas{2},:,2),numel(vals(roimas{2},:,2)),1),'FaceColor',[0 0 1]);
%%
lVal = cvntransfertosubject('fsaverage6','fsaverage', p(1:40962)', 'lh', 'nearest', 'orig', 'orig');
rVal = cvntransfertosubject('fsaverage6','fsaverage', p(40963:end)', 'rh', 'nearest', 'orig', 'orig');
ppphv = [lVal;rVal];
%%
total = [aa{3}' aa{4}'];
myP = zeros(81924,1);
for kk = 1:81924   
    myP(kk) = sum(pcChance(:,kk)>total(1,kk))/1000;   
end
ppv=myP;
lVal = cvntransfertosubject('fsaverage6','fsaverage', myP(1:40962,1), 'lh', 'nearest', 'orig', 'orig');
rVal = cvntransfertosubject('fsaverage6','fsaverage', myP(40963:end,1), 'rh', 'nearest', 'orig', 'orig');
ppv = [lVal;rVal];
%%
total = [aa{1}' aa{2}'];
myP = zeros(81924,1);
for kk = 1:81924   
    myP(kk) = sum(pcChance(:,kk)>total(1,kk))/1000;   
end
pph = myP;
lVal = cvntransfertosubject('fsaverage6','fsaverage', myP(1:40962,1), 'lh', 'nearest', 'orig', 'orig');
rVal = cvntransfertosubject('fsaverage6','fsaverage', myP(40963:end,1), 'rh', 'nearest', 'orig', 'orig');
pph = [lVal;rVal];
%%
v1l = cvntransfertosubject('fsaverage','fsaverage6', double(roimas{1}(1:163842,1)), 'lh', 'nearest', 'orig', 'orig');
v1r = cvntransfertosubject('fsaverage','fsaverage6', double(roimas{1}(163843:end,1)), 'rh', 'nearest', 'orig', 'orig');
roimas{1} = [v1l;v1r];
%%
v1ph_v=p(roimas{1});
v1pv=ppv(roimas{1});
v1ph=pph(roimas{1});
dV1h = data(:,roimas{1},1);
dV1v = data(:,roimas{1},2);
forreview = find(v1ph<0.05&v1pv>0.05&v1ph_v>0.05);
%%
i = forreview(randperm(385,1)); %356 511 852 44 251
[v1ph(i) v1pv(i) v1ph_v(i)]
%%
i = 210; %496
ddhv = [dV1h(:,i) dV1v(:,i) dV1h(:,i)-dV1v(:,i)].*100
mv = total(roimas{1},1);
[mean(dV1h(:,i)) mean(dV1v(:,i))]
[~, x, ~,~] = ttest(dV1h(:,i),dV1v(:,i));
[v1ph(i) v1pv(i) v1ph_v(i) x]
%%
i = 210
[v1ph(i) v1pv(i) v1ph_v(i) x]
close all
figure; hold on;
h1 = histogram(chancev1(:,210),20);h1.FaceAlpha = 0.2;
plot([dV1v(:,i) dV1v(:,i)]',repmat([0 150],9,1)','b--','linewidth',1);
plot([dV1h(:,i) dV1h(:,i)]',repmat([0 150],9,1)','r--','linewidth',1);
plot([mean(dV1v(:,i)) mean(dV1v(:,i))],[0 150],'b-','linewidth',3);
plot([mean(dV1h(:,i)) mean(dV1h(:,i))],[0 150],'r-','linewidth',3);
set(gca,'fontsize',15)
%%
val11 = [];
for i = 1:9
lVal = cvntransfertosubject('fsaverage6','fsaverage', data(i,1:40962,1)', 'lh', 'nearest', 'orig', 'orig');
rVal = cvntransfertosubject('fsaverage6','fsaverage', data(i,40963:end,1)', 'rh', 'nearest', 'orig', 'orig');
val11(:,i,1) = [lVal;rVal];
lVal = cvntransfertosubject('fsaverage6','fsaverage', data(i,1:40962,2)', 'lh', 'nearest', 'orig', 'orig');
rVal = cvntransfertosubject('fsaverage6','fsaverage', data(i,40963:end,2)', 'rh', 'nearest', 'orig', 'orig');
val11(:,i,2) = [lVal;rVal];
end
%%
i = 1
[v1ph_c(i) v1pv_c(i) v1ph_v(i)]

%%
subplot(2,1,1)
hold on
histogram(valsH1(roimas{1}),'FaceColor',[1 0 0])
histogram(valsV1(roimas{1}),'FaceColor',[0 0 1])
xlim([0.125 0.25])
subplot(2,1,2)
hold on
histogram(valsH1(roimas{2}),'FaceColor',[1 0 0])
histogram(valsV1(roimas{2}),'FaceColor',[0 0 1])
xlim([0.125 0.25])
drawnow


%% p 
load('searchlightChance');

total = [aa{3}' aa{4}'];
myP = zeros(81924,1);
for kk = 1:81924
   
    myP(kk) = sum(pcChance(:,kk)>total(1,kk))/1000;
    
end

lVal = cvntransfertosubject('fsaverage6','fsaverage', myP(1:40962,1), 'lh', 'nearest', 'orig', 'orig');
rVal = cvntransfertosubject('fsaverage6','fsaverage', myP(40963:end,1), 'rh', 'nearest', 'orig', 'orig');

close all

vals = [lVal;rVal];

valsalpha = ~isnan(vals)&vals<=0.1;
vals(isnan(vals)) = 0;

bins = [logspace(-3,-2,30) logspace(-2,-1.3010,30) logspace(-1.3010,-1,30)];%0.001:0.001:0.05; 
cmap0 = cmaplookup(bins,min(bins),max(bins),[],(cmapsign4));

%cvnlookup('fsaverage',13,vals,[min(bins) max(bins)],flipud(jet),[],[],[],{'roiname',{'Kastner2015'},'roicolor',{'w'},'drawrpoinames',0,'roiwidth',{3},'fontsize',20,'overlayalpha',valsalpha});


[rawimg,Lookup,rgbimg] = cvnlookup('fsaverage',13,vals,[min(bins) max(bins)],flipud(jet),[],[],0,{'roiname',{'Kastner2015'},'roicolor',{'w'},'drawroinames',1,'roiwidth',{0.5},'fontsize',20,'overlayalpha',valsalpha});
color = [0.5];
[r,c,t] = size(rgbimg);
[i j] = find(all(rgbimg == repmat(color,r,c,3),3));

for ii = 1: length(i)
    rgbimg(i(ii),j(ii),:) = ones(1,3);
end

a = imagesc(rgbimg); axis image tight;

    
set(gcf,'Position',[277 119 1141 898])
axis off
hold on

plot(0,0);
colormap(flipud(jet));
hcb=colorbar('SouthOutside');
hcb.Ticks = [0 0.33 0.66 1];
hcb.TickLabels = {'0.001';'0.01';'0.05';'0.1'}; 

%hcb.TickLabels = {sprintf('%.3f',min(bins));sprintf('%.3f',bins(end/2));sprintf('%.3f',max(bins))}; 
hcb.FontSize = 25
hcb.Label.String = 'Searchlight Decoding p-val (horizontal motion)'
hcb.TickLength = 0.001;

%% h - v diff
path =  '/Users/pw1246/Dropbox (RVL)/pw1246/Home/Decoding/motion';
 load('surfNeighbour.mat');
  %load('/Users/pw1246/Desktop/dd/final/trackSize.mat');
sub = {'sub-0201','sub-0202','sub-0204','sub-0205','sub-0206','sub-0228','sub-0229','sub-0248','sub-0903'};                                          % subject ID
ses = {'0102','0304'};
surfSize = 40962;
data = zeros(numel(sub),surfSize*2,numel(ses));
for iSub = 1:numel(sub)
    for iSes = 1:numel(ses)
            lh = gifti([path,'/maps/l.',sub{iSub},'_',ses{iSes},'.mgz']);
            rh = gifti([path,'/maps/r.',sub{iSub},'_',ses{iSes},'.mgz']);

        data(iSub,:,iSes) = cat(1,[lh.cdata',rh.cdata']);
    end
end
[~,p] = ttest(data(:,:,1),data(:,:,2));
 neighbour = surfNeighbour(2:7,:);
difference = [aa{5}' aa{6}'];
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

gitDir = '~/Documents/Github';
fsDir = '/Applications/freesurfer/7.2.0';
bidsDir = projectDir;
setenv('FREESURFER_HOME',fsDir)
setenv('SUBJECTS_DIR',[bidsDir '/derivatives/freesurfer'])

%%
maxSize = 1;
%p(p>0.1)=NaN;
pp =  nan(size(p));
pp(trackSize>maxSize) = p(trackSize>maxSize);
pp = p;
lVal = cvntransfertosubject('fsaverage6','fsaverage', pp(1,1:40962)', 'lh', 'nearest', 'orig', 'orig');
rVal = cvntransfertosubject('fsaverage6','fsaverage', pp(1,40963:end)', 'rh', 'nearest', 'orig', 'orig');
%
close all

vals = [lVal;rVal];

valsalpha = ~isnan(vals);
vals(isnan(vals)) = 0;
bins = [logspace(-3,-2,30) logspace(-2,-1.3010,30) logspace(-1.3010,-1,30)];%0.001:0.001:0.05; 
cmap0 = cmaplookup(bins,min(bins),max(bins),[],(cmapsign4));
[rawimg,Lookup,rgbimg] = cvnlookup('fsaverage',13,vals,[min(bins) max(bins)],flipud(jet),[],[],0,{'roiname',{'Kastner2015'},'roicolor',{'w'},'drawroinames',1,'roiwidth',{0.5},'fontsize',20,'overlayalpha',valsalpha});

color = [0.5];
[r,c,t] = size(rgbimg);
[i j] = find(all(rgbimg == repmat(color,r,c,3),3));

for ii = 1: length(i)
    rgbimg(i(ii),j(ii),:) = ones(1,3);
end

a = imagesc(rgbimg); axis image tight;

    
set(gcf,'Position',[277 119 1141 898])
axis off
hold on

plot(0,0);
colormap(flipud(jet));
hcb=colorbar('SouthOutside');
 hcb.Ticks = [0 0.33 0.66 1];
hcb.TickLabels = {'0.001','0.01','0.05','0.1'}; 
hcb.FontSize = 25
hcb.Label.String = 'Difference between horizontal and vertical motion'
hcb.TickLength = 0.001;


%% for bas

close all
gitDir = '~/Documents/Github';
fsDir = '/Applications/freesurfer/7.2.0';
bidsDir = projectDir;
setenv('FREESURFER_HOME',fsDir)
setenv('SUBJECTS_DIR',[bidsDir '/derivatives/freesurfer'])

lVal = cvntransfertosubject('fsaverage6','fsaverage', aa{9}, 'lh', 'nearest', 'orig', 'orig');
rVal = cvntransfertosubject('fsaverage6','fsaverage', aa{10}, 'rh', 'nearest', 'orig', 'orig');
%%
figure(2)
vals = [lVal;rVal];
valsalpha = ~isnan(vals);
vals(isnan(vals)) = 0;

vals(vals<=0.15) = 0;

bins = 0.125:0.01:0.3;
cmap0 = cmaplookup(bins,min(bins),max(bins),[],(cmapsign4));

%cvnlookup('fsaverage',13,vals,[min(bins) max(bins)],hot,[],[],[],{'roiname',{'Kastner2015'},'roicolor',{'w'},'drawrpoinames',0,'roiwidth',{3},'fontsize',20,'overlayalpha',valsalpha});
%'HCP-MMP1'
[rawimg,Lookup,rgbimg] = cvnlookup('fsaverage',13,vals,[min(bins) max(bins)],cmap0,min(bins),[],0,{'roiname',{'Glasser2016'},'roicolor',{'w'},'drawroinames',1,'roiwidth',{1},'fontsize',20}); %MT_exvivo %Kastner2015

color = [0.5];
[r,c,t] = size(rgbimg);
[i j] = find(all(rgbimg == repmat(color,r,c,3),3));

for ii = 1: length(i)
    rgbimg(i(ii),j(ii),:) = ones(1,3);
end

a = imagesc(rgbimg); axis image tight;

    
set(gcf,'Position',[277 119 1141 898])
axis off
hold on

plot(0,0);
colormap(cmap0);
hcb=colorbar('SouthOutside');
hcb.Ticks = [0 1];
hcb.TickLabels = {'12.5%';'30%'}
hcb.FontSize = 25
hcb.Label.String = 'Searchlight Decoding Accuracy (vertical motion)'
hcb.TickLength = 0.001;


%%