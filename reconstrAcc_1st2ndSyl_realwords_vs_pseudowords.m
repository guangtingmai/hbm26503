function reconstrAcc_1st2ndSyl_realwords_vs_pseudowords

pathfolder = 'C:\Users\mszgm1\OneDrive - The University of Nottingham\LEL 2013';
addpath(genpath(pathfolder));
folder = sprintf('%s//analyses 2022',pathfolder);
decimate_factor = 4;

subjectID1 = '03_20am';  subjectID2 = '03_20eve'; subjectID3 = '03_22pm';   subjectID4 = '03_22eve'; subjectID5 = '03_24am';  
subjectID6 = '03_24pm';  subjectID7 = '03_24eve'; subjectID8 = '03_26pm';   subjectID9 = '03_29am';  subjectID10 = '03_29pm'; 
subjectID11 = '03_30am'; subjectID12 = '03_30pm'; subjectID13 = '03_31am';  subjectID14 = '03_31pm'; subjectID15 = '04_02eve'; 
subjectID16 = '04_06am'; subjectID17 = '04_06pm'; subjectID18 = '04_06eve'; subjectID19 = '04_07am'; subjectID20 = '04_07pm';
subjectID = {subjectID1, subjectID2, subjectID3, subjectID4, subjectID5, subjectID6, subjectID7, subjectID8, subjectID9, subjectID10, ...
     subjectID11, subjectID12, subjectID13, subjectID14, subjectID15, subjectID16, subjectID17, subjectID18, subjectID19, subjectID20};

stimfilename = sprintf('%s//dataStim_RealPseudo_norm.mat',folder);
a = load(stimfilename);
stim = a.stim;
stim1stSyl = stim.data(10,:);
stim2ndSyl = stim.data(11,:);

fs = 64;
corr_dur = round(0.1*fs); % 200 ms correlation window

% {'Fp1', 'AF3', 'F7', 'F3', 'FC1', 'FC5', 'T7', 'C3', ...
%    1      2     3     4      5      6     7     8  
%  'CP1', 'CP5', 'P7', 'P3', 'Pz', 'PO3', 'O1', 'Oz', ...
%    9      10    11    12    13    14     15    16
%  'O2', 'PO4', 'P4', 'P8', 'CP6', 'CP2', 'C4', 'T8', ...
%   17     18    19    20     21    22     23    24
%  'FC6', 'FC2', 'F4', 'F8', 'AF4', 'Fp2', 'Fz', 'Cz'};
%   25     26     27    28    29     30     31    32

electrodes = [1:10,21:32];

freqRange = 2;
switch freqRange
    case 1 % Delta
        FreqRange = 'Delta';
    case 2 % Theta
        FreqRange = 'Theta';
    otherwise % DeltaTheta
        FreqRange = 'DeltaTheta';
end
features = 'SpectrogramPhoneticsPhonemes';
trfname = sprintf('mTRF_encoding_%s_RealPseudoBack_%s',FreqRange,features);
partial_or_not = 0;
trfname2 = sprintf('%s_0-500ms_eegraw_64HzFs',trfname);
trfname = sprintf('%s_subperiods',trfname);
bb = load(sprintf('%s//%s.mat',folder,trfname2));

for s = 1:20
    aa = load(sprintf('%s//%s_pred_%s.mat',folder,trfname,subjectID{s}));
    lag_pts = length(aa.pred);
    ddd = [];
    if partial_or_not, ddd = bb.eeg_partial{s}; end
    [r_real_1stSyl(:,:,s), r_real_2ndSyl(:,:,s), r_pseudo_1stSyl(:,:,s), r_pseudo_2ndSyl(:,:,s)] = ...
        reconstrAcc_1st2ndSyl_realwords_vs_pseudowords_singleSubject...
        (subjectID{s},aa.pred,ddd,freqRange,stim1stSyl,stim2ndSyl,partial_or_not);
    s
end

p1 = [];
for ttt = 1:size(r_real_1stSyl,1)
    p1(ttt) = bootstrapping(squeeze(mean(r_real_1stSyl(ttt,electrodes,:),2)-mean(r_pseudo_1stSyl(ttt,electrodes,:),2)),10000);
end

p2 = [];
for ttt = 1:size(r_real_2ndSyl,1)
    p2(ttt) = bootstrapping(squeeze(mean(r_real_2ndSyl(ttt,electrodes,:),2)-mean(r_pseudo_2ndSyl(ttt,electrodes,:),2)),10000);
end
[~,~,~,adj_p]=fdr_bh([p1,p2]);
adj_p1 = adj_p(1:length(adj_p)/2);
adj_p2 = adj_p(length(adj_p)/2+1:end);

% r_real_2ndSyl 
figure(1); clf;
r4plot = [];
r4plot(:,:,1) = squeeze(mean(r_real_1stSyl(:,electrodes,:),2)).';
r4plot(:,:,2) = squeeze(mean(r_pseudo_1stSyl(:,electrodes,:),2)).';
r4plot(:,:,3) = squeeze(mean(r_real_2ndSyl(:,electrodes,:),2)).';
r4plot(:,:,4) = squeeze(mean(r_pseudo_2ndSyl(:,electrodes,:),2)).';
interpFactor = 2;
LineColor = {[1 0 0 0.2],[0 0 1 0.2],'r','b'};
xyRange = {[0,520],[0.01,0.14]};
xyLabel = {'Time lag (ms)','Reconstruction accuracy'};
% FaceColor = {'r','b','r','b'};
% FaceAlpha = [0.00,0.00,0.25,0.25];
varargin1 = {};
p = [adj_p1;adj_p2];
sigLineVal = [0.025,0.02];
sigColor = {'#6E2C00','#6E2C00'};
varargin2 = {p,sigLineVal,sigColor};

plotLineGraph_with_seShape(r4plot,fs,interpFactor,LineColor,xyRange,xyLabel,...
    varargin1,varargin2);

r_real_vs_pseudo_2ndSyl = r_real_2ndSyl - r_pseudo_2ndSyl;
period = 1:round(0.2*fs);
r_real_vs_pseudo_2ndSyl = (squeeze(mean(r_real_vs_pseudo_2ndSyl(period,:,:))))';
p = [];
for n = 1:size(r_real_vs_pseudo_2ndSyl,2)
    p(n) = bootstrapping(r_real_vs_pseudo_2ndSyl(:,n),10000);
end
[~,~,~,adj_p]=fdr_bh(p);

figure(2); clf;
cmin = -0.05; cmax = 0.05;
topoplot(mean(mean(r_real_vs_pseudo_2ndSyl(period,:,:)),3), 'channel_loc.elp', ...
    'maplimits', [cmin cmax],'emarker2',{find(adj_p<0.05),'*','w',16,3}); 
colormap(jet); caxis([cmin cmax]);

savename = sprintf('%s//r_real_vs_pseudo_%s_%s_2ndSyl.mat',folder,features,FreqRange);
save(savename,'r_real_vs_pseudo_2ndSyl');




function [r_real_1stSyl, r_real_2ndSyl, r_pseudo_1stSyl, r_pseudo_2ndSyl] = ...
    reconstrAcc_1st2ndSyl_realwords_vs_pseudowords_singleSubject...
    (subjectId, pred, eeg_partial, freqRange, stim1stSyl, stim2ndSyl, partial_or_not)

trlNum = 64;
eegfilename = sprintf('%s//%s_eeg_%s.mat',folder,subjectId,'raw');
b = load(eegfilename);
eeg = b.eeg;

% get IDs of trials with artefacts
artefactID = [];
for tr = 1:length(eeg.data)
    eeg_check4artefact = eeg.data(3,:);
    if abs(max(max(eeg_check4artefact{tr})))>35
        artefactID = [artefactID, tr];
    end
end
artefactID
artefactID_real = artefactID(artefactID<trlNum+1);
artefactID_pseudo = artefactID((artefactID>trlNum)&(artefactID<2*trlNum+1))-trlNum;

% TRF
% Normalising EEG data
eeg = eeg.data(freqRange,:);

for tr = 1:length(eeg) % decimation
    eeg{tr} = eeg{tr}(256*0.25+1:end,:);
    eeg_tr = [];
    for c = 1:size(eeg{tr},2)
        eeg_tr_c = decimate(eeg{tr}(:,c),decimate_factor,'fir');
        eeg_tr(:,c) = eeg_tr_c-mean(eeg_tr_c);
         %eeg_tr(:,c) = (eeg_tr_c-mean(eeg_tr_c))/std(eeg_tr_c);
    end
    eeg{tr} = eeg_tr;
end

splitGroupNum = 5;
realPos = 1:trlNum;
pseudoPos = trlNum+1:2*trlNum;

% real
stim1stSyl_real = stim1stSyl(realPos);
stim1stSyl_real(artefactID_real) = [];
stim2ndSyl_real = stim2ndSyl(realPos);
stim2ndSyl_real(artefactID_real) = [];
if partial_or_not
    eeg_real = eeg_partial.real;
else
    eeg_real = eeg(realPos);
    eeg_real(artefactID_real) = [];
end
trl4testId = splitgroup(1:length(eeg_real),splitGroupNum);
r_real_1stSyl = [];
r_real_2ndSyl = [];
for tt = 1:length(trl4testId)
    eeg_real4test = eeg_real(trl4testId{tt});  
    stim1stSyl_real4test = stim1stSyl_real(trl4testId{tt});
    stim2ndSyl_real4test = stim2ndSyl_real(trl4testId{tt});
    r_1stSyl = [];
    r_2ndSyl = [];
    for ss = 1:length(eeg_real4test) % loop testing trials
        pos_1stSyl = round((find(stim1stSyl_real4test{ss}==1))/decimate_factor);        
        for pp = 1:length(pos_1stSyl)
            for lag = 1:lag_pts
                pred_real4test = pred{lag}(tt).real;
                %pred_real4test_ss = (pred_real4test{ss}-repmat(mean(pred_real4test{ss}),[size(pred_real4test{ss},1),1]))./repmat(std(pred_real4test{ss}),[size(pred_real4test{ss},1),1]);
                pred_real4test_ss = pred_real4test{ss}-repmat(mean(pred_real4test{ss}),[size(pred_real4test{ss},1),1]);
                corr_period = pos_1stSyl(pp)+lag:pos_1stSyl(pp)+lag-1+corr_dur;
                if max(corr_period)>129,continue;end
                r_1stSyl(lag,ss,pp,:) = Fisher(PearsonCorr(eeg_real4test{ss}(corr_period,:),pred_real4test_ss(corr_period,:)));
                % [(lags)x(trials)x(syllable pos)x(electrodes)]
            end
        end
        pos_2ndSyl = round((find(stim2ndSyl_real4test{ss}==1))/decimate_factor);        
        for pp = 1:length(pos_2ndSyl)
            for lag = 1:lag_pts
                pred_real4test = pred{lag}(tt).real;
                %pred_real4test_ss = (pred_real4test{ss}-repmat(mean(pred_real4test{ss}),[size(pred_real4test{ss},1),1]))./repmat(std(pred_real4test{ss}),[size(pred_real4test{ss},1),1]);
                pred_real4test_ss = pred_real4test{ss}-repmat(mean(pred_real4test{ss}),[size(pred_real4test{ss},1),1]);
                corr_period = pos_2ndSyl(pp)+lag:pos_2ndSyl(pp)+lag-1+corr_dur;
                if max(corr_period)>129,continue;end
                r_2ndSyl(lag,ss,pp,:) = Fisher(PearsonCorr(eeg_real4test{ss}(corr_period,:),pred_real4test_ss(corr_period,:)));
                % [(lags)x(trials)x(syllable pos)x(electrodes)]
            end
        end
    end
    r_real_1stSyl(:,:,tt) = squeeze(mean(mean(r_1stSyl,2),3));
    r_real_2ndSyl(:,:,tt) = squeeze(mean(mean(r_2ndSyl,2),3));
end
r_real_1stSyl = mean(r_real_1stSyl,3);
r_real_2ndSyl = mean(r_real_2ndSyl,3);

% pseudo
stim1stSyl_pseudo = stim1stSyl(pseudoPos);
stim1stSyl_pseudo(artefactID_pseudo) = [];
stim2ndSyl_pseudo = stim2ndSyl(pseudoPos);
stim2ndSyl_pseudo(artefactID_pseudo) = [];
if partial_or_not
    eeg_pseudo = eeg_partial.pseudo;
else
    eeg_pseudo = eeg(pseudoPos);
    eeg_pseudo(artefactID_pseudo) = [];
end
trl4testId = splitgroup(1:length(eeg_pseudo),splitGroupNum);
r_pseudo_1stSyl = [];
r_pseudo_2ndSyl = [];
for tt = 1:length(trl4testId)
    eeg_pseudo4test = eeg_pseudo(trl4testId{tt});  
    stim1stSyl_pseudo4test = stim1stSyl_pseudo(trl4testId{tt});
    stim2ndSyl_pseudo4test = stim2ndSyl_pseudo(trl4testId{tt});
    r_1stSyl = [];
    r_2ndSyl = [];
    for ss = 1:length(eeg_pseudo4test)
        pos_1stSyl = round((find(stim1stSyl_pseudo4test{ss}==1))/decimate_factor);        
        for pp = 1:length(pos_1stSyl)
            for lag = 1:lag_pts
                pred_pseudo4test = pred{lag}(tt).pseudo;
                %pred_pseudo4test_ss = (pred_pseudo4test{ss}-repmat(mean(pred_pseudo4test{ss}),[size(pred_pseudo4test{ss},1),1]))./repmat(std(pred_pseudo4test{ss}),[size(pred_pseudo4test{ss},1),1]);
                pred_pseudo4test_ss = pred_pseudo4test{ss}-repmat(mean(pred_pseudo4test{ss}),[size(pred_pseudo4test{ss},1),1]);
                corr_period = pos_1stSyl(pp)+lag:pos_1stSyl(pp)+lag-1+corr_dur;
                if max(corr_period)>129,continue;end
                r_1stSyl(lag,ss,pp,:) = Fisher(PearsonCorr(eeg_pseudo4test{ss}(corr_period,:),pred_pseudo4test_ss(corr_period,:)));
                % [(lags)x(trials)x(syllable pos)x(electrodes)]
            end
        end
        pos_2ndSyl = round((find(stim2ndSyl_pseudo4test{ss}==1))/decimate_factor);        
        for pp = 1:length(pos_2ndSyl)
            for lag = 1:lag_pts
                pred_pseudo4test = pred{lag}(tt).pseudo;
                %pred_pseudo4test_ss = (pred_pseudo4test{ss}-repmat(mean(pred_pseudo4test{ss}),[size(pred_pseudo4test{ss},1),1]))./repmat(std(pred_pseudo4test{ss}),[size(pred_pseudo4test{ss},1),1]);
                pred_pseudo4test_ss = pred_pseudo4test{ss}-repmat(mean(pred_pseudo4test{ss}),[size(pred_pseudo4test{ss},1),1]);
                corr_period = pos_2ndSyl(pp)+lag:pos_2ndSyl(pp)+lag-1+corr_dur;
                if max(corr_period)>129,continue;end
                r_2ndSyl(lag,ss,pp,:) = Fisher(PearsonCorr(eeg_pseudo4test{ss}(corr_period,:),pred_pseudo4test_ss(corr_period,:)));
                % [(lags)x(trials)x(syllable pos)x(electrodes)]
            end
        end
    end
    r_pseudo_1stSyl(:,:,tt) = squeeze(mean(mean(r_1stSyl,2),3));
    r_pseudo_2ndSyl(:,:,tt) = squeeze(mean(mean(r_2ndSyl,2),3));
end
r_pseudo_1stSyl = mean(r_pseudo_1stSyl,3);
r_pseudo_2ndSyl = mean(r_pseudo_2ndSyl,3);


end

function subgroupIDs = splitgroup(groupIDs,splitGroupNum)
    subgroupIDs = {};
    lenOfEachSubGroup = round(length(groupIDs)/splitGroupNum);
    for g = 1:splitGroupNum
        if g==splitGroupNum
            subgroupOrder = (g-1)*lenOfEachSubGroup+1:length(groupIDs);
        else
            subgroupOrder = (g-1)*lenOfEachSubGroup+1:g*lenOfEachSubGroup;
        end
        subgroupIDs{g} = groupIDs(subgroupOrder);
    end
end

function r = PearsonCorr(aaa,bbb)
%aaa = aaa - repmat(mean(aaa),[size(aaa,1),1]);
%bbb = bbb - repmat(mean(bbb),[size(bbb,1),1]);
r = (sum(aaa.*bbb))./(((sum(aaa.^2)).*(sum(bbb.^2))).^0.5);
% for n = 1:size(aaa,2)
%     r(:,n) = corr(aaa(:,n),bbb(:,n));
% end
end


end