function [model,bestlambda,pred,predstats,stimFeaName,tmin,tmax] = mTRF_SingleSubject...
    (folder, subjectID, freqRange, comparison, dirTRF, stimIdx, ica_or_raw, downsamplefactor,trialpercent)


fs = 256/downsamplefactor;

% TRF hyperparameters
tmin = 0;
tmax = 500;
lambdapower = -4:0.5:4;
lambdas = [];
for p = 1:length(lambdapower)
    lambdas(p) = 10^(lambdapower(p));
end

trlNum = 64;
itnum = 100;

% lambdas = [1e-6,1e-3,1e-4,1e-3,1e-2,1e-1,1e0,1,1e2,1e3,1e4];
% dirTRF = 1; % Forward TRF model
% stimIdx = 1; % 1: env; 2: word onset

% Loading Stim and EEG data
switch comparison
    case 1 % real vs pseudo vs back
        Comparison = 'RealPseudoBack';
    otherwise % real vs pseudo
        Comparison = 'RealPseudo';
end

stimfilename = sprintf('%s//dataStim_%s_norm.mat',folder,Comparison);
a = load(stimfilename);
stim = a.stim;
stimfea = stim.data(stimIdx,:);
size(stimfea)
stimFea = {};
for tr = 1:size(stimfea,2) % decimation
    stimFea{tr} = [];
    for fea = 1:size(stimfea,1)
        if (stimIdx(fea)>8)&(stimIdx(fea)<12)&(comparison==2) % syllable/word onsets
            onesPos = ceil((find(stimfea{fea,tr}==1))/downsamplefactor);
            stimfea_fea_tr = zeros(length(stimfea{fea,tr}),1);
            stimfea_fea_tr = stimfea_fea_tr(1:downsamplefactor:end);
            stimfea_fea_tr(onesPos)=1;
        elseif (stimIdx(fea)>6)&(stimIdx(fea)<9) % phonetics and phonemes
            stimfea_fea_tr = [];
            for c = 1:size(stimfea{fea,tr},2)
               stimfea_fea_tr(:,c) = stimfea{fea,tr}(1:downsamplefactor:end,c);
            end
        elseif ((stimIdx(fea)==9)|(stimIdx(fea)==10))&(comparison==1)
            stimfea_fea_tr = stimfea{fea,tr};
        else
            stimfea_fea_tr = [];
            for c = 1:size(stimfea{fea,tr},2)
               stimfea_fea_tr(:,c) = decimate(stimfea{fea,tr}(:,c),downsamplefactor,'fir');
            end      
        end
        stimFea{tr} = [stimFea{tr},stimfea_fea_tr];
    end
end
size(stimFea{tr})
stimFeaName = '';
for st = 1:length(stimIdx)
    stimFeaName = sprintf('%s%s',stimFeaName,stim.names{stimIdx(st)});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

eegfilename = sprintf('%s//%s_eeg_%s.mat',folder,subjectID,ica_or_raw);
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
artefactID_back = artefactID(artefactID>2*trlNum)-2*trlNum;

% TRF
% Normalising EEG data
eeg = eeg.data(freqRange,:);

% cateeg = eeg{1};
% for tr = 2:length(eeg) % getting all values
%     cateeg = cat(1,cateeg,eeg{tr});
% end
% normFactorEeg = std(cateeg(:));
for tr = 1:length(eeg) % decimation
    %eeg{tr} = eeg{tr}/normFactorEeg;
    eeg{tr} = eeg{tr}(256*0.25+1:end,:);
    eeg_tr = [];
    for c = 1:size(eeg{tr},2)
        eeg_tr(:,c) = decimate(eeg{tr}(:,c),downsamplefactor,'fir');
    end
    eeg{tr} = eeg_tr;
end

%% cross validation and testing
splitGroupNum = 5;
realPos = 1:trlNum;
pseudoPos = trlNum+1:2*trlNum;
backPos = 2*trlNum+1:3*trlNum;

model = [];
predstats = [];

% real
stimFea_real = stimFea(realPos);
stimFea_real(artefactID_real) = [];
eeg_real = eeg(realPos);
eeg_real(artefactID_real) = [];
[eeg_real,stimFea_real] = choose_eeg_stim_based_on_percent(eeg_real,stimFea_real,trialpercent);
trl4testId = splitgroup(1:length(eeg_real),splitGroupNum);
for tt = 1:length(trl4testId)
    % (1) determine the training and testing set
    trl4trainId = setdiff(1:length(eeg_real),trl4testId{tt});
    eeg_real4train = eeg_real(trl4trainId);
    eeg_real4test = eeg_real(trl4testId{tt});    
    stimFea_real4train = stimFea_real(trl4trainId);
    stimFea_real4test = stimFea_real(trl4testId{tt});
    % (2) cross-validation using the training set to determine optimal lambda
    [stats_real,t] = mTRFcrossval(stimFea_real4train,eeg_real4train,fs,dirTRF,tmin,tmax,lambdas,'verbose',0);
    [~,bestLambda_real] = max(squeeze(mean(mean(stats_real.r),3)));
    bestlambda(tt).real = lambdas(bestLambda_real);
    % (3) use the optimal lambda to obtain mTRF model
    model(tt).real.forTest = mTRFtrain(stimFea_real4train,eeg_real4train,fs,dirTRF,tmin,tmax,lambdas(bestLambda_real),'verbose',0);
    model(tt).real.forIllustration = mTRFtrain(stimFea_real4train,eeg_real4train,fs,dirTRF,-200,600,lambdas(bestLambda_real),'verbose',0);
    % (4) use the testing set and the mTRF model to obtain predictive r
    [pred(tt).real,predstats(tt).real] = mTRFpredict(stimFea_real4test,eeg_real4test,model(tt).real.forTest);
%      % (5) obtain the shuffled predictive r
%     predstats_shuffled_r = [];
%     predstats_shuffled_fisher_r = [];
%     predstats_shuffled_err = [];
%     for it = 1:itnum
%         while 1
%             randomseq = randperm(length(trl4testId{tt}));
%             if all(randomseq-trl4testId{tt}),break;end % if nonzero then skip
%         end
%         [~,predstats_shuffled_it] = mTRFpredict(stimFea_real4test,eeg_real4test(randomseq),model(tt).real.forTest,'verbose',0);
%         predstats_shuffled_r(:,:,it) = predstats_shuffled_it.r;
%         predstats_shuffled_fisher_r(:,:,it) = Fisher(predstats_shuffled_it.r);
%         predstats_shuffled_err(:,:,it) = predstats_shuffled_it.err;
%     end
%     predstats(tt).real_shuffled.r = mean(predstats_shuffled_r,3);
%     predstats(tt).real_shuffled.r_std = std(predstats_shuffled_r,0,3);
%     predstats(tt).real_shuffled.fisher_r = mean(predstats_shuffled_fisher_r,3);
%     predstats(tt).real_shuffled.fisher_r_std = std(predstats_shuffled_fisher_r,0,3);
%     predstats(tt).real_shuffled.err = mean(predstats_shuffled_err,3);
%     predstats(tt).real_shuffled.err_std = std(predstats_shuffled_err,0,3);
    tt
end
    
% pseudo   
stimFea_pseudo = stimFea(pseudoPos);
stimFea_pseudo(artefactID_pseudo) = [];
eeg_pseudo = eeg(pseudoPos);
eeg_pseudo(artefactID_pseudo) = [];
[eeg_pseudo,stimFea_pseudo] = choose_eeg_stim_based_on_percent(eeg_pseudo,stimFea_pseudo,trialpercent);
trl4testId = splitgroup(1:length(eeg_pseudo),splitGroupNum);
for tt = 1:length(trl4testId)
    % (1) determine the training and testing set
    trl4trainId = setdiff(1:length(eeg_pseudo),trl4testId{tt});
    eeg_pseudo4train = eeg_pseudo(trl4trainId);
    eeg_pseudo4test = eeg_pseudo(trl4testId{tt});
    stimFea_pseudo4train = stimFea_pseudo(trl4trainId);
    stimFea_pseudo4test = stimFea_pseudo(trl4testId{tt});
    % (2) cross-validation using the training set to determine optimal lambda
    [stats_pseudo,t] = mTRFcrossval(stimFea_pseudo4train,eeg_pseudo4train,fs,dirTRF,tmin,tmax,lambdas,'verbose',0);
    [~,bestLambda_pseudo] = max(squeeze(mean(mean(stats_pseudo.r),3)));
    bestlambda(tt).pseudo = lambdas(bestLambda_pseudo);
    % (3) use the optimal lambda to obtain mTRF model
    model(tt).pseudo.forTest = mTRFtrain(stimFea_pseudo4train,eeg_pseudo4train,fs,dirTRF,tmin,tmax,lambdas(bestLambda_pseudo),'verbose',0);
    model(tt).pseudo.forIllustration = mTRFtrain(stimFea_pseudo4train,eeg_pseudo4train,fs,dirTRF,-200,600,lambdas(bestLambda_pseudo),'verbose',0);
    % (4) use the testing set and the mTRF model to obtain predictive r
    [pred(tt).pseudo,predstats(tt).pseudo] = mTRFpredict(stimFea_pseudo4test,eeg_pseudo4test,model(tt).pseudo.forTest);
%      % (5) obtain the shuffled predictive r
%     predstats_shuffled_r = [];
%     predstats_shuffled_fisher_r = [];
%     predstats_shuffled_err = [];
%     for it = 1:itnum
%         while 1
%             randomseq = randperm(length(trl4testId{tt}));
%             if all(randomseq-trl4testId{tt}),break;end % if nonzero then skip
%         end
%         [~,predstats_shuffled_it] = mTRFpredict(stimFea_pseudo4test,eeg_pseudo4test(randomseq),model(tt).pseudo.forTest,'verbose',0);
%         predstats_shuffled_r(:,:,it) = predstats_shuffled_it.r;
%         predstats_shuffled_fisher_r(:,:,it) = Fisher(predstats_shuffled_it.r);
%         predstats_shuffled_err(:,:,it) = predstats_shuffled_it.err;
%     end
%     predstats(tt).pseudo_shuffled.r = mean(predstats_shuffled_r,3);
%     predstats(tt).pseudo_shuffled.r_std = std(predstats_shuffled_r,0,3);
%     predstats(tt).pseudo_shuffled.fisher_r = mean(predstats_shuffled_fisher_r,3);
%     predstats(tt).pseudo_shuffled.fisher_r_std = std(predstats_shuffled_fisher_r,0,3);
%     predstats(tt).pseudo_shuffled.err = mean(predstats_shuffled_err,3);
%     predstats(tt).pseudo_shuffled.err_std = std(predstats_shuffled_err,0,3);
    tt
end
    
if strcmp(Comparison,'RealPseudoBack')
    stimFea_back = stimFea(backPos);
    stimFea_back(artefactID_back) = [];
    eeg_back = eeg(backPos);
    eeg_back(artefactID_back) = [];
    [eeg_back,stimFea_back] = choose_eeg_stim_based_on_percent(eeg_back,stimFea_back,trialpercent);
    trl4testId = splitgroup(1:length(eeg_back),splitGroupNum);
    for tt = 1:length(trl4testId)
        % (1) determine the training and testing set
        trl4trainId = setdiff(1:length(eeg_back),trl4testId{tt});
        eeg_back4train = eeg_back(trl4trainId);
        eeg_back4test = eeg_back(trl4testId{tt});
        stimFea_back4train = stimFea_back(trl4trainId);
        stimFea_back4test = stimFea_back(trl4testId{tt});
        % (2) cross-validation using the training set to determine optimal lambda
        [stats_back,t] = mTRFcrossval(stimFea_back4train,eeg_back4train,fs,dirTRF,tmin,tmax,lambdas,'verbose',0);
        [~,bestLambda_back] = max(squeeze(mean(mean(stats_back.r),3)));
        bestlambda(tt).back = lambdas(bestLambda_back);
        % (3) use the optimal lambda to obtain mTRF model
        model(tt).back.forTest = mTRFtrain(stimFea_back4train,eeg_back4train,fs,dirTRF,tmin,tmax,lambdas(bestLambda_back),'verbose',0);
        model(tt).back.forIllustration = mTRFtrain(stimFea_back4train,eeg_back4train,fs,dirTRF,-200,600,lambdas(bestLambda_back),'verbose',0);
        % (4) use the testing set and the mTRF model to obtain predictive r
        [pred(tt).back,predstats(tt).back] = mTRFpredict(stimFea_back4test,eeg_back4test,model(tt).back.forTest);
%         % (5) obtain the shuffled predictive r
%         predstats_shuffled_r = [];
%         predstats_shuffled_fisher_r = [];
%         predstats_shuffled_err = [];
%         for it = 1:itnum
%             while 1
%                 randomseq = randperm(length(trl4testId{tt}));
%                 if all(randomseq-trl4testId{tt}),break;end % if nonzero then skip
%             end
%             [~,predstats_shuffled_it] = mTRFpredict(stimFea_back4test,eeg_back4test(randomseq),model(tt).back.forTest,'verbose',0);
%             predstats_shuffled_r(:,:,it) = predstats_shuffled_it.r;
%             predstats_shuffled_fisher_r(:,:,it) = Fisher(predstats_shuffled_it.r);
%             predstats_shuffled_err(:,:,it) = predstats_shuffled_it.err;
%         end
%         predstats(tt).back_shuffled.r = mean(predstats_shuffled_r,3);
%         predstats(tt).back_shuffled.r_std = std(predstats_shuffled_r,0,3);
%         predstats(tt).back_shuffled.fisher_r = mean(predstats_shuffled_fisher_r,3);
%         predstats(tt).back_shuffled.fisher_r_std = std(predstats_shuffled_fisher_r,0,3);
%         predstats(tt).back_shuffled.err = mean(predstats_shuffled_err,3);
%         predstats(tt).back_shuffled.err_std = std(predstats_shuffled_err,0,3);
        tt
    end
end

% modelAll(sub) = model;
%     
%     % Plot average TRF
%     normFlag = 1;
%     el = 1;
%     avgModel = mTRFmodelAvg(modelAll,normFlag);
%     % Plot avg TRF model
%     subplot(1,2,1)
%     % mTRFplot(model,'trf',[],el);
%     plot(avgModel.t,squeeze(avgModel.w))
%     title('Envelope avgTRF')
%     ylabel('Magnitude (a.u.)')
%     xlim([tmin+50,tmax-50])
%     ylim([-4,4])
%     axis square
%     run prepExport.m


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

function [eeg,stim] = choose_eeg_stim_based_on_percent(eeg,stim,trialpercentage)
    randomID = randperm(length(eeg));
    trialrejectID = randomID(1:round(length(eeg)*(1-trialpercentage)));
    eeg(trialrejectID) = [];
    stim(trialrejectID) = [];
end

end
