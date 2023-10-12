function [pred,predstats] = mTRF_subperiod_SingleSubject...
    (folder, subjectID, bestlambda, comparison, freqRange, stimIdx, fs, downsamplefactor, dirTRF, subperiods, ica_or_raw)

trlNum = 64;

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
stimFea = {};
for tr = 1:size(stimfea,2) % decimation
    stimFea{tr} = [];
    for fea = 1:size(stimfea,1)     
        if (stimIdx(fea)>8)&(stimIdx(fea)<12) % syllable/word onsets
            onesPos = ceil((find(stimfea{fea,tr}==1))/downsamplefactor);
            stimfea_fea_tr = zeros(length(stimfea{fea,tr}),1);
            stimfea_fea_tr = stimfea_fea_tr(1:downsamplefactor:end);
            stimfea_fea_tr(onesPos)=1;
        else
            stimfea_fea_tr = [];
            for c = 1:size(stimfea{fea,tr},2)
               stimfea_fea_tr(:,c) = decimate(stimfea{fea,tr}(:,c),downsamplefactor,'fir');
            end      
        end
        stimFea{tr} = [stimFea{tr},stimfea_fea_tr];
    end
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
realPos = 1:trlNum;
pseudoPos = trlNum+1:2*trlNum;
backPos = 2*trlNum+1:3*trlNum;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
splitGroupNum = 5;
pred = {};
predstats = {};

% real
stimFea_real = stimFea(realPos);
stimFea_real(artefactID_real) = [];
eeg_real = eeg(realPos);
eeg_real(artefactID_real) = [];
trl4testId = splitgroup(1:length(eeg_real),splitGroupNum);
for tt = 1:length(trl4testId)
    % (1) determine the training and testing set
    trl4trainId = setdiff(1:length(eeg_real),trl4testId{tt});
    eeg_real4train = eeg_real(trl4trainId);
    eeg_real4test = eeg_real(trl4testId{tt});
    stimFea_real4train = stimFea_real(trl4trainId);
    stimFea_real4test = stimFea_real(trl4testId{tt});
    % (2) use the optimal lambda to obtain mTRF model for the subperiods
    for p = 1:length(subperiods)
        tmin = subperiods{p}(1);
        tmax = subperiods{p}(2);
        model_real = mTRFtrain(stimFea_real4train,eeg_real4train,fs,dirTRF,tmin,tmax,bestlambda(tt).real,'verbose',0);
        % (3) use the testing set and the mTRF model to obtain predictive r
        [pred{p}(tt).real,predstats{p}(tt).real] = mTRFpredict(stimFea_real4test,eeg_real4test,model_real,'verbose',0);
    end
    tt
end

% pseudo
stimFea_pseudo = stimFea(pseudoPos);
stimFea_pseudo(artefactID_pseudo) = [];
eeg_pseudo = eeg(pseudoPos);
eeg_pseudo(artefactID_pseudo) = [];
trl4testId = splitgroup(1:length(eeg_pseudo),splitGroupNum);
for tt = 1:length(trl4testId)
    % (1) determine the training and testing set
    trl4trainId = setdiff(1:length(eeg_pseudo),trl4testId{tt});
    eeg_pseudo4train = eeg_pseudo(trl4trainId);
    eeg_pseudo4test = eeg_pseudo(trl4testId{tt});
    stimFea_pseudo4train = stimFea_pseudo(trl4trainId);
    stimFea_pseudo4test = stimFea_pseudo(trl4testId{tt});
    % (2) use the optimal lambda to obtain mTRF model for the subperiods
    for p = 1:length(subperiods)
        tmin = subperiods{p}(1);
        tmax = subperiods{p}(2);
        model_pseudo = mTRFtrain(stimFea_pseudo4train,eeg_pseudo4train,fs,dirTRF,tmin,tmax,bestlambda(tt).pseudo,'verbose',0);
        % (3) use the testing set and the mTRF model to obtain predictive r
        [pred{p}(tt).pseudo,predstats{p}(tt).pseudo] = mTRFpredict(stimFea_pseudo4test,eeg_pseudo4test,model_pseudo,'verbose',0);
    end
    tt
end
   
if strcmp(Comparison,'RealPseudoBack')
    % back
    stimFea_back = stimFea(backPos);
    stimFea_back(artefactID_back) = [];
    eeg_back = eeg(backPos);
    eeg_back(artefactID_back) = [];
    trl4testId = splitgroup(1:length(eeg_back),splitGroupNum);
    for tt = 1:length(trl4testId)
        % (1) determine the training and testing set
        trl4trainId = setdiff(1:length(eeg_back),trl4testId{tt});
        eeg_back4train = eeg_back(trl4trainId);
        eeg_back4test = eeg_back(trl4testId{tt});
        stimFea_back4train = stimFea_back(trl4trainId);
        stimFea_back4test = stimFea_back(trl4testId{tt});
        % (2) use the optimal lambda to obtain mTRF model for the subperiods
        for p = 1:length(subperiods)
            tmin = subperiods{p}(1);
            tmax = subperiods{p}(2);
            model_back = mTRFtrain(stimFea_back4train,eeg_back4train,fs,dirTRF,tmin,tmax,bestlambda(tt).back,'verbose',0);
            % (3) use the testing set and the mTRF model to obtain predictive r
            [pred{p}(tt).back,predstats{p}(tt).back] = mTRFpredict(stimFea_back4test,eeg_back4test,model_back,'verbose',0);
        end
        tt
    end
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

end
