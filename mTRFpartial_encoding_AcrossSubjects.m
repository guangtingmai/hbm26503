function mTRFpartial_encoding_AcrossSubjects(freqRange, comparison, TargetStimIdx, PartialStimIdx, ica_or_raw, ver)

folder = 'C:\Users\mszgm1\OneDrive - The University of Nottingham\LEL 2013';
addpath(genpath(folder));
subfolder = sprintf('%s//analyses 2022',folder);

subjectID1 = '03_20am';  subjectID2 = '03_20eve'; subjectID3 = '03_22pm';   subjectID4 = '03_22eve'; subjectID5 = '03_24am';  
subjectID6 = '03_24pm';  subjectID7 = '03_24eve'; subjectID8 = '03_26pm';   subjectID9 = '03_29am';  subjectID10 = '03_29pm'; 
subjectID11 = '03_30am'; subjectID12 = '03_30pm'; subjectID13 = '03_31am';  subjectID14 = '03_31pm'; subjectID15 = '04_02eve'; 
subjectID16 = '04_06am'; subjectID17 = '04_06pm'; subjectID18 = '04_06eve'; subjectID19 = '04_07am'; subjectID20 = '04_07pm';
subjectID = {subjectID1, subjectID2, subjectID3, subjectID4, subjectID5, subjectID6, subjectID7, subjectID8, subjectID9, subjectID10, ...
     subjectID11, subjectID12, subjectID13, subjectID14, subjectID15, subjectID16, subjectID17, subjectID18, subjectID19, subjectID20};


downsamplefactor = 4;
new_fs = 256/downsamplefactor;

switch comparison
    case 1 % real vs pseudo vs back
        Comparison = 'RealPseudoBack';
    otherwise % real vs pseudo
        Comparison = 'RealPseudo';
end
switch freqRange
    case 1 % Delta
        FreqRange = 'Delta';
    case 2 % Theta
        FreqRange = 'Theta';
    otherwise % DeltaTheta
        FreqRange = 'DeltaTheta';
end

model = {};
bestlambda = {};
pred = {};
predstats = {};
eeg_partial = {};
targetStimFea = {};
dirTRF = 1;
tmin = 0;
tmax = 500;
shuffle_or_not = 1;
for s = 1:length(subjectID)
    [model{s},bestlambda{s},pred{s},predstats{s},eeg_partial{s},targetStimFea{s},TargetStimname,PartialStimname] = ...
        mTRFpartial_SingleSubject(subfolder, subjectID{s}, freqRange, comparison, dirTRF, TargetStimIdx, PartialStimIdx, ...
        downsamplefactor, tmin, tmax, ica_or_raw, ver, shuffle_or_not);
    s
end

savename = sprintf('%s//mTRF_encoding_%s_%s_%s(target)_%s(partialed)%d_%d-%dms_eeg%s_%dHzFs.mat',subfolder,FreqRange,...
    Comparison,TargetStimname,PartialStimname,ver,tmin,tmax,ica_or_raw,new_fs);
save(savename,'model','bestlambda','pred','predstats','eeg_partial','targetStimFea');