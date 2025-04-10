clear
pathName = 'C:\Users\jch83\OneDrive - University of Pittsburgh\Documents\Research - Proprioception - VOP\Data\Thalamic Density Mapping\';


% fileName = 'MER_23_39_LT_mergedData.mat';
fileName = 'MER_23_39_RT_mergedData.mat';

%macro stim depth
stimDepth = [0 -3]; 

singleStruct = load([pathName fileName]);
nameStruct = (fieldnames(singleStruct));

data = singleStruct.(nameStruct{1});

%% find stim depths
close all

fields = fieldnames(data);
ind = find(contains(fields,'ONSTIM'));
count = 0;
timeStamps = {};
freqLims = [.5 4
    4 8
    8 13
    13 30
    30 100]; %Delta 0.5-4Hz, Theta 4-8Hz, alpha 8-13Hz, beta 14-30, gamma 31+
freqLabs = {'Delta','Theta','Alpha','Beta','Gamma'};

TimeRes = 1;
Overlap = 20;
preLength = 10;
stimLength = 30;
postLength = 20;

preWindowPower = {};
stimWindowPower = {};
postWindowPower = {};

p2pWindow = {};
meanPower = {};

windowLabs = {'Pre Stim', 'Stim', 'Post Stim'};
trackLabs = {'Center','Posterior','Medial'};
depthLabs = {};
trackct = [0 0 0];

for f = 1%:size(freqLims,1)
    count = 0;
    for t = ind'
        if t==28
            continue
        end

        count = count+1;
        singleDat = data.(fields{t});
        subFields = fieldnames(singleDat);
        Fs = singleDat.SR_kHz*1000;
        begin = singleDat.timeBegin;
        allTimePoints = (singleDat.stimMarker/(44*1000))-begin-1;
        allTimePoints = allTimePoints+1;
        depthLabs{count} = convertStringsToChars(string(singleDat.depth));

        %get timestamps
        startPoints = [allTimePoints(1)];
        endPoints = [allTimePoints(end)];
        sep = diff(allTimePoints);
        minSeps = min(unique(sep))*2;


        startPoints = sort([startPoints allTimePoints(find(sep>minSeps)+1)]);
        endPoints = sort([allTimePoints(find(sep>minSeps)) endPoints]);
        timeStamps{count} = [startPoints; endPoints];

        time = (1:length(singleDat.CEEG_1___02___F3))/Fs;

        ceegInds = find(contains(subFields,'CEEG'));

        for aa = 1:size(timeStamps{count},2)
            startPt = timeStamps{count}(1,aa)*44000;
            endPt = timeStamps{count}(2,aa)*44000;
            lenghtPt = endPt-startPt;

            centralSig = notchFilter_LFP(singleDat.('CMacro_RAW_01___Central'));
            posteriorSig = notchFilter_LFP(singleDat.('CMacro_RAW_02___Posterior'));
            medialSig = notchFilter_LFP(singleDat.('CMacro_RAW_03___Medial'));

            macroTrack{count}(1,aa) = peak2peak(centralSig(startPt:startPt+lenghtPt/2));
            macroTrack{count}(2,aa) = peak2peak(posteriorSig(startPt:startPt+lenghtPt/2));
            macroTrack{count}(3,aa) = peak2peak(medialSig(startPt:startPt+lenghtPt/2));

        end

        count2 = 0;
        for tt = [4:7]%ceegInds'
            trackct = [0 0 0];
            count2 = count2+1;
            cEEG = double(singleDat.(subFields{tt}));

            for aa = 1:size(timeStamps{count},2)
                startTime = floor(timeStamps{count}(1,aa)*Fs);
                endTime = floor(timeStamps{count}(2,aa)*Fs);

                [minV, idx] = min(macroTrack{count}(:,aa));

                if startTime+(stimLength*Fs) > endTime
                    display('STIM WINDOW TOO LONG')
                    continue
                end

                if endTime+(postLength*Fs) > length(cEEG)
                    display('POST WINDOW TOO LONG')
                    continue
                    %                     postWindow = cEEG(endTime:length(cEEG));
                else
                    postWindow = cEEG(endTime:endTime+(postLength*Fs));

                end

                stimWindow = cEEG(startTime:startTime+(stimLength*Fs));
                preWindow = cEEG(startTime-(preLength*Fs):startTime);


                if idx ==1
                    trackct(1) = trackct(1) + 1;
                elseif idx ==2
                    trackct(2) = trackct(2) + 1;
                elseif idx == 3
                    trackct(3) = trackct(3) + 1;
                end

                figure;

                ax1 = subplot (1,3,1);
                [plot_s1,plot_f1,plot_t1] =   pspectrum(preWindow,Fs,'spectrogram',TimeResolution=TimeRes,OverlapPercent=Overlap);
                imagesc(plot_t1, plot_f1, 10*log10(plot_s1));
                axis xy;
                title([trackLabs{(idx)} ' ' windowLabs{1}])
                ylim([0 100])
                xlabel('Time (s)');
                ylabel('Frequency (Hz)');
                colorbar;
                
                
                ax2 = subplot(1,3,2);
                [plot_s2,plot_f2,plot_t2] =  pspectrum(stimWindow,Fs,'spectrogram',TimeResolution=TimeRes,OverlapPercent=Overlap);
                imagesc(plot_t2, plot_f2, 10*log10(plot_s2));
                axis xy;
                title([trackLabs{(idx)} ' ' windowLabs{2}])
                ylim([0 100])
                xlabel('Time (s)');
                ylabel('Frequency (Hz)');
                colorbar;
                
                ax3 = subplot(1,3,3);
                [plot_s3,plot_f3,plot_t3] = pspectrum(postWindow,Fs,'spectrogram',TimeResolution=TimeRes,OverlapPercent=Overlap);
                imagesc(plot_t3, plot_f3, 10*log10(plot_s3));
                axis xy;
                title([trackLabs{(idx)} ' ' windowLabs{3}])
                ylim([0 100])
                xlabel('Time (s)');
                ylabel('Frequency (Hz)');
                colorbar;
                colormapeditor

                minPower = min([10*log10(plot_s1(:)); 10*log10(plot_s2(:)); 10*log10(plot_s3(:))]);
                maxPower = max([10*log10(plot_s1(:)); 10*log10(plot_s2(:)); 10*log10(plot_s3(:))]);
                caxis(ax1, [minPower maxPower]);
                caxis(ax2, [minPower maxPower]);
                caxis(ax3, [minPower maxPower]);

                linkaxes([ax1 ax2 ax3],'y')
                title([trackLabs{trackct(idx)} ' ' windowLabs{3}])
                sgtitle([convertStringsToChars(depthLabs{count}) ' ' convertStringsToChars(string(subFields{tt}))])
                
            end
        end
    end
end

%% FUNCTIONS

function [ci95, rejectNull] = bootstrapCompMeans(dataSet1, dataSet2, bootstrapReps,alpha)

sampMeans1 = nan(1,bootstrapReps);
sampMeans2 = nan(1,bootstrapReps);
diffSampMeans = nan(1,bootstrapReps);
for i=1:bootstrapReps
    % Resample from each dataset with replacement
    bootstrapSamp1 = randsample(dataSet1, length(dataSet1), true);
    bootstrapSamp2 = randsample(dataSet2, length(dataSet2), true);

    % Get means of both samples
    meanSamp1 = mean(bootstrapSamp1);
    sampMeans1(i) = meanSamp1;
    meanSamp2 = mean(bootstrapSamp2);
    sampMeans2(i) = meanSamp2;

    % Get the difference of the means
    diffMeans = meanSamp1 - meanSamp2;
    diffSampMeans(i) = diffMeans;
end

% Calculate confidence interval of difference of means
ci95 = quantile(diffSampMeans, [alpha/2, 1-(alpha/2)]);

% If ci95 contains 0 then don't reject null
if (ci95(1) <= 0) && (ci95(2) >= 0)
    rejectNull = false;
else
    rejectNull = true;
end

%     Plot histogram to confirm
%     figure;
%     hist(diffSampMeans, 100)

end

function filtered_LFP = notchFilter_LFP(LFP)
fs_chan = 44000;

filterbands_Line=[58,62];
[z, p, k] = butter(2, filterbands_Line/(fs_chan/2), 'stop');
[sos_line,g_line] = zp2sos(z, p, k, 'down', 'two');

% filterbands=[300,3000];
% [z,p,k] = butter(2, filterbands/(fs_chan/2), 'bandpass');
% [sos,g] = zp2sos(z, p, k, 'down', 'two');


% filtered_LFP=filtfilt(sos,g,double(LFP'))';
filtered_LFP=filtfilt(sos_line,g_line,double(LFP'))';
end