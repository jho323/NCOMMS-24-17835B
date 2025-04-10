clear
pathName = 'C:\Users\jch83\OneDrive - University of Pittsburgh\Documents\Research - Proprioception - VOP\Data\Thalamic Density Mapping\';


% fileName = 'MER_23_36_LT_mergedData.mat';
fileName = 'MER_23_36_RT_mergedData.mat';


singleStruct = load([pathName fileName]);
nameStruct = (fieldnames(singleStruct));

data = singleStruct.(nameStruct{1});
%% find stim depths
close all

fields = fieldnames(data);
ind = find(contains(fields,'ONSTIM'));
count = 0;
timeStamps = {};
freqLims = [.5 4];
%     4 8
%     8 13
%     13 30
%     30 100]; %Delta 0.5-4Hz, Theta 4-8Hz, alpha 8-13Hz, beta 14-30, gamma 31+
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

for f = 1:size(freqLims,1)
    count = 0;
    for t = ind'
%         if t==28
%             continue
%         end

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

            centralSig = singleDat.('CMacro_RAW_01___Central');
            posteriorSig = singleDat.('CMacro_RAW_02___Posterior');
            medialSig = singleDat.('CMacro_RAW_03___Medial');

            macroTrack{count}(1,aa) = peak2peak(centralSig(startPt:startPt+lenghtPt/2));
            macroTrack{count}(2,aa) = peak2peak(posteriorSig(startPt:startPt+lenghtPt/2));
            macroTrack{count}(3,aa) = peak2peak(medialSig(startPt:startPt+lenghtPt/2));

        end

        count2 = 0;
        for tt = [1:14]%ceegInds'
            trackct = [0 0 0];
            count2 = count2+1;
            cEEG = double(singleDat.(subFields{tt}));

            %         figure(t)
            %         subplot(length([1:14]),1,count2)
            %         hold on
            %         xline(startPoints,'--g')
            %         xline(endPoints,'--r')
            %         plot(time,cEEG)
            %         title(string(subFields{tt}))
            %         hold off
            %         sgtitle(depthLabs{count})

            %         preWindowPower = [];
            %         stimWindowPower = [];
            %         postWindowPower = [];

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

                [s1,f1,t1] =  pspectrum(preWindow,Fs,'spectrogram','FrequencyLimits',freqLims(f,:),TimeResolution=TimeRes,OverlapPercent=Overlap);
                preWindowPower{f}{count2,count}{idx}(:,:,trackct(idx)) = 10*log10(s1);

                [s2,f2,t2] =  pspectrum(stimWindow,Fs,'spectrogram','FrequencyLimits',freqLims(f,:),TimeResolution=TimeRes,OverlapPercent=Overlap);
                stimWindowPower{f}{count2,count}{idx}(:,:,trackct(idx)) = 10*log10(s2);
                %             p2pWindow{count2,count}(aa) = peak2peak(stimWindow);

                [s3,f3,t3] = pspectrum(postWindow,Fs,'spectrogram','FrequencyLimits',freqLims(f,:),TimeResolution=TimeRes,OverlapPercent=Overlap);
                postWindowPower{f}{count2,count}{idx}(:,:,trackct(idx)) = 10*log10(s3);

%                   [s1,f1,t1] =  pspectrum(preWindow,Fs,'spectrogram','FrequencyLimits',freqLims(f,:),'FrequencyResolution',750,OverlapPercent=Overlap);
%                 preWindowPower{f}{count2,count}{idx}(:,:,trackct(idx)) = 10*log10(s1);
% 
%                 [s2,f2,t2] =  pspectrum(stimWindow,Fs,'spectrogram','FrequencyLimits',freqLims(f,:),'FrequencyResolution',750,OverlapPercent=Overlap);
%                 stimWindowPower{f}{count2,count}{idx}(:,:,trackct(idx)) = 10*log10(s2);
%             
%                 [s3,f3,t3] = pspectrum(postWindow,Fs,'spectrogram','FrequencyLimits',freqLims(f,:),'FrequencyResolution',750,OverlapPercent=Overlap);
%                 postWindowPower{f}{count2,count}{idx}(:,:,trackct(idx)) = 10*log10(s3);

            end

            for ttt = 1:length(trackLabs)
                preMean = mean(preWindowPower{f}{count2,count}{ttt},3);
                stimMean = mean(stimWindowPower{f}{count2,count}{ttt},3);
                postMean = mean(postWindowPower{f}{count2,count}{ttt},3);

                preComp = mean(preMean,1);
                stimComp = mean(stimMean,1);
                postComp = mean(postMean,1);

                meanPower{f}{count2,count}{ttt,1} = preComp';
                meanPower{f}{count2,count}{ttt,2} = stimComp';
                meanPower{f}{count2,count}{ttt,3} = postComp';

                 powerMat = [preComp';stimComp';postComp'];
                p = [ones(size(preComp')); ones(size(stimComp'))*2; ones(size(postComp'))*3];


                %             figure(count2)
                %
                %             subplot(length(trackLabs),length(ind),count+2*(ttt-1))
                %             boxplot(powerMat,p,'Labels',windowLabs,'symbol','+r')
                %             title([trackLabs{ttt} ' ' convertStringsToChars(string(singleDat.depth))])
                %             sgtitle(convertStringsToChars(string(subFields{tt})))
                %

                trialMat = [1 3];

                %             alpha = .05;
                %             offset = 0;
                %             for comb = 1:size(trialMat,1)
                %                 trial1 = trialMat(comb,1);
                %                 trial2 = trialMat(comb,2);
                %
                %                 [CI,sig]=bootstrapCompMeans(meanPower{count2,count}{ttt}(:,trial1),meanPower{count2,count}{ttt}(:,trial2),10000,alpha/size(trialMat,1));
                %
                %                 if sig
                %                     yt = get(gca, 'YTick');
                %                     axis([xlim    0  ceil(max(yt)*1.3)])
                %                     annotation('textbox',[.1 .9 .1 .1],'String',"alpha = "+alpha)
                %                     xt = get(gca, 'XTick');
                %                     hold on
                %                     plot(xt([trial1 trial2]), [1 1]*max(yt)*(1.1+offset), '-k',  mean(xt([trial1 trial2])), max(yt)*(1.15+offset), '*k')
                %                     hold off
                %                 end
                %                 offset = offset+0.001;
                %             end
            end

        end
    end
end

%% Comp between depths 
close all
tracker = {};
percentVar = [];
freqsToPlot= [4 5];

% ipsilateral leads: 2 3 6 7 8 11 12 13 14
eegToPlot= [1:14];
count2 = 0;
count3 = 0;

%1 - single Depth, 2 - across freqs 3- line plots
plotStyle = 2;

if plotStyle==2 || plotStyle==3
    freqsToPlot = [1 2 3 4 5];
end

for ff = freqsToPlot
    for tt = 1:size(meanPower{ff}{1,1},2) %compare tract
        if plotStyle==1
            windowsToPlot = 1:size(meanPower{ff}{1,1}{1,tt},2);
        else
            windowsToPlot=1; %1: pre, 2: stim, 3: post
        end

        for aa = windowsToPlot %compare windows
            tracker{tt,ff} =  [0 0 0];
            singlepercentVar = [];
            singleMat = [];
            s = [];
            compDepthsVals = [];
            count3 = 0;
        
            for t = eegToPlot
                count3 = count3+1;
                for ttt = 1:size(meanPower{ff},2)
                    singleTrial = meanPower{ff}{t,ttt}{tt,aa};

%                      if aa~=1 
% %                     if ttt==1 && aa~=1   
% %                         singleTrial = singleTrial-powerDif;
%                         singleTrial = singleTrial-median(meanPower{ff}{t,ttt}{tt}(:,1));
%                         display('hi')
%                     end

                    singleMat = [singleMat;singleTrial];
                    s = [s; ones(size(singleTrial))*ttt];
                    compDepthsVals(count3,ttt) = median(singleTrial);
                end

                [maxV, idxD] = max(compDepthsVals(count3,:));

                singlepercentVar(end+1) = (compDepthsVals(1)-compDepthsVals(2))/compDepthsVals(1)*100;

              
%                 if idxD==1
%                     tracker{tt,ff}(idxD) = tracker{tt,ff}(idxD)+1;
%                 elseif idxD==2
%                     tracker{tt,ff}(idxD) = tracker{tt,ff}(idxD)+1;
%                 end
               
                lineColor = [.7 .7 .7];


                if plotStyle==1
                    figure(count3+16*count2)
                    subplot(size(meanPower{ff}{t},2),size(meanPower{ff}{t}{1,tt},2),aa+3*(tt-1))
                    boxplot(singleMat,s,'Labels',depthLabs,'symbol','+r')
                    title([trackLabs{tt} ' ' windowLabs{aa}])
                    sgtitle([convertStringsToChars(string(subFields{t})) ' ' freqLabs{ff}])
                
                elseif plotStyle ==3
%                     trialMat = [1 2];
%                     alpha = .05;
%                     offset = 0;
%                     for comb = 1:size(trialMat,1)
%                         trial1 = trialMat(comb,1);
%                         trial2 = trialMat(comb,2);
%     
%                         [CI,sig]=bootstrapCompMeans(meanPower{ff}{t,trial1}{tt}(:,aa),meanPower{ff}{t,trial2}{tt}(:,aa),10000,alpha/size(trialMat,1));
%     
%                         if sig
%                             lineColor = [0 0 0];
%                         else
%                             
%                         end
%                     end

                    figure(1)
                    subplot(3,length(freqsToPlot),ff+length(freqsToPlot)*(tt-1))
                    hold on 
                    plot([1 2], [compDepthsVals(count3,1) compDepthsVals(count3,2)],'-o','color',lineColor)
                    
                    xticks([1 2])
                    xticklabels({'4' '0'})
                    title([trackLabs{tt}  ' ' freqLabs{ff}])
                    sgtitle(['MER 36 RT ' windowLabs{aa}])
                    hold off

                else
                    figure(count3)
                    subplot(3,length(freqsToPlot),ff+length(freqsToPlot)*(tt-1))
                    boxplot(singleMat,s,'Labels',depthLabs,'symbol','+r')
                    title([trackLabs{tt} ' ' windowLabs{aa} ' ' freqLabs{ff}])
                    sgtitle([convertStringsToChars(string(subFields{t})) ])
                end

                trialMat = [1 2
                    ];

                alpha = .05;
                offset = 0;
                for comb = 1:size(trialMat,1)
                    trial1 = trialMat(comb,1);
                    trial2 = trialMat(comb,2);

                    [CI,sig]=bootstrapCompMeans(meanPower{ff}{t,trial1}{tt,aa},meanPower{ff}{t,trial2}{tt,aa},10000,alpha/size(trialMat,1));

                    if sig
                        yt = get(gca, 'YTick');
%                         axis([xlim    -1  ceil(max(yt)*1.3)])
                        annotation('textbox',[.1 .9 .1 .1],'String',"alpha = "+alpha)
                        xt = get(gca, 'XTick');
                        hold on
                        plot(xt([trial1 trial2]), [1 1]*max(yt)*(1.1+offset), '-k',  mean(xt([trial1 trial2])), max(yt)*(1.15+offset), '*k')
                        hold off
                        
                        if idxD==1
                            tracker{tt,ff}(idxD) = tracker{tt,ff}(idxD)+1;
                        elseif idxD==2
                            tracker{tt,ff}(idxD) = tracker{tt,ff}(idxD)+1;
                        end

                    else
                        tracker{tt,ff}(3) = tracker{tt,ff}(3)+1;
                    end
                    
                    offset = offset+0.001;
                end
                percentVar(tt,aa) =  mean(singlepercentVar);
            end

            if plotStyle==3
                subplot(3,length(freqsToPlot),ff+length(freqsToPlot)*(tt-1))
                hold on 
                plot([1 2],[mean(compDepthsVals,1)],'-o','LineWidth',3,'color',[1 0 0])
                hold off
        
            end
            
        end
        
    end
    count2 = count2+1;
end


%% Plot Cell Density

firingRateLabs = {};
firingRateMat = [];
skip = 1;
count = 1;

param = 1;
axes = [];
spikeFs = 44*1000;

stimDepth = [0 4]; 

for t = 1:length(fields)
    close all
    depth = fields{t};

    if contains(depth,'ONSTIM')
        continue
    end

    singleDepth = data.(depth);
    subFields = fieldnames(singleDepth);
    spkInds = find(contains(subFields,'SPK'));
    firingDat = [];


    for aa = 1:length(spkInds)
        chName = subFields{spkInds(aa)};
        if length(singleDepth.(chName)) < 3
            continue
            skip = 0;
        end
        firingDat(aa,:) = singleDepth.(chName);
    end

    if ~isempty(firingDat)
        depthName = singleDepth.depth;
        firingRateLabs{count,1} = depthName(5:end);

        depthLab = round(str2double(firingRateLabs{count,1}));

        if ismember(depthLab,stimDepth)

            figure;
            axes = [];
            for ii = 1:size(firingDat,1)
                % Filter and Derivative
                filtDat = (gradient(notchFilter_LFP(firingDat(ii,:))));
                spikeThres = std(filtDat);
                spikeProm = spikeThres/2;
                
                findpeakVals= [spikeProm spikeThres*4 50
                               spikeProm spikeThres*5 50
                               spikeProm*2 spikeThres*3 50];

                if ismember(count, [18 21 25 26 27 28])
                    param = 2;
                elseif ismember(count, [11 13 18])
                    param = 3;
                elseif ismember(count, [])
                    param = 4;
                end
    
                [pks,loc] = findpeaks(filtDat,'MinPeakProminence',findpeakVals(param,1),'MinPeakHeight',findpeakVals(param,2),'MinPeakDistance',findpeakVals(param,3));
                fireRate = length(pks)/(length(filtDat)/spikeFs);
                firingRateMat{ii}(count,3) = fireRate;
                firingRateMat{ii}(count,1) = depthLab;
                firingRateMat{ii}(count,2) = str2double(firingRateLabs{count,1});
                
                time = 0:1/spikeFs:(length(filtDat)-1)/spikeFs;
    
                ax = subplot(size(firingDat,1),1,ii);
                hold on
                plot(loc/(spikeFs),pks,'o')
                plot(time,(filtDat))
                sgtitle(string(count) + ' ' +depthName)
                annotation('textbox',[.1 .9/(ii) .1 .1],'String',"fire Rate = "+fireRate)
                hold off
                axes(end+1) = ax;
            end
            %         linkaxes(axes,'xy')
            %         clear axes

                param = 1;
                count = count+1;
            
        end
    end
    
end

 firingRateMat{ii}(20,:) = [];

% close all
% channelLab = {'Center','Posterior','Medial'};
% figure;
% for aa = 1:size(firingRateMat,2)
%     singleElec = firingRateMat{aa};
% 
%     dep = unique(singleElec(:,1));
%     for tt = 1:length(dep)
% 
%         vals = singleElec(singleElec(:,1)==dep(tt),3);
%         depthVals(tt,:) = [dep(tt),mean(vals)];
% 
%     end
%     sortedMat{aa} = sortrows(depthVals,1,'ascend');
% 
% 
%     subplot(1,size(firingRateMat,2),aa)
%     barh(sortedMat{aa}(:,2)')
%     yticks(1:length(sortedMat{aa}(:,2)))
%     yticklabels(string(num2cell(sortedMat{aa}(:,1))))
%     title([channelLab{aa} ' Density Map'])
%     ylabel('Depth (mm)')
%     xlabel('spike/sec')
% 
% end

% close all

subVals={};
meanVals = [];


figure;
for aa = 1:size(firingRateMat,2)
%     singleMeanElec = sortedMat{aa};
    singelElecAll = firingRateMat{aa};
    error = [];


    for tt = 1:length(stimDepth)
        meanVals(aa,tt) = mean(singelElecAll(singelElecAll(:,1)==stimDepth(tt),3));
        subVals{aa}{tt} = singelElecAll(singelElecAll(:,1)==stimDepth(tt),3);
        error(tt) = std(subVals{aa}{tt})/(length(subVals{aa}{tt})^0.5);

    end

    subplot(1,size(firingRateMat,2),aa)
    bar(meanVals(aa,:))
    xticklabels(string(num2cell(stimDepth)))


    x = [1 2];

    hold on
    errorbar([1 2],meanVals(aa,:),error,'k','linestyle','none');
    for tt = 1:length(stimDepth)
        f = swarmchart(tt,subVals{aa}{tt},1000,'red','.');
    end
    hold off

    trialMat = [1 2];

    alpha = .05;
    offset = 0;
    for comb = 1:size(trialMat,1)
        trial1 = trialMat(comb,1);
        trial2 = trialMat(comb,2);

        [CI,sig]=bootstrapCompMeans(subVals{aa}{trial1},subVals{aa}{trial2},10000,alpha/size(trialMat,1));

        if sig
            yt = get(gca, 'YTick');
            axis([xlim    0  ceil(max(yt)*1.3)])
            annotation('textbox',[.1 .9 .1 .1],'String',"alpha = "+alpha)
            xt = get(gca, 'XTick');
            hold on
            plot(xt([trial1 trial2]), [1 1]*max(yt)*(1.1+offset), '-k',  mean(xt([trial1 trial2])), max(yt)*(1.15+offset), '*k')
            hold off
        end
        offset = offset+0.001;
    end
end

%% FUNCTIONS

function waterplot(s,f,t)
% Waterfall plot of spectrogram
%     waterfall(f,t,abs(s)'.^2)
%     set(gca,XDir="reverse",View=[30 50])

    imagesc(t,f,pow2db(s))
    axis xy;
    ylabel("Frequency (Hz)")
    xlabel("Time (s)")
    h = colorbar;
    h.Label.String = 'Power (db)';
end


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

filterbands=[300,3000];
[z,p,k] = butter(2, filterbands/(fs_chan/2), 'bandpass');
[sos,g] = zp2sos(z, p, k, 'down', 'two');


filtered_LFP=filtfilt(sos,g,double(LFP'))';
filtered_LFP=filtfilt(sos_line,g_line,double(filtered_LFP'))';
end   