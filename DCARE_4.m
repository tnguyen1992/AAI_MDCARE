function DCARE_4
%% Extract hbo time series from dc variable in subject file
for dyad=[1:9 11:31 33:48 50:52 54:68] 
    load(sprintf(['P:/projects/DCARE/DCARE/MATLAB/procData/hmrData/DCARE_%02d_sub1.mat'],dyad));
    hbo_ch=hbo;
    load(sprintf(['P:/projects/DCARE/DCARE/MATLAB/procData/hmrData/DCARE_%02d_sub2.mat'],dyad));
    hbo_cg=hbo;
   
%% prepare matrix for wtc 
    wtc_all=zeros(16,3);
%     wtc_sep=zeros(16,8);    
    
%% Extract markes 
    [row,col] = find(s); %find start time frames of each condition from s (trigger variable)
    
% define durations of conditions (in time frames)
    duration_task=937;
    duration_rest=625;
%     duration_talk=2344;

%%  find period of interest
    hbo1=[t, hbo_ch(:,1)];
    hbo2=[t, hbo_cg(:,1)];

    [~,period,~,~,~]=wtc(hbo1,hbo2,'mcc',0,'ms',128); 
    period_low = find(period>10.0);                                          %period of interest for CARE: 10 - 50 s
    period_low = period_low(1);
    period_high = find(period>50.0);
    period_high = period_high(1);
    


%% calculate coherences for every channel
    for ch=1:16
        hbo1=[t, hbo_ch(:,ch)];
        hbo2=[t, hbo_cg(:,ch)];
        [Rsq, period,coi]=wtc(hbo1,hbo2,'mcc',0,'ms',128);
        %set values outside of coi to NaN
        for j=1:1:length(coi)
        Rsq(period >= coi(j), j) = NaN;
        end
% calculate mean activation in frequency band of interest
% collaboration condition
        %wtc_collaboration1 = mean(Rsq(:, markers(1):markers(1)+duration_task),2);
        wtc_collaboration1 = mean(mean(Rsq(period_low:period_high, row(1):row(1)+duration_task)));
        wtc_collaboration2 = mean(mean(Rsq(period_low:period_high, row(2):row(2)+duration_task)));

% individual condition
        wtc_individual1 = mean(mean(Rsq(period_low:period_high, row(3):row(3)+duration_task)));
        wtc_individual2 = mean(mean(Rsq(period_low:period_high, row(4):row(4)+duration_task)));

% resting phase
        wtc_rest1 = mean(mean(Rsq(period_low:period_high, row(5):row(5)+duration_rest)));
        wtc_rest2 = mean(mean(Rsq(period_low:period_high, row(6):row(6)+duration_rest)));
        wtc_rest3 = mean(mean(Rsq(period_low:period_high, row(7):row(7)+duration_rest)));

% talking condition
%         if length(row)<=15;
%             wtc_talk = NaN;
%         elseif length(s)-row(16)>=duration_talk;
%             wtc_talk = mean(mean(Rsq(period_low:period_high, row(16):row(16)+duration_talk)));
%         else 
%             wtc_talk = NaN;
%         end

%% calculate mean coherences

        wtc_collaboration_m = (wtc_collaboration1+wtc_collaboration2)/2;
        wtc_individual_m = (wtc_individual1+wtc_individual2)/2;
        wtc_rest_m = (wtc_rest1+wtc_rest2+wtc_rest3)/3;
        %% save values
        wtc_all(ch,1:3)=[wtc_collaboration_m, wtc_individual_m, wtc_rest_m];
%         wtc_sep(ch,1:8)=[wtc_collaboration1,wtc_collaboration2,wtc_individual1,wtc_individual2,wtc_rest1,wtc_rest2,wtc_rest3,wtc_talk];
    end

%%  save cohin for dyad in file
save(sprintf(['P:/projects/DCARE/DCARE/MATLAB/procData/wtcData/DCARE_%02d_wtc.mat'],dyad),'wtc_all', 'Rsq')
end
end