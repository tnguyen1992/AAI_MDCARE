%% Extract SPM information
% wo sind unsere Daten?
srcPath = 'P:\projects\DCARE\DCARE\rawData\';


%% welche Probanden sind in dem Ordner?
sourceList    = dir([srcPath, 'DCARE_*']);
sourceList    = struct2cell(sourceList);
sourceList    = sourceList(1,:);
% wie viele Probanden gibt es?
numOfSources  = length(sourceList);
numOfPart       = zeros(1, numOfSources);

  for i=1:1:numOfSources
    numOfPart(i)  = sscanf(sourceList{i}, ['DCARE_%d']);
  end
  
  %% GLM
  
    for i = numOfPart
  srcFolder   = strcat(srcPath, sprintf(['DCARE_%02d/'], i));
  
  Sub1SrcDir = strcat(srcFolder, sprintf(['Subject1/HbO/']));
  Sub2SrcDir = strcat(srcFolder, sprintf(['Subject2/HbO/']));
  
% -------------------------------------------------------------------------
% Wo sind die Daten im jeweiligen Ordner? 
% -------------------------------------------------------------------------

  Sub1_SPMFile = strcat(Sub1SrcDir, sprintf(['SPM.mat']));
  Sub2_SPMFile = strcat(Sub2SrcDir, sprintf(['SPM.mat']));
      
  beta_id(i,1)=i;

    load(Sub1_SPMFile)
    beta_coop(i,1:16)=SPM.beta(2,:);
    beta_ind(i,1:16)=SPM.beta(3,:);
    beta_rest(i,1:16)=SPM.beta(4,:);
    beta_preform(i,1:16)=SPM.beta(6,:);
    
    clear SPM
    load(Sub2_SPMFile)
    beta_coop(i,17:32)=SPM.beta(2,:);
    beta_ind(i,17:32)=SPM.beta(3,:);
    beta_rest(i,17:32)=SPM.beta(4,:);
    beta_preform(i,17:32)=SPM.beta(6,:);
    
    end
    
      %% WTC 
i=1;
for id = numOfPart
srcFolder   = strcat(srcPath, sprintf(['DCARE_%02d/'], id));
    
% -------------------------------------------------------------------------
% Wo sind die Daten im jeweiligen Ordner? 
% -------------------------------------------------------------------------

  WTCFile = strcat(srcFolder, sprintf(['WTC.mat']));
  load(WTCFile)

  for condition=1:1:4
        
      if condition==1 
%             for trial=1:1:2
                for channel=1:1:16
                wtc(i,1)=id;
                wtc(i,2)=condition;
%                 wtc(i,3)=trial;
                wtc(i,3)=channel;
                wtc(i,4)=nanmean([coherence(channel,1,condition),...
                coherence(channel,2,condition)]);
                i=i+1;
                end
%             end
      elseif condition==2
%           for trial=1:1:2
                for channel=1:1:16
                wtc(i,1)=id;
                wtc(i,2)=condition;
%                 wtc(i,3)=trial;
                wtc(i,3)=channel;
                wtc(i,4)=nanmean([coherence(channel,1,condition),...
                coherence(channel,2,condition)]);
                i=i+1;
                end
%           end
      elseif condition==3
%           for trial=1:1:3
                for channel=1:1:16
                wtc(i,1)=id;
                wtc(i,2)=condition;
%                 wtc(i,3)=trial;
                wtc(i,3)=channel;
                wtc(i,4)=nanmean([coherence(channel,1,condition),...
                coherence(channel,2,condition),...
                coherence(channel,3,condition)]);
                i=i+1;
                end
%           end
       elseif condition==4
%           for trial=1
                for channel=1:1:16
                wtc(i,1)=id;
                wtc(i,2)=condition;
%                 wtc(i,3)=trial;
                wtc(i,3)=channel;
                wtc(i,4)=coherence(channel,1,condition);
                i=i+1;
                end
%           end
          end
  end
end
      
dlmwrite('P:\projects\DCARE\DCARE\data_frames\wtc_df_10_50.csv',wtc)
 %% extract RPA data
clear all 
clc

srcPath = 'P:\projects\DCARE\DCARE\MATLAB\procData\rpaData\';

%% welche Probanden sind in dem Ordner?
sourceList    = dir([srcPath, 'rpa_*']);
sourceList    = struct2cell(sourceList);
sourceList    = sourceList(1,:);
% wie viele Probanden gibt es?
numOfSources  = length(sourceList);
numOfPart       = zeros(1, numOfSources);

  for i=1:1:numOfSources
    numOfPart(i)  = sscanf(sourceList{i}, ['rpa_%d']);
  end
%%
i=1;
for id = numOfPart
    
% -------------------------------------------------------------------------
% Wo sind die Daten im jeweiligen Ordner? 
% -------------------------------------------------------------------------

  rpaFile = strcat(srcPath, sprintf(['rpa_%02d/'], id));
  load(rpaFile)

  for condition=1:1:4
        
      if condition==1 
%             for trial=1:1:2
                for channel=1:1:16
                rpa(i,1)=id;
                rpa(i,2)=condition;
%                 rpa(i,3)=trial;
                rpa(i,3)=channel;
                rpa(i,4)=nanmean([coherence(channel,1,condition),...
                coherence(channel,2,condition)]);
                i=i+1;
                end
%             end
      elseif condition==2
%           for trial=1:1:2
                for channel=1:1:16
                rpa(i,1)=id;
                rpa(i,2)=condition;
%                 rpa(i,3)=trial;
                rpa(i,3)=channel;
                rpa(i,4)=nanmean([coherence(channel,1,condition),...
                coherence(channel,2,condition)]);
                i=i+1;
                end
%           end
      elseif condition==3
%           for trial=1:1:3
                for channel=1:1:16
                rpa(i,1)=id;
                rpa(i,2)=condition;
%                 rpa(i,3)=trial;
                rpa(i,3)=channel;
                rpa(i,4)=nanmean([coherence(channel,1,condition),...
                coherence(channel,2,condition),...
                coherence(channel,3,condition)]);
                i=i+1;
                end
%           end
       elseif condition==4
%           for trial=1
                for channel=1:1:16
                rpa(i,1)=id;
                rpa(i,2)=condition;
%                 rpa(i,3)=trial;
                rpa(i,3)=channel;
                try
                rpa(i,4)=coherence(channel,1,condition);
                catch 
                rpa(i,4)=NaN;
                end
                i=i+1;
                end
%           end
          end
  end
end
 

   dlmwrite('P:\projects\DCARE\DCARE\data_frames\rpa_df_10_50.csv',rpa)
