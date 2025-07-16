function DCARE_6

for i=[1:9 11:31 33:48 50:52 54:68] 
    load(sprintf(['P:/projects/DCARE/DCARE/MATLAB/procData/wtcData/DCARE_%02d_wtc.mat'],i));
  for cond=1:3   
  % wtc
  CARE_wtc((length(wtc_all)*cond-15+((i-1)*48)):...
      (length(wtc_all)*cond+((i-1)*48)),1)  =  wtc_all(1:16,cond);
   % dyad ID
  CARE_wtc((length(wtc_all)*cond-15+((i-1)*48)):...
      (length(wtc_all)*cond+((i-1)*48)),2)  =  [i i i i i i i i i i i i i i i i];
  % cond
  CARE_wtc((length(wtc_all)*cond-15+((i-1)*48)):...
      (length(wtc_all)*cond+((i-1)*48)),3)  = [cond cond cond cond cond cond cond cond cond cond cond cond cond cond cond cond];
  % ch
  CARE_wtc((length(wtc_all)*cond-15+((i-1)*48)):...
      (length(wtc_all)*cond+((i-1)*48)),5)  = [1:16];
  % roi
  CARE_wtc((length(wtc_all)*cond-15+((i-1)*48)):...
      (length(wtc_all)*cond+((i-1)*48)),6)  = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4];
    end
end

dlmwrite('P:/projects/DCARE/DCARE/MATLAB/procData/Export/DCARE_wtc_hbo.csv',CARE_wtc) 

for i=[1:9 11:31 33:48 50:52 54:68] 
    load(sprintf(['P:/projects/DCARE/DCARE/MATLAB/procData/wtcData/DCARE_%02d_wtc_hbr.mat'],i));
  for cond=1:3   
  % wtc
  CARE_wtc_hbr((length(wtc_all)*cond-15+((i-1)*48)):...
      (length(wtc_all)*cond+((i-1)*48)),1)  =  wtc_all(1:16,cond);
   % dyad ID
  CARE_wtc_hbr((length(wtc_all)*cond-15+((i-1)*48)):...
      (length(wtc_all)*cond+((i-1)*48)),2)  =  [i i i i i i i i i i i i i i i i];
  % cond
  CARE_wtc_hbr((length(wtc_all)*cond-15+((i-1)*48)):...
      (length(wtc_all)*cond+((i-1)*48)),3)  = [cond cond cond cond cond cond cond cond cond cond cond cond cond cond cond cond];
  % ch
  CARE_wtc_hbr((length(wtc_all)*cond-15+((i-1)*48)):...
      (length(wtc_all)*cond+((i-1)*48)),5)  = [1:16];
  % roi
  CARE_wtc_hbr((length(wtc_all)*cond-15+((i-1)*48)):...
      (length(wtc_all)*cond+((i-1)*48)),6)  = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4];
    end
end

dlmwrite('P:/projects/DCARE/DCARE/MATLAB/procData/Export/DCARE_wtc_hbr.csv',CARE_wtc_hbr) 


end