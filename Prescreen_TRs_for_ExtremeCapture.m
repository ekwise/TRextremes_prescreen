% Prescreen for extreme capture in tree-ring records

% Example data and code to accompany this manuscript:
% Wise, Erika K. and Matthew P. Dannenberg (2019). Climate Factors Leading  
% to Asymmetric Extreme Capture in the Tree-Ring Record, Geophysical  
% Research Letters. (to be updated post-publication)

%% load example data
% this is a structure with eight ITRDB sites from AZ (site name, species,
% latitude, longitude, elevation, residual chronology, residual chronology
% years, the climate data we will test against (in this case, water-year
% precipitation from the nearest CRU TS4.01 gridpoint) and the year vector
% associated with the climate data.
load('InputData.mat')

% Example data credit:
% ITRDB data from NOAA NCEI
% (https://www.ncdc.noaa.gov/data-access/paleoclimatology-data/datasets)
% CRU data from British Atmospheric Data Centre (http://badc.nerc.ac.uk/)

% After running this code, you will have additional fields in the InputData
% structure that report the correlation with climate, the percentage of low 
% and high extreme years captured, and whether that represents significant 
% extreme capture at your p-level.  You will also have a new structure, 
% ExtremesCaptured, that contains only chronologies significantly capturing 
% both low and high extremes in the climate data. 

%% Set the extreme value threshold (EVT) and significance level:
EVT=0.2; %0.2 used in Wise & Dannenberg 2019.
sig=0.05;

%% Correlation
% Determine whether each chronology is positively or negative correlated 
% with the climate data.

for i=1:length(InputData)
    % match years of overlap:
    [~, ia, ib] = intersect(InputData(i).YEAR, InputData(i).ClimateYears); 
    % correlate:
    [InputData(i).TR_Clim_R, InputData(i).TR_Clim_p] = corr(InputData(i).RES(ia), InputData(i).ClimateData(ib), 'rows','pairwise');
end
clear i ia ib

%% Calculate extreme response and significance

for i = 1:length(InputData)
     % match years of overlap:
    [~, i1, i2] = intersect(InputData(i).YEAR, InputData(i).ClimateYears); 
     % extract the tree-ring and climate data in those years of overlap: 
    treering=InputData(i).RES(i1);
    climate=InputData(i).ClimateData(i2);
    % determine the upper and lower bounds defining extremes in TR and
    % climate data:
Y1 = quantile(treering, [EVT 1-EVT]);
Y2 = quantile(climate, [EVT 1-EVT]);

% determine which years fall aboce or below those thresholds in both the
% tree ring and climate data.  If there was a negative correlation between
% tree growth and climate, we are looking for opposite extremes (when one
% is high, the other low); if positively correlated we seek the extremes in
% the same direction:
if InputData(i).TR_Clim_R <0 % in the case of negative correlation with climate
lowYears = intersect(find(treering>=Y1(2)), find(climate<=Y2(1)));
highYears = intersect(find(treering<=Y1(1)), find(climate>=Y2(2)));
else % for all other cases (positively correlation between tree rings and climate)
    lowYears = intersect(find(treering<=Y1(1)), find(climate<=Y2(1)));
    highYears = intersect(find(treering>=Y1(2)), find(climate>=Y2(2)));
end

% determine how many extremes match versus total years evaluated:
lowResponse = length(lowYears) / length(find(climate<=Y2(1)));
highResponse = length(highYears) / length(find(climate>=Y2(2)));

% test for significance against a binomial distribution:
pLow = 2*binopdf(length(lowYears), length(find(climate<=Y2(1))), EVT);
pHigh = 2*binopdf(length(highYears), length(find(climate>=Y2(2))), EVT);

% determine the percent of low and high extremes captured and the
% associated p values:
extremePct = [lowResponse highResponse];
p = [pLow pHigh];

% add the above to your data structure:
InputData(i).Extremes.LowPct=extremePct(1);
InputData(i).Extremes.HighPct=extremePct(2);
InputData(i).Extremes.pLow=p(1);
InputData(i).Extremes.pHigh=p(2);
end

%% Which chronologies significantly capture both high and low extremes?
% Make a new structure with just those chronologies

% retains chronologies that are significantly capturing both high and low 
% extremes, based on the significance level you set above:
j=1;
for i=1:length(InputData)
    if InputData(i).Extremes.pLow <=sig && InputData(i).Extremes.pHigh <=sig
      ExtremesCaptured(j)=InputData(i);
      j=j+1;
    end
end

%% remove unneeded variables 
clearvars -except InputData ExtremesCaptured 

        

