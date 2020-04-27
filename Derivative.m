close all;
clear all;
clc;

set(0,'defaultaxesfontname','cmu serif')
set(0,'DefaultAxesFontSize',20)


%% Import data from text file.
[filename, folderpath] = uigetfile('*.txt');
cd(folderpath)
delimiter = '\t';
startRow = 2;

formatSpec = '%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

SteadyStateGrowthTimes = dataArray{:, 1};
DeconvolvedDividedSignal = dataArray{:, 2};
Length = dataArray{:, 3};
GrowthRate = dataArray{:, 4};

clearvars filename delimiter startRow formatSpec fileID dataArray ans;



for i=1:max(size(Length))
    if i == 1
        pointLocation = 'First';
        derivout(i) = threePointDerivative(SteadyStateGrowthTimes(i), Length(i), SteadyStateGrowthTimes(i+1), Length(i+1), SteadyStateGrowthTimes(i+2), Length(i+2), pointLocation);
        
    elseif i == max(size(Length))
        pointLocation = 'Last';
        derivout(i) = threePointDerivative(SteadyStateGrowthTimes(i-2), Length(i-2), SteadyStateGrowthTimes(i-1), Length(i-1), SteadyStateGrowthTimes(i), Length(i), pointLocation);
    
    else
        pointLocation = 'Interior';
        derivout(i) = threePointDerivative(SteadyStateGrowthTimes(i-1), Length(i-1), SteadyStateGrowthTimes(i), Length(i), SteadyStateGrowthTimes(i+1), Length(i+1), pointLocation);
    end
end

figure();
plot(SteadyStateGrowthTimes, derivout, 'ko-')
%title('Something','Interpreter','latex')
xlabel('$$t$$ (s)','interpreter','latex')
ylabel('Growth Rate (arb units)','interpreter','latex')
axis tight
drawnow

