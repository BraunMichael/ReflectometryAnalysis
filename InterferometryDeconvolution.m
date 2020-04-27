%% Interferometry Deconvolution
%Michael Braun

clear; close all; clc; format shortG; warning('off','all')
set(0,'defaultaxesfontname','arial')
set(0,'DefaultAxesFontSize',24)

%Save new text files? (0 no, 1 yes)
Save_file=1;

%Start of anneal/seeding time in original seconds
AnnealStartTime=45;

%GeH4 on time, only records in file, doesn't split
GermaneTime=285;

%Start of cooling time in original seconds
CoolingStartTime=644;

%Start of steady state growth time in original seconds
GrowthStartTime=705;

%End time in original seconds
EndTime=3405;

%Time check statements
if AnnealStartTime>CoolingStartTime
    fprintf('Annealing start time must be before cooling start time! Fix and try again\n')
    return
end
if AnnealStartTime>GrowthStartTime
    fprintf('Annealing start time must be before growth start time! Fix and try again\n')
    return
end
if AnnealStartTime>EndTime
    fprintf('Annealing start time must be before end time! Fix and try again\n')
    return
end
if CoolingStartTime>GrowthStartTime
    fprintf('Cooling start time must be before growth start time! Fix and try again\n')
    return
end
if CoolingStartTime>EndTime
    fprintf('Cooling start time must be before end time! Fix and try again\n')
    return
end
if GrowthStartTime>EndTime
    fprintf('Growth start time must be before end time! Fix and try again\n')
    return
end

%File Headers
header1 = 'Time (s)';
header2 = 'Raw Photovoltage (V)';
header3 = 'Laser Deconvolved Photovoltage (V)';
header4 = 'Normalized Laser Deconvolved Photovoltage';
header5 = 'Exponential Normalized Deconvolved Photovoltage';
header6 = 'Temperature (°C)';
header7 = 'Setpoint (°C)';
header8 = 'Heater Power (%)';
header9 = 'Germane Start Time (relative s)';

%% Import data from text file.
addpath(genpath('E:\Michael\Stanford\Research\Data\Interferometry'))
cd 'E:\Michael\Stanford\Research\Data\Interferometry'
%addpath(genpath('C:\Spectre Working Folder\Interferometry'))
%cd 'C:\Spectre Working Folder\Interferometry'
[interferomfilename, folderpath] = uigetfile('*.txt;*.dat');
cd(folderpath)
interferomMatrix=importdata(interferomfilename);
[~, interferomfilename_only, ~]=fileparts(interferomfilename);
datamatrix=interferomMatrix.data;

if size(datamatrix,2)==6;
    Times = datamatrix(:, 1);
    SamplePhotovoltageVtemp = datamatrix(:, 2);
    SampleStandardDeviationVtemp = datamatrix(:, 3);
    ReferencePhotovoltageVtemp = datamatrix(:, 4);
    ReferenceStandardDeviationVtemp = datamatrix(:, 5);
    Numberofmeasurementstemp = datamatrix(:, 6);
elseif size(datamatrix,2)==10;
    Times = datamatrix(:, 1);
    Temperaturetemp = datamatrix(:, 2);
    Setpointtemp = datamatrix(:, 3);
    HeaterPowertemp = datamatrix(:, 4);
    StepTimes = datamatrix(:, 5);
    SamplePhotovoltageVtemp = datamatrix(:, 6);
    SampleStandardDeviationVtemp = datamatrix(:, 7);
    ReferencePhotovoltageVtemp = datamatrix(:, 8);
    ReferenceStandardDeviationVtemp = datamatrix(:, 9);
    Numberofmeasurementstemp = datamatrix(:, 10);
else
    fprintf('Check the number of data columns and modify the code accordingly\n')
    return
end


%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%% Full Graph

CutTimes=Times(Times>AnnealStartTime & Times<EndTime);
SamplePhotovoltageV = SamplePhotovoltageVtemp(Times>AnnealStartTime & Times<EndTime);
SampleStandardDeviationV = SampleStandardDeviationVtemp(Times>AnnealStartTime & Times<EndTime);
ReferencePhotovoltageV = ReferencePhotovoltageVtemp(Times>AnnealStartTime & Times<EndTime);
ReferenceStandardDeviationV = ReferenceStandardDeviationVtemp(Times>AnnealStartTime & Times<EndTime);
Numberofmeasurements = Numberofmeasurementstemp(Times>AnnealStartTime & Times<EndTime);
if size(datamatrix,2)==10;
    Temperature = Temperaturetemp(Times>AnnealStartTime & Times<EndTime);
    Setpoint = Setpointtemp(Times>AnnealStartTime & Times<EndTime);
    HeaterPower = HeaterPowertemp(Times>AnnealStartTime & Times<EndTime);
    CutStepTimes = StepTimes(Times>AnnealStartTime & Times<EndTime);
end

t=CutTimes-CutTimes(1);
FullGermaneTime=(GermaneTime-CutTimes(1)).*ones(size(t));
R=ReferencePhotovoltageV;
S=SamplePhotovoltageV;

deconv=(S./R).*mean(R);
figure()
plot(t,deconv,'k')
ylabel('Signal (arb units)','Interpreter','latex')
title('Interferometry Laser Deconvolution','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex')
axis([0 inf -inf inf])

norm=deconv./max(deconv);
figure()
plot(t,norm,'k')
ylabel('Signal (arb units)','Interpreter','latex')
title('Normalized Interferometry Laser Deconvolution','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex')
axis([0 inf -inf inf])


figure()
yyaxis left
plot(t,S,'k')
ylabel('Sample Signal (Volts)','Interpreter','latex')
yyaxis right
plot(t,R,'b')
title('Raw Interferometry Data','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex')
ylabel('Reference Signal (Volts)','Interpreter','latex')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
legend('Sample','Reference','Location','Best')
axis([0 inf -inf inf])


%% Annealing and Seeding Section
AnnealCutTimes=Times(Times>AnnealStartTime & Times<CoolingStartTime);
AnnealSamplePhotovoltageV = SamplePhotovoltageVtemp(Times>AnnealStartTime & Times<CoolingStartTime);
AnnealSampleStandardDeviationV = SampleStandardDeviationVtemp(Times>AnnealStartTime & Times<CoolingStartTime);
AnnealReferencePhotovoltageV = ReferencePhotovoltageVtemp(Times>AnnealStartTime & Times<CoolingStartTime);
AnnealReferenceStandardDeviationV = ReferenceStandardDeviationVtemp(Times>AnnealStartTime & Times<CoolingStartTime);
AnnealNumberofmeasurements = Numberofmeasurementstemp(Times>AnnealStartTime & Times<CoolingStartTime);
if size(datamatrix,2)==10;
    AnnealTemperature = Temperaturetemp(Times>AnnealStartTime & Times<CoolingStartTime);
    AnnealSetpoint = Setpointtemp(Times>AnnealStartTime & Times<CoolingStartTime);
    AnnealHeaterPower = HeaterPowertemp(Times>AnnealStartTime & Times<CoolingStartTime);
    AnnealCutStepTimes = StepTimes(Times>AnnealStartTime & Times<CoolingStartTime);
end


Annealt=AnnealCutTimes-AnnealCutTimes(1);
AnnealGermaneTime=(GermaneTime-AnnealCutTimes(1)).*ones(size(Annealt));
AnnealR=AnnealReferencePhotovoltageV;
AnnealS=AnnealSamplePhotovoltageV;

Annealdeconv=(AnnealS./AnnealR).*mean(AnnealR);
figure()
plot(Annealt,Annealdeconv,'k')
ylabel('Signal (arb units)','Interpreter','latex')
title('Annealing and Seeding Interferometry Laser Deconvolution','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex')
axis([0 inf -inf inf])

Annealnorm=Annealdeconv./max(Annealdeconv);
figure()
plot(Annealt,Annealnorm,'k')
ylabel('Signal (arb units)','Interpreter','latex')
title('Annealing and Seeding Normalized Interferometry Laser Deconvolution','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex')
axis([0 inf -inf inf])


figure()
yyaxis left
plot(Annealt,AnnealS,'k')
ylabel('Sample Signal (Volts)','Interpreter','latex')
yyaxis right
plot(Annealt,AnnealR,'b')
title('Annealing and Seeding Raw Interferometry Data','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex')
ylabel('Reference Signal (Volts)','Interpreter','latex')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
legend('Sample','Reference','Location','Best')
axis([0 inf -inf inf])


%% Cooling Section
CoolingCutTimes=Times(Times>CoolingStartTime & Times<GrowthStartTime);
CoolingSamplePhotovoltageV = SamplePhotovoltageVtemp(Times>CoolingStartTime & Times<GrowthStartTime);
CoolingSampleStandardDeviationV = SampleStandardDeviationVtemp(Times>CoolingStartTime & Times<GrowthStartTime);
CoolingReferencePhotovoltageV = ReferencePhotovoltageVtemp(Times>CoolingStartTime & Times<GrowthStartTime);
CoolingReferenceStandardDeviationV = ReferenceStandardDeviationVtemp(Times>CoolingStartTime & Times<GrowthStartTime);
CoolingNumberofmeasurements = Numberofmeasurementstemp(Times>CoolingStartTime & Times<GrowthStartTime);
if size(datamatrix,2)==10;
    CoolingTemperature = Temperaturetemp(Times>CoolingStartTime & Times<GrowthStartTime);
    CoolingSetpoint = Setpointtemp(Times>CoolingStartTime & Times<GrowthStartTime);
    CoolingHeaterPower = HeaterPowertemp(Times>CoolingStartTime & Times<GrowthStartTime);
    CoolingCutStepTimes = StepTimes(Times>CoolingStartTime & Times<GrowthStartTime);
end

Coolingt=CoolingCutTimes-CoolingCutTimes(1);
CoolingR=CoolingReferencePhotovoltageV;
CoolingS=CoolingSamplePhotovoltageV;

Coolingdeconv=(CoolingS./CoolingR).*mean(CoolingR);
figure()
plot(Coolingt,Coolingdeconv,'k')
ylabel('Signal (arb units)','Interpreter','latex')
title('Cooling Interferometry Laser Deconvolution','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex')
axis([0 inf -inf inf])

Coolingnorm=Coolingdeconv./max(Coolingdeconv);
figure()
plot(Coolingt,Coolingnorm,'k')
ylabel('Signal (arb units)','Interpreter','latex')
title('Cooling Normalized Interferometry Laser Deconvolution','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex')
axis([0 inf -inf inf])


%% Growth Section
GrowthCutTimes=Times(Times>GrowthStartTime & Times<EndTime);
GrowthSamplePhotovoltageV = SamplePhotovoltageVtemp(Times>GrowthStartTime & Times<EndTime);
GrowthSampleStandardDeviationV = SampleStandardDeviationVtemp(Times>GrowthStartTime & Times<EndTime);
GrowthReferencePhotovoltageV = ReferencePhotovoltageVtemp(Times>GrowthStartTime & Times<EndTime);
GrowthReferenceStandardDeviationV = ReferenceStandardDeviationVtemp(Times>GrowthStartTime & Times<EndTime);
GrowthNumberofmeasurements = Numberofmeasurementstemp(Times>GrowthStartTime & Times<EndTime);
if size(datamatrix,2)==10;
    GrowthTemperature = Temperaturetemp(Times>GrowthStartTime & Times<EndTime);
    GrowthSetpoint = Setpointtemp(Times>GrowthStartTime & Times<EndTime);
    GrowthHeaterPower = HeaterPowertemp(Times>GrowthStartTime & Times<EndTime);
    GrowthCutStepTimes = StepTimes(Times>GrowthStartTime & Times<EndTime);
end

Growtht=GrowthCutTimes-GrowthCutTimes(1);
GrowthR=GrowthReferencePhotovoltageV;
GrowthS=GrowthSamplePhotovoltageV;

Growthdeconv=(GrowthS./GrowthR).*mean(GrowthR);
figure()
plot(Growtht,Growthdeconv,'k')
ylabel('Signal (arb units)','Interpreter','latex')
title('Growth Interferometry Laser Deconvolution','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex')
axis([0 inf -inf inf])

Growthnorm=Growthdeconv./max(Growthdeconv);
figure()
plot(Growtht,Growthnorm,'k')
ylabel('Signal (arb units)','Interpreter','latex')
title('Growth Normalized Interferometry Laser Deconvolution','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex')
axis([0 inf -inf inf])


Growthf=fit(Growtht,Growthdeconv,'exp2')
Growthf_coeff=coeffvalues(Growthf);
a=Growthf_coeff(1);
b=Growthf_coeff(2);
c=Growthf_coeff(3);
d=Growthf_coeff(4);

Growthexpfit = @(x) a.*exp(b.*x)+c.*exp(d.*x);
Growthf_values=Growthexpfit(Growtht);

Growthdiv=Growthdeconv./Growthf_values;

figure()
plot(Growtht,Growthdiv,'k')
title('Growth Interferometry Deconvolved Exponential Division','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex')
ylabel('Signal (Volts)','Interpreter','latex')
axis([0 inf -inf inf])

figure()
yyaxis left
plot(Growtht,GrowthS,'k',Growtht,Growthf_values,'r')
ylabel('Sample Signal (Volts)','Interpreter','latex')
yyaxis right
plot(Growtht,GrowthR,'b')
title('Growth Raw Interferometry Data','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex')
ylabel('Reference Signal (Volts)','Interpreter','latex')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
legend('Sample','Exponential fit','Reference','Location','Best')
axis([0 inf -inf inf])

%Save file
if Save_file==0
else
    fullfilename=strcat(interferomfilename_only,'_full.txt');
    Annealingfilename=strcat(interferomfilename_only,'_annealing.txt');
    Coolingfilename=strcat(interferomfilename_only,'_cooling.txt');
    Growthfilename=strcat(interferomfilename_only,'_growth.txt');
    
    delete(fullfilename)
    delete(Annealingfilename)
    delete(Coolingfilename)
    delete(Growthfilename)
    
    if size(datamatrix,2)==6;
        fid=fopen(fullfilename,'wt');
        fprintf(fid, [ header1 '\t' header2 '\t' header3 '\t' header4 '\t' header9 '\r\n']);
        fprintf(fid, '%f\t%f\t%f\t%f\t%f\r\n', [t S deconv norm FullGermaneTime]');
        fclose(fid);true
        
        fid=fopen(Annealingfilename,'wt');
        fprintf(fid, [ header1 '\t' header2 '\t' header3 '\t' header4 '\t' header9 '\r\n']);
        fprintf(fid, '%f\t%f\t%f\t%f\t%f\r\n', [Annealt AnnealS Annealdeconv Annealnorm AnnealGermaneTime]');
        fclose(fid);true
        
        fid=fopen(Coolingfilename,'wt');
        fprintf(fid, [ header1 '\t' header2 '\t' header3 '\t' header4 '\r\n']);
        fprintf(fid, '%f\t%f\t%f\t%f\r\n', [Coolingt CoolingS Coolingdeconv Coolingnorm]');
        fclose(fid);true
        
        fid=fopen(Growthfilename,'wt');
        fprintf(fid, [ header1 '\t' header2 '\t' header3 '\t' header4 '\t' header5 '\r\n']);
        fprintf(fid, '%f\t%f\t%f\t%f\t%f\r\n', [Growtht GrowthS Growthdeconv Growthnorm Growthdiv]');
        fclose(fid);true
        
    elseif size(datamatrix,2)==10;
        fid=fopen(fullfilename,'wt');
        fprintf(fid, [ header1 '\t' header2 '\t' header3 '\t' header4 '\t' header6 '\t' header7 '\t' header8 '\t' header9 '\r\n']);
        fprintf(fid, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\r\n', [t S deconv norm Temperature Setpoint HeaterPower FullGermaneTime]');
        fclose(fid);true
        
        fid=fopen(Annealingfilename,'wt');
        fprintf(fid, [ header1 '\t' header2 '\t' header3 '\t' header4 '\t' header6 '\t' header7 '\t' header8 '\t' header9 '\r\n']);
        fprintf(fid, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\r\n', [Annealt AnnealS Annealdeconv Annealnorm AnnealTemperature AnnealSetpoint AnnealHeaterPower AnnealGermaneTime]');
        fclose(fid);true
        
        fid=fopen(Coolingfilename,'wt');
        fprintf(fid, [ header1 '\t' header2 '\t' header3 '\t' header4 '\t' header6 '\t' header7 '\t' header8 '\r\n']);
        fprintf(fid, '%f\t%f\t%f\t%f\t%f\t%f\t%f\r\n', [Coolingt CoolingS Coolingdeconv Coolingnorm CoolingTemperature CoolingSetpoint CoolingHeaterPower]');
        fclose(fid);true
        
        fid=fopen(Growthfilename,'wt');
        fprintf(fid, [ header1 '\t' header2 '\t' header3 '\t' header4 '\t' header5 '\t' header6 '\t' header7 '\t' header8 '\r\n']);
        fprintf(fid, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\r\n', [Growtht GrowthS Growthdeconv Growthnorm Growthdiv GrowthTemperature GrowthSetpoint GrowthHeaterPower]');
        fclose(fid);true
        
    else
        fprintf('Remember to modify the file output code as well for different number/order of data columns!\n')
        return
    end
    
end