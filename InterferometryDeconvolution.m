
%% Interferometry Deconvolution
%Michael Braun

clear; close all; format shortG; warning('off','all')
set(0,'defaultaxesfontname','arial')
set(0,'DefaultAxesFontSize',24)

%Save new text files? (0 no, 1 yes, true false also works)
Save_file = true;

%Load from saved .mat file? (0 no, 1 yes)
Load_mat = true;

%Show annealing/cooling graphs? (0 for no, 1 for yes)
AnnealCoolGraphShow = false;

%Show intermediate graphs? (0 for no, 1 for yes)
IntermediateShow = false;

%Update smoothing factors? (0 for no, 1 for yes)
Update_smooth = false;

if Update_smooth
    %Growth Loess Smoothing factor
    temp_Growth_loess = 0.03;
    
    %Derivative Loess Smoothing factor
    temp_Derivative_loess = 0.03;
end

if ~Load_mat
    %Normal Growth splits? (0 for a constant temperature/step file, 1 for
    %normal growth)
    %Basically only does laser deconvolution and normalization of file
    %Only uses Anneal StartTime and EndTime if 0
    GrowthSplit = 1;
    
    %Growth Loess Smoothing factor
    Growth_loess = 0.07;
    
    %Derivative Loess Smoothing factor
    Derivative_loess = 0.07;
    
    %Start of anneal/seeding time in original seconds
    AnnealStartTime = 70.8;
    
    %GeH4 on time, only records in file, doesn't split
    GermaneTime = 310.8;
    
    %Start of cooling time in original seconds
    CoolingStartTime = 551.2;
    
    %Start of steady state growth time in original seconds
    GrowthStartTime = 611.5;
    
    %End time in original seconds
    EndTime = 2411;
    
    %Time check statements
    if AnnealStartTime>EndTime
        fprintf('Annealing start time must be before end time! Fix and try again\n')
        return
    end
    if GrowthSplit==1;
        if AnnealStartTime>CoolingStartTime
            fprintf('Annealing start time must be before cooling start time! Fix and try again\n')
            return
        end
        if AnnealStartTime>GrowthStartTime
            fprintf('Annealing start time must be before growth start time! Fix and try again\n')
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
    end
end

%File Headers
header1 = 'Time (s)';
header2 = 'Raw Photovoltage (V)';
headerlaser = 'Laser Photovoltage (V)';
header3 = 'Laser Deconvolved Photovoltage (V)';
header4 = 'Normalized Laser Deconvolved Photovoltage';
header5 = 'Exponential Normalized Deconvolved Photovoltage';
header6 = 'Temperature (°C)';
header7 = 'Setpoint (°C)';
header8 = 'Heater Power (%%)';
header9 = 'Germane Start Time (relative s)';


%% Import data from text file.
addpath(genpath('E:\Michael\Stanford\Research\Data\Reflectometry'))
cd 'E:\Michael\Stanford\Research\Data\Reflectometry'
%addpath(genpath('C:\Spectre Working Folder\Reflectometry'))
%cd 'C:\Spectre Working Folder\Reflectometry'
[interferomfilename, folderpath] = uigetfile('*.txt;*.dat');
cd(folderpath)
[~, interferomfilename_only, ~]=fileparts(interferomfilename);
clc;
fprintf('Working on %s\n', interferomfilename_only);
matfilename = strcat(interferomfilename_only, '.mat');
if Load_mat
    if ~exist(matfilename, 'file') == 2
        errorMessage = sprintf('Error: The following folder does not exist:\n%s', folderpath);
        uiwait(warndlg(errorMessage));
        return;
    end
    load(matfilename)
end


if Update_smooth
    %Growth Loess Smoothing factor
    Growth_loess = temp_Growth_loess;
    
    %Derivative Loess Smoothing factor
    Derivative_loess = temp_Derivative_loess;
end

[Times,Temperaturetemp,Setpointtemp,HeaterPowertemp,StepTimes,SamplePhotovoltageVtemp,SampleStandardDeviationVtemp,ReferencePhotovoltageVtemp,ReferenceStandardDeviationVtemp,~,Numberofmeasurementstemp,Timestamp] = ImportFunction(interferomfilename);


%% Full Graph

CutTimes=Times(Times>AnnealStartTime & Times<EndTime);
SamplePhotovoltageV = SamplePhotovoltageVtemp(Times>AnnealStartTime & Times<EndTime);
SampleStandardDeviationV = SampleStandardDeviationVtemp(Times>AnnealStartTime & Times<EndTime);
ReferencePhotovoltageV = ReferencePhotovoltageVtemp(Times>AnnealStartTime & Times<EndTime);
ReferenceStandardDeviationV = ReferenceStandardDeviationVtemp(Times>AnnealStartTime & Times<EndTime);
Numberofmeasurements = Numberofmeasurementstemp(Times>AnnealStartTime & Times<EndTime);
Temperature = Temperaturetemp(Times>AnnealStartTime & Times<EndTime);
Setpoint = Setpointtemp(Times>AnnealStartTime & Times<EndTime);
HeaterPower = HeaterPowertemp(Times>AnnealStartTime & Times<EndTime);
CutStepTimes = StepTimes(Times>AnnealStartTime & Times<EndTime);


t=CutTimes-CutTimes(1);
FullGermaneTime=(GermaneTime-CutTimes(1)).*ones(size(t));
R=ReferencePhotovoltageV;
S=SamplePhotovoltageV;

deconv=(S./R).*mean(R);
if IntermediateShow == 1
    figure()
    plot(t,deconv,'k')
    ylabel('Signal (arb units)','Interpreter','latex')
    title('Interferometry Laser Deconvolution','Interpreter','latex')
    xlabel('Time (s)','Interpreter','latex')
    axis([0 inf -inf inf])
end

norm=deconv./max(deconv);
if IntermediateShow == 1
    figure()
    plot(t,norm,'k')
    ylabel('Signal (arb units)','Interpreter','latex')
    title('Normalized Interferometry Laser Deconvolution','Interpreter','latex')
    xlabel('Time (s)','Interpreter','latex')
    axis([0 inf -inf inf])
end

if IntermediateShow == 1
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
end

if GrowthSplit==1;
    %% Annealing and Seeding Section
    AnnealCutTimes=Times(Times>AnnealStartTime & Times<CoolingStartTime);
    AnnealSamplePhotovoltageV = SamplePhotovoltageVtemp(Times>AnnealStartTime & Times<CoolingStartTime);
    AnnealSampleStandardDeviationV = SampleStandardDeviationVtemp(Times>AnnealStartTime & Times<CoolingStartTime);
    AnnealReferencePhotovoltageV = ReferencePhotovoltageVtemp(Times>AnnealStartTime & Times<CoolingStartTime);
    AnnealReferenceStandardDeviationV = ReferenceStandardDeviationVtemp(Times>AnnealStartTime & Times<CoolingStartTime);
    AnnealNumberofmeasurements = Numberofmeasurementstemp(Times>AnnealStartTime & Times<CoolingStartTime);
    AnnealTemperature = Temperaturetemp(Times>AnnealStartTime & Times<CoolingStartTime);
    AnnealSetpoint = Setpointtemp(Times>AnnealStartTime & Times<CoolingStartTime);
    AnnealHeaterPower = HeaterPowertemp(Times>AnnealStartTime & Times<CoolingStartTime);
    AnnealCutStepTimes = StepTimes(Times>AnnealStartTime & Times<CoolingStartTime);
    
    
    Annealt=AnnealCutTimes-AnnealCutTimes(1);
    AnnealGermaneTime=(GermaneTime-AnnealCutTimes(1)).*ones(size(Annealt));
    AnnealR=AnnealReferencePhotovoltageV;
    AnnealS=AnnealSamplePhotovoltageV;
    
    Annealdeconv=(AnnealS./AnnealR).*mean(AnnealR);
    if AnnealCoolGraphShow && IntermediateShow
        figure()
        plot(Annealt,Annealdeconv,'k')
        ylabel('Signal (arb units)','Interpreter','latex')
        title('Annealing and Seeding Interferometry Laser Deconvolution','Interpreter','latex')
        xlabel('Time (s)','Interpreter','latex')
        axis([0 inf -inf inf])
    end
    
    Annealnorm=Annealdeconv./max(Annealdeconv);
    if AnnealCoolGraphShow && IntermediateShow
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
    end
    
    %% Cooling Section
    CoolingCutTimes=Times(Times>CoolingStartTime & Times<GrowthStartTime);
    CoolingSamplePhotovoltageV = SamplePhotovoltageVtemp(Times>CoolingStartTime & Times<GrowthStartTime);
    CoolingSampleStandardDeviationV = SampleStandardDeviationVtemp(Times>CoolingStartTime & Times<GrowthStartTime);
    CoolingReferencePhotovoltageV = ReferencePhotovoltageVtemp(Times>CoolingStartTime & Times<GrowthStartTime);
    CoolingReferenceStandardDeviationV = ReferenceStandardDeviationVtemp(Times>CoolingStartTime & Times<GrowthStartTime);
    CoolingNumberofmeasurements = Numberofmeasurementstemp(Times>CoolingStartTime & Times<GrowthStartTime);
    CoolingTemperature = Temperaturetemp(Times>CoolingStartTime & Times<GrowthStartTime);
    CoolingSetpoint = Setpointtemp(Times>CoolingStartTime & Times<GrowthStartTime);
    CoolingHeaterPower = HeaterPowertemp(Times>CoolingStartTime & Times<GrowthStartTime);
    CoolingCutStepTimes = StepTimes(Times>CoolingStartTime & Times<GrowthStartTime);
    
    Coolingt=CoolingCutTimes-CoolingCutTimes(1);
    CoolingR=CoolingReferencePhotovoltageV;
    CoolingS=CoolingSamplePhotovoltageV;
    
    Coolingdeconv=(CoolingS./CoolingR).*mean(CoolingR);
    if AnnealCoolGraphShow && IntermediateShow
        figure()
        plot(Coolingt,Coolingdeconv,'k')
        ylabel('Signal (arb units)','Interpreter','latex')
        title('Cooling Interferometry Laser Deconvolution','Interpreter','latex')
        xlabel('Time (s)','Interpreter','latex')
        axis([0 inf -inf inf])
    end
    
    Coolingnorm=Coolingdeconv./max(Coolingdeconv);
    if AnnealCoolGraphShow && IntermediateShow
        figure()
        plot(Coolingt,Coolingnorm,'k')
        ylabel('Signal (arb units)','Interpreter','latex')
        title('Cooling Normalized Interferometry Laser Deconvolution','Interpreter','latex')
        xlabel('Time (s)','Interpreter','latex')
        axis([0 inf -inf inf])
    end
    
    %% Growth Section
    GrowthCutTimes=Times(Times>GrowthStartTime & Times<EndTime);
    GrowthSamplePhotovoltageV = SamplePhotovoltageVtemp(Times>GrowthStartTime & Times<EndTime);
    GrowthSampleStandardDeviationV = SampleStandardDeviationVtemp(Times>GrowthStartTime & Times<EndTime);
    GrowthReferencePhotovoltageV = ReferencePhotovoltageVtemp(Times>GrowthStartTime & Times<EndTime);
    GrowthReferenceStandardDeviationV = ReferenceStandardDeviationVtemp(Times>GrowthStartTime & Times<EndTime);
    GrowthNumberofmeasurements = Numberofmeasurementstemp(Times>GrowthStartTime & Times<EndTime);
    GrowthTemperature = Temperaturetemp(Times>GrowthStartTime & Times<EndTime);
    GrowthSetpoint = Setpointtemp(Times>GrowthStartTime & Times<EndTime);
    GrowthHeaterPower = HeaterPowertemp(Times>GrowthStartTime & Times<EndTime);
    GrowthCutStepTimes = StepTimes(Times>GrowthStartTime & Times<EndTime);
    
    
    Growtht=GrowthCutTimes-GrowthCutTimes(1);
    GrowthR=GrowthReferencePhotovoltageV;
    GrowthS=GrowthSamplePhotovoltageV;
    
    Growthdeconv=(GrowthS./GrowthR).*mean(GrowthR);
    if IntermediateShow
        figure()
        plot(Growtht,Growthdeconv,'k')
        ylabel('Signal (arb units)','Interpreter','latex')
        title('Growth Interferometry Laser Deconvolution','Interpreter','latex')
        xlabel('Time (s)','Interpreter','latex')
        axis([0 inf -inf inf])
    end
    
    Growthnorm=Growthdeconv./max(Growthdeconv);
    if IntermediateShow
        figure()
        plot(Growtht,Growthnorm,'k')
        ylabel('Signal (arb units)','Interpreter','latex')
        title('Growth Normalized Interferometry Laser Deconvolution','Interpreter','latex')
        xlabel('Time (s)','Interpreter','latex')
        axis([0 inf -inf inf])
    end
    
    Growthf=fit(Growtht,Growthdeconv,'exp2');
    Growthf_coeff=coeffvalues(Growthf);
    a=Growthf_coeff(1);
    b=Growthf_coeff(2);
    c=Growthf_coeff(3);
    d=Growthf_coeff(4);
    
    %Save .mat file
    save(matfilename, 'GrowthSplit', 'AnnealStartTime', 'GermaneTime', 'CoolingStartTime', 'GrowthStartTime', 'EndTime', 'Growth_loess', 'Derivative_loess', 'a', 'b', 'c', 'd')
    
    Growthexpfit = @(x) a.*exp(b.*x)+c.*exp(d.*x);
    Growthf_values=Growthexpfit(Growtht);
    
    Growthdiv=Growthdeconv./Growthf_values;
    
    figure()
    plot(Growtht,Growthdiv,'k')
    title('Growth Interferometry Deconvolved Exponential Division','Interpreter','latex')
    xlabel('Time (s)','Interpreter','latex')
    ylabel('Signal (Volts)','Interpreter','latex')
    axis([0 inf -inf inf])
    
    if IntermediateShow
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
    end
    
    Growth_smoothed = smooth(Growtht,Growthdiv,Growth_loess,'rloess');
    Deriv1_Growth_smoothed = nderiv_fornberg(1, Growtht, Growth_smoothed);
    Deriv2_Growth_smoothed = nderiv_fornberg(2, Growtht, Growth_smoothed);
    Deriv2_Growth_doublesmoothed = smooth(Growtht,Deriv2_Growth_smoothed,Derivative_loess,'rloess');
    [zeros_1stderiv_indices, zeros_1stderiv_xvalues] = NumericalRootsFunction(Growtht, Deriv1_Growth_smoothed);
    [zeros_2ndderiv_indices, zeros_2ndderiv_xvalues] = NumericalRootsFunction(Growtht, Deriv2_Growth_doublesmoothed);
    zeros_xindices = sort([zeros_1stderiv_indices, zeros_2ndderiv_indices]);
    zeros_xvalues_exact = sort([zeros_1stderiv_xvalues, zeros_2ndderiv_xvalues])';
    
    zeros_xvalues = Growtht(zeros_xindices);
    zeros_yvalues = Growthdiv(zeros_xindices);
    
    if IntermediateShow
        figure()
        plot(Growtht,Growth_smoothed-1,'k',Growtht,Deriv1_Growth_smoothed,'r',Growtht,Deriv2_Growth_smoothed,'b')
        title('Growth Interferometry Deconvolved Exponential Division Smoothed Derivatives','Interpreter','latex')
        xlabel('Time (s)','Interpreter','latex')
        ylabel('Signal (Volts)','Interpreter','latex')
        axis([0 inf -inf inf])
        
        figure()
        plot(Growtht,Deriv2_Growth_doublesmoothed,'k',Growtht,Deriv2_Growth_smoothed,'b')
        title('Growth Interferometry Deconvolved Exponential Division Smoothed Derivatives','Interpreter','latex')
        xlabel('Time (s)','Interpreter','latex')
        ylabel('Signal (Volts)','Interpreter','latex')
        axis([0 inf -inf inf])
    end
    
    figure()
    plot(Growtht,Growthdiv,'k',Growtht,Growth_smoothed,'r', zeros_xvalues, zeros_yvalues, 'o')
    title('Growth Interferometry Deconvolved Exponential Division Smoothed','Interpreter','latex')
    xlabel('Time (s)','Interpreter','latex')
    ylabel('Signal (Volts)','Interpreter','latex')
    axis([0 inf -inf inf])
    
    length_segments = (1:max(size(zeros_xvalues)));
    %segmentoffset = ((length_segments(2)*zeros_xvalues(1) - zeros_xvalues(2)*length_segments(1))/(zeros_xvalues(2)-zeros_xvalues(1)));
    segmentoffset = ((length_segments(2)*zeros_xvalues_exact(1) - zeros_xvalues_exact(2)*length_segments(1))/(zeros_xvalues_exact(2)-zeros_xvalues_exact(1)));
    %segmentoffset2 = length_segments(1)-zeros_xvalues(1)*((length_segments(2) - length_segments(1))/(zeros_xvalues(2)-zeros_xvalues(1)));
    
    figure()
    plot(zeros_xvalues_exact, length_segments+segmentoffset, 'o')
    title('Length vs Time','Interpreter','latex')
    xlabel('Time (s)','Interpreter','latex')
    ylabel('Length (arb units)','Interpreter','latex')
    axis([0 inf 0 inf])
    
    Growth_rate = nderiv_fornberg(1, zeros_xvalues_exact, length_segments+segmentoffset);
    figure()
    plot(zeros_xvalues_exact, Growth_rate, 'o')
    title('Growth Rate vs Time','Interpreter','latex')
    xlabel('Time (s)','Interpreter','latex')
    ylabel('Growth Rate (arb units)','Interpreter','latex')
    axis([0 inf 0 inf])
    
    %     zeros_xindices = sort([zeros_1stderiv_indices, zeros_2ndderiv_indices]);
    %     peaks_yvalues = Growthdiv(zeros_1stderiv_indices);
    %     [top, middle, bottom] = PeakIndiciesSplit(zeros_xindices, zeros_1stderiv_indices, zeros_2ndderiv_indices, peaks_yvalues);
    %     Osc_middle_fit=polyfit(Growtht(middle),Growthdiv(middle),4);
    %
    %     Osc_cleaned=Growthdiv./polyval(Osc_middle_fit, Growtht);
    %     figure()
    %     plot(Growtht,Growthdiv,'k',Growtht,Osc_cleaned,'r')
    %     title('Osciallation Cleaning','Interpreter','latex')
    %     xlabel('Time (s)','Interpreter','latex')
    %     ylabel('Signal (Volts)','Interpreter','latex')
    %     axis([0 inf -inf inf])
    
    
    
    
    
end




%Save file
if Save_file==0
else
    samplefolderpath = strcat(folderpath,interferomfilename_only,'_Processed');
    mkdir(samplefolderpath)
    cd(samplefolderpath)
    fullfilename=strcat(interferomfilename_only,'_full.txt');
    delete(fullfilename)
    
    if GrowthSplit==1;
        Annealingfilename=strcat(interferomfilename_only,'_annealing.txt');
        Coolingfilename=strcat(interferomfilename_only,'_cooling.txt');
        Growthfilename=strcat(interferomfilename_only,'_growth.txt');
        Paramfilename=strcat(interferomfilename_only,'_timeparameters.txt');
        Processedfilename=strcat(interferomfilename_only,'_processed.txt');
        
        
        delete(Annealingfilename)
        delete(Coolingfilename)
        delete(Growthfilename)
        delete(Paramfilename)
        delete(Processedfilename)
    end
    
    
    fid=fopen(fullfilename,'wt');
    fprintf(fid, [ header1 '\t' header2 '\t' headerlaser '\t' header3 '\t' header4 '\t' header6 '\t' header7 '\t' header8 '\t' header9 '\r\n']);
    fprintf(fid, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\r\n', [t S R deconv norm Temperature Setpoint HeaterPower FullGermaneTime]');
    fclose(fid);true;
    
    
    if GrowthSplit==1;
        fid=fopen(Paramfilename,'wt');
        fprintf(fid, [ 'AnnealStartTime' '\t' 'GermaneTime' '\t' 'CoolingStartTime' '\t' 'GrowthStartTime' '\t' 'EndTime' '\r\n']);
        fprintf(fid, '%6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f\r\n', [AnnealStartTime GermaneTime CoolingStartTime GrowthStartTime EndTime]');
        fclose(fid);true;
        
        fid=fopen(Processedfilename,'wt');
        fprintf(fid, [ 'Steady State Growth Time (s)' '\t' 'Deconvolved Divided Signal' '\t'  'Length (arb units)' '\t' 'Growth Rate (arb units/s)' '\t' '\r\n']);
        fprintf(fid, '%f\t%f\t%f\t%f\r\n', [zeros_xvalues_exact zeros_yvalues (length_segments+segmentoffset)' Growth_rate']');
        fclose(fid);true;
        fid=fopen(Annealingfilename,'wt');
        fprintf(fid, [ header1 '\t' header2 '\t' headerlaser '\t' header3 '\t' header4 '\t' header6 '\t' header7 '\t' header8 '\t' header9 '\r\n']);
        fprintf(fid, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\r\n', [Annealt AnnealS AnnealR Annealdeconv Annealnorm AnnealTemperature AnnealSetpoint AnnealHeaterPower AnnealGermaneTime]');
        fclose(fid);true;
        
        fid=fopen(Coolingfilename,'wt');
        fprintf(fid, [ header1 '\t' header2 '\t' headerlaser '\t' header3 '\t' header4 '\t' header6 '\t' header7 '\t' header8 '\r\n']);
        fprintf(fid, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\r\n', [Coolingt CoolingS CoolingR Coolingdeconv Coolingnorm CoolingTemperature CoolingSetpoint CoolingHeaterPower]');
        fclose(fid);true;
        
        fid=fopen(Growthfilename,'wt');
        fprintf(fid, [ header1 '\t' header2 '\t' headerlaser '\t' header3 '\t' header4 '\t' header5 '\t' 'Smoothed ExpNormDecon Signal' '\t' '1st Derivative of Smoothed' '\t' '2nd Derivative of Smoothed' '\t' 'Smoothed 2nd Derivative of Smoothed' '\t' header6 '\t' header7 '\t' header8 '\r\n']);
        fprintf(fid, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\r\n', [Growtht GrowthS GrowthR Growthdeconv Growthnorm Growthdiv Growth_smoothed Deriv1_Growth_smoothed Deriv2_Growth_smoothed Deriv2_Growth_doublesmoothed GrowthTemperature GrowthSetpoint GrowthHeaterPower]');
        fclose(fid);true;
    end
    cd(folderpath)
end