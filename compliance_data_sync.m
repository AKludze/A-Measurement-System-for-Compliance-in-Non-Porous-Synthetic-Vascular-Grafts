clear all; 
close all; 

%% STEP 1: LOAD DATA (Arduino & Micrometer Data)

filename_arduino = 'logging_compliance_data.csv' ; %type-in-save-file-name for arduino data
arduino_data = readtable(filename_arduino); %raw arduino data from PUTTY 
arduino_data.Properties.VariableNames([1 2 3]) = {'Pressure (mmHG)' 'Valve Step Function' 'Time (seconds)'}; 
%rename data and remove un-needed columns (e.g., insert Pressure, Water Valve-Step Function, and Time (seconds) in place of
% PuTTY log 2023.05.03 12:35:16, VAR 2, and VAR 3)

filename_micrometer = 'COMPLIANCE_TEST_2.xlsx';
micrometer_data = readtable(filename_micrometer);
micrometer_data.Properties.VariableNames([1 2]) = {'Time (seconds)' 'Diameter (mm)'}; 

%% STEP 2: CLEAN DATA (Arduino & Micrometer Data)

arduino_data(:,4:width(arduino_data)) = []; %delete the un-needed columns (e.g., 4-14) and keep columns 1-3)
arduino_data = table2array(arduino_data); %convert data table to matrix/array
arduino_data(isnan(arduino_data(:,1)),:)=[]; %remove those rows w/ NaN in column 1 (Pressure (mmHG))
arduino_data(isnan(arduino_data(:,2)),:)=[]; %remove those rows w/ NaN in column 2 (Valve Step Function)
arduino_data(isnan(arduino_data(:,3)),:)=[]; %remove those rows w/ NaN in column 3 (Time (seconds))

micrometer_data(:,3:width(micrometer_data)) = []; %delete the un-needed columns (e.g., 4-14) and keep columns 1-3)
micrometer_data = table2array(micrometer_data); %convert data table to matrix/array
micrometer_data(isnan(micrometer_data(:,1)),:)=[]; %remove those rows w/ NaN in column 1 (Time (seconds))
micrometer_data(isnan(micrometer_data(:,2)),:)=[]; %remove those rows w/ NaN in column 2 (Diameter (mm))
for n = 1:height(micrometer_data)
   micrometer_data(n,2) = micrometer_data(n,2);  % Copy original value
   micrometer_data(n,2) = micrometer_data(n,2) + 10; %Add 10 because laser micrometer measures ("number" - 10)
end

delay_index_1 = 105; %approximately 5 second delay from arduino data-collection
delay_index_2 = 50; %approximately 5 second delay from micrometer data-collection

arduino_data = arduino_data(delay_index_1:end,:);
micrometer_data = micrometer_data(delay_index_2:end,:);

%% STEP 3: FIND INITIAL PEAKS (Pressure, Diameter vs. Time) 
%Not Necessary but could help: Valve Step Function vs. Time 

arduino_time = arduino_data(:,3); %load arduino_data into separate components
arduino_valve = arduino_data(:,2);
arduino_pressure = arduino_data(:,1);

micrometer_time = micrometer_data(:,1);
micrometer_diameter = micrometer_data(:,2);

[pressure_max,time_max_1] = findpeaks(arduino_pressure,arduino_time,'MinPeakDistance',1.20); %find max pressure at specific time
[diameter_max,time_max_2] = findpeaks(micrometer_diameter,micrometer_time); %find max diameter at specific time 

[pressure_min,time_min_1] = findpeaks(-arduino_pressure,arduino_time,'MinPeakDistance',1.20); %find min pressure at specific time
[diameter_min,time_min_2] = findpeaks(-micrometer_diameter,micrometer_time); %find min diameter at specific time

diameter_min = abs(diameter_min); %get-real result
pressure_min = abs(pressure_min); %get-real result

pressure_max_n_min = cat(1,pressure_max,pressure_min); %combine data
time_max_n_min_1 = cat(1,time_max_1,time_min_1);

diameter_max_n_min = cat(1,diameter_max,diameter_min); %combine data
time_max_n_min_2 = cat(1,time_max_2,time_min_2);

pressure_time_max_min = cat(2,pressure_max_n_min,time_max_n_min_1); %combine pressure and time data
diameter_time_max_min = cat(2,diameter_max_n_min,time_max_n_min_2); %combine pressure and time data

pressure_time_max_min = sortrows(pressure_time_max_min,2); %sort data based on time
diameter_time_max_min = sortrows(diameter_time_max_min,2); %sort data based on time 

for n = 1:length(pressure_time_max_min) %if duplicate datapoints within the same interval 
    if n == length(pressure_time_max_min)
        break
    end 
   
    if abs((pressure_time_max_min(n,1)-pressure_time_max_min(n+1,1)))/(pressure_time_max_min(n,1)) < 0.0150
        temp_avg_pressure = mean(pressure_time_max_min(n:n+1,1));
        temp_avg_time = mean(pressure_time_max_min(n:n+1,2));
        pressure_time_max_min(n,:) = cat(2,temp_avg_pressure,temp_avg_time);
        pressure_time_max_min(n+1,:) = [];
    end 
end



%% STEP 4: SYNC DATA, ADJUST FOR TIME DELAY  (Pressure, Diameter vs. Time) BASED ON PEAKS

% ADJUST FOR TIME DELAY IN Pressure vs. Time
for m = 1:length(arduino_data)
% for each data point, if the time-point is less than the first peak-time
% point, then delete the data point
        if arduino_time(m) < pressure_time_max_min(1,2) 
            arduino_time(m) = NaN;
            arduino_pressure(m) = NaN;
%if the time-point is greater than the first peak-time correct by
%subtracting the first-peak time (time correct) 
        elseif arduino_time(m) == pressure_time_max_min(1,2)
            arduino_time(m) = arduino_time(m) - pressure_time_max_min(1,2);
%if the time-point is greater than the first peak-time correct by
%subtracting the first-peak time (time correct) 
        elseif arduino_time(m) > pressure_time_max_min(1,2)
            arduino_time(m) = arduino_time(m) - pressure_time_max_min(1,2);
        end
end

arduino_time(isnan(arduino_time(:,1)),:)=[]; %remove those rows w/ NaN in column 1 (Time (seconds))
arduino_pressure(isnan(arduino_pressure(:,1)),:)=[]; %remove those rows w/ NaN in column 1 (Time (seconds))

for j = 1:length(pressure_time_max_min)
    pressure_time_max_min(j,2) = pressure_time_max_min(j,2) - pressure_time_max_min(1,2);
end 

% ADJUST FOR TIME DELAY IN Diameter vs. Time
for m = 1:length(micrometer_data)
% for each data point, if the time-point is less than the first peak-time
% point, then delete the data point
        if micrometer_time(m) < diameter_time_max_min(1,2) 
            micrometer_time(m) = NaN;
            micrometer_diameter(m) = NaN;
%if the time-point is greater than the first peak-time correct by
%subtracting the first-peak time (time correct) 
        elseif micrometer_time(m) == diameter_time_max_min(1,2)
            micrometer_time(m) = micrometer_time(m) - diameter_time_max_min(1,2);
%if the time-point is greater than the first peak-time correct by
%subtracting the first-peak time (time correct) 
        elseif micrometer_time(m) > diameter_time_max_min(1,2)
            micrometer_time(m) = micrometer_time(m) - diameter_time_max_min(1,2);
        end
end

micrometer_time(isnan(micrometer_time(:,1)),:)=[]; %remove those rows w/ NaN in column 1 (Time (seconds))
micrometer_diameter(isnan(micrometer_diameter(:,1)),:)=[]; %remove those rows w/ NaN in column 1 (Time (seconds))

for j = 1:length(diameter_time_max_min)
    diameter_time_max_min(j,2) = diameter_time_max_min(j,2) - diameter_time_max_min(1,2);
end 

%% STEP 5: FIND PEAKS FOR SYNCED DATA (Pressure, Diameter vs. Time) 

%%re-find peaks based on adjusted-time
[pressure_max,time_max_1] = findpeaks(arduino_pressure,arduino_time,'MinPeakDistance',1.20); %find max pressure at specific time
[diameter_max,time_max_2] = findpeaks(micrometer_diameter,micrometer_time); %find max diameter at specific time 

[pressure_min,time_min_1] = findpeaks(-arduino_pressure,arduino_time,'MinPeakDistance',1.20); %find min pressure at specific time
[diameter_min,time_min_2] = findpeaks(-micrometer_diameter,micrometer_time); %find min diameter at specific time

diameter_min = abs(diameter_min); %get-real result
pressure_min = abs(pressure_min); %get-real result

pressure_max_n_min = cat(1,pressure_max,pressure_min); %combine data
time_max_n_min_1 = cat(1,time_max_1,time_min_1);

diameter_max_n_min = cat(1,diameter_max,diameter_min); %combine data
time_max_n_min_2 = cat(1,time_max_2,time_min_2);

pressure_time_max_min = cat(2,pressure_max_n_min,time_max_n_min_1); %combine pressure and time data
diameter_time_max_min = cat(2,diameter_max_n_min,time_max_n_min_2); %combine pressure and time data

pressure_time_max_min = sortrows(pressure_time_max_min,2); %sort data based on time
diameter_time_max_min = sortrows(diameter_time_max_min,2); %sort data based on time 
%%re-find peaks based on adjusted-time

for n = 1:length(pressure_time_max_min) %if duplicate datapoints within the same interval (e.g., 
    if n == length(pressure_time_max_min)
        break
    end 
   
    if abs((pressure_time_max_min(n,1)-pressure_time_max_min(n+1,1)))/(pressure_time_max_min(n,1)) < 0.0150
        temp_avg_pressure = mean(pressure_time_max_min(n:n+1,1));
        temp_avg_time = mean(pressure_time_max_min(n:n+1,2));
        pressure_time_max_min(n,:) = cat(2,temp_avg_pressure,temp_avg_time);
        pressure_time_max_min(n+1,:) = [];
    end 
end

%Organize Peak-Data to START with D_high and P_high

if diameter_time_max_min(1,1) < diameter_time_max_min(2,1)
    diameter_time_max_min(1,1) = NaN;
    diameter_time_max_min(1,2)= NaN;
end

diameter_time_max_min(isnan(diameter_time_max_min(:,1)),:)=[]; 

if pressure_time_max_min(1,1) < pressure_time_max_min(2,1)
    pressure_time_max_min(1,1) = NaN;
    pressure_time_max_min(1,2)= NaN;
end

pressure_time_max_min(isnan(pressure_time_max_min(:,1)),:)=[]; 

%Organize Peak-Data to END with D_low and P_low

if diameter_time_max_min(end,1) > diameter_time_max_min(end-1,1)
    diameter_time_max_min(end,1) = NaN;
    diameter_time_max_min(end,2)= NaN;
end

diameter_time_max_min(isnan(diameter_time_max_min(:,1)),:)=[]; 

if pressure_time_max_min(end,1) > pressure_time_max_min(end-1,1)
    pressure_time_max_min(end,1) = NaN;
    pressure_time_max_min(end,2)= NaN;
end

pressure_time_max_min(isnan(pressure_time_max_min(:,1)),:)=[]; 

if pressure_time_max_min(1,2) < diameter_time_max_min(1,2) %pressure start peak time is smaller than diameter start peak time (adjust diameter time)
    time_shift = diameter_time_max_min(1,2) - pressure_time_max_min(1,2); %calculate the time shift/delay between the two peaks

    for n = 1:length(diameter_time_max_min)
        diameter_time_max_min(n,2) =  diameter_time_max_min(n,2) - time_shift;
    end 

    for m = 1:length(micrometer_time)
        micrometer_time(m) = micrometer_time(m) - time_shift;
    end

%delete time and pressure data that have negative time

    temp_diameter_time_merge = cat(2,micrometer_time,micrometer_diameter);
    temp_diameter_time_merge(temp_diameter_time_merge(:,1)<0,:)=[];
    micrometer_time =  temp_diameter_time_merge(:,1);
    micrometer_diameter = temp_diameter_time_merge(:,2); 
end 

if pressure_time_max_min(1,2) > diameter_time_max_min(1,2) %pressure start peak time is greater than diameter start peak time (adjust pressure time)
    time_shift = pressure_time_max_min(1,2) - diameter_time_max_min(1,2); %calculate the time shift/delay between the two peaks

    for n = 1:length(pressure_time_max_min)
        pressure_time_max_min(n,2) =  pressure_time_max_min(n,2) - time_shift;
    end 

    for m = 1:length(arduino_time)
        arduino_time(m) = arduino_time(m) - time_shift;
    end

%delete time and pressure data that have negative time
    temp_pressure_time_merge = cat(2,arduino_time,arduino_pressure);
    temp_pressure_time_merge(temp_pressure_time_merge(:,1)<0,:)= [];
    arduino_time =  temp_pressure_time_merge(:,1);
    arduino_pressure = temp_pressure_time_merge(:,2); 
end 


%% STEP 6: DATA ANALYSIS: COMPLIANCE CALCULATION (Wu. W. et al (2012) Nature Medicine)


%Ensure vectors are same size
if length(diameter_time_max_min) ~= length(pressure_time_max_min)
    if length(diameter_time_max_min) > length(pressure_time_max_min)
        size_diff = length(diameter_time_max_min) - length(pressure_time_max_min);
        diameter_time_max_min = diameter_time_max_min(1:end-size_diff,:);
    end 

    else if length(diameter_time_max_min) < length(pressure_time_max_min)
        size_diff = length(pressure_time_max_min) - length(diameter_time_max_min); 
        pressure_time_max_min = pressure_time_max_min(1:end-size_diff,:);
    end 
end 


Compliance_Value = [];
Compliance_Time_Intervals = [];

for n = 1:2:length(diameter_time_max_min)
   if n == length(diameter_time_max_min) && rem(length(diameter_time_max_min),2) == 1 
       break %if the last index is ODD, quit the loop
   end 
   D_high = diameter_time_max_min(n,1);
   D_low =  diameter_time_max_min(n+1,1);
   P_high = pressure_time_max_min(n,1);
   P_low = pressure_time_max_min(n+1,1);
   t_i = diameter_time_max_min(n,2);
   t_f = diameter_time_max_min(n+1,2);
   C = (10000* (D_high-D_low))/((D_low)*(P_high-P_low)); 
   Compliance_Value(end+1) = C;
   Compliance_Time_Intervals(end+1) = t_f;
end


%% STEP 7: PLOT SYNC-ED DATA WITH PEAKS
figure(1);
plot(arduino_time,arduino_pressure, pressure_time_max_min(:,2), pressure_time_max_min(:,1), 'pg')
xlabel('Time (seconds)')
ylabel('Pressure (mmHG)')
title('Pressure (mmHG) vs. Time (seconds)')


figure(2);
plot(micrometer_time,micrometer_diameter,diameter_time_max_min(:,2), diameter_time_max_min(:,1),'pg')
xlabel('Time (seconds)')
ylabel('Diameter (mm)')
title('Diameter (mm) vs. Time (seconds)')

figure(3);
subplot(3,1,1,'align');
plot(arduino_time,arduino_pressure, pressure_time_max_min(:,2), pressure_time_max_min(:,1), 'pg')
ylim([min(pressure_time_max_min(:,1))-5 max(pressure_time_max_min(:,1))+5])
xlabel('Time (seconds)')
ylabel('Pressure (mmHG)')
%title('Pressure (mmHG) vs. Time (seconds)')

subplot(3,1,2, 'align');
plot(micrometer_time,micrometer_diameter,diameter_time_max_min(:,2), diameter_time_max_min(:,1),'pg')
ylim([min(diameter_time_max_min(:,1))-0.01 max(diameter_time_max_min(:,1))+0.01])
xlabel('Time (seconds)')
ylabel('Diameter (mm)')
%title('Diameter (mm) vs. Time (seconds)')

subplot(3,1,3, 'align');
%plot(Compliance_Time_Intervals,Compliance_Value)
plot(Compliance_Time_Intervals,Compliance_Value,'o')
ylim([min(Compliance_Value)-0.1 max(Compliance_Value)+0.1])
xlabel('Time (seconds)')
ylabel('Compliance (%/100 mmHg)')
%title('Compliance (%/100 mmHg) vs. Time (seconds)')