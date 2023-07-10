%%% Signal analysis, lowpass filter 
%%% Nan->0 and Nan removed look equal. Withcam filtered looks as
%%% withoutcam_filtered
clc 
close all

%% Parameters to modify

% M45 C4
% mouseID = 'M567-C4'; % Animal ID_experiment
% tmp_start = 12.8; % Led ON start of the experiment in seconds

% M186 C4
% mouseID = 'M186-C4'; % Animal ID_experiment
% tmp_start = 11.2; % Led ON start of the experiment in seconds

% M423 C4 * pretty graphs
mouseID = 'M423-C4'; % Animal ID_experiment
tmp_start = 98.76; % Led ON start of the experiment in seconds

% M451 C4
% mouseID = 'M451-C4'; % Animal ID_experiment
% tmp_start = 46;%46.52; % Led ON start of the experiment in seconds

% M567 C4
% mouseID = 'M567-C4'; % Animal ID_experiment
% tmp_start = 14.04; % Led ON start of the experiment in seconds

% M674 C4
% mouseID = 'M674-C4-2'; % Animal ID_experiment
% tmp_start = 25.12; % Led ON start of the experiment in seconds

% M749 C4
% mouseID = 'M749-C4'; % Animal ID_experiment
% tmp_start = 18.28; % Led ON start of the experiment in seconds

% M716 C4
% mouseID = 'M716-C4'; % Animal ID_experiment
% tmp_start = 15.72; % Led ON start of the experiment in seconds

% M763 C4
% mouseID = 'M763-C4'; % Animal ID_experiment
% tmp_start = 25.48; % Led ON start of the experiment in seconds

% M328 C4
% mouseID = 'M328-C4'; % Animal ID_experiment
% tmp_start = 71.56; % Led ON start of the experiment in seconds




%% Compute instant speed

% Colors for plots
yellow = '#EDB120';
orange = "#D95319";
purple = 	"#7E2F8E";
green = "#77AC30";
light_blue = "#4DBEEE";

Cc_array=table2array(Cc);
Sry_array=table2array(Sry);
nb_total_frame=size(Cc_array,1);

spd=zeros(nb_total_frame,4);
spd(1,:)=missing;
tilt=[];

% Compute instant speed 
for frame=2:nb_total_frame
    spd(frame,1)=((frame-1)*40)/1000; % time in s
    dx = Cc_array(frame,3)-Cc_array(frame-1,3); % x2-x1
    dy = Cc_array(frame,4)-Cc_array(frame-1,4); % y2-y1
    dt= (Cc_array(frame,2)-Cc_array(frame-1,2))/1000; % dt in s
    d= sqrt((dx^2)+(dy^2)); % Eucledian distance
    spd(frame,2)=d/dt;% Instantaneous velocity (x2-x1) / dt  in Pixels per second

     if Sry_array(frame,5)~=Sry_array(frame-1,5) % s'il y a un changement de camera
         tilt(end+1)=(frame-1)*40; % ajout dans tilt le timestamp du frame 
     end

end

x_wcam=spd(:,1);
y_wcam=spd(:,2);

%% 1)  Instant speed with camera switches 
figure(1)
plot(x_wcam,y_wcam, 'k')
title('Instant speed - raw data with camera switches')
xlabel('Time (s)')
ylabel('Pixels/s')

%% 2)  Correction for camera switch
% Adds NaN to frames of changing channels (camera switch)
spd_wo_cam = spd;
spd_wo_cam(1+(tilt/40),:)=missing; 

x = x_wcam;
y= spd_wo_cam(:,2);

figure()

plot(x,y,'Color', yellow, 'DisplayName'," Camera switch removed") % After camera switch removed
hold on 
plot(x_wcam,y_wcam, 'k',  'DisplayName'," Without camera switch removed" ) % Before camera switch remove

hold off

title('Instant speed')
xlabel('Time (s)')
ylabel('Pixels/s')
legend show


%% 2b) Camera swith + Outliers removal 
m = 21;
[y_wo_outliers] = movmean(y,m, 'omitnan');

figure()

plot(x, y, 'Color',light_blue, 'DisplayName', "Wo Camera switch")
hold on
plot(x,y_wo_outliers,'Color', purple,  'DisplayName'," Movmean" ) 
hold off

title('Instant speed')
xlabel('Time (s)')
ylabel('Pixels/s')
legend show


%% 2c) Camera switch + Median Filter 
%Filtered, Nan not removed 
n= 21;%18; % Window 720 ms
y_medfilt=medfilt1(y_wo_outliers,n);
y_medfilt_omitnan=medfilt1(y_wo_outliers,n, 'omitnan','truncate'); 

figure()
plot(x,y_wo_outliers,'Color', purple,  'DisplayName'," Movmean" ) % After moving filter
hold on
plot(x,y_medfilt, 'Color', orange, 'DisplayName'," MedFilt21") % After median filter
hold on
plot(x,y_medfilt_omitnan, 'Color', green, 'DisplayName'," MedFilt21 omitnan") % After median filter omitnan
hold off

title('Instant speed - Camera swtich removed')
xlabel('Time (s)')
ylabel('Pixels/s')

legend show



%% For Editrice figure only 

figure()
A = y_medfilt_omitnan;
B = fillmissing(A,'movmedian',50);
F = smoothdata(B,'gaussian',15);
plot(x,F,  'DisplayName'," Movmean" ) % Before camera switch removed
hold off

title('Instant speed (for editor)')
xlabel('Time (s)')
ylabel('Pixels/s')
legend show


%% 3) Extract experimental data 

t = x;
signal= F; 
filter_opt = 0; % No behaviour remove
b_list = [];
plot_opt = 1; 

[OFFs, ONs, Mean_on, Mean_off, NanRateONs, NanRateOFFs ] = slmg_average_speed(t, signal, tmp_start, filter_opt, b_list, mouseID, plot_opt);


%% Standard deviation of 3 min blocks 
t = x;
signal = F; 


filter_opt = 0; % No behaviour remove
threshold = std(signal,'omitnan') ;

figure () 
plot(x,F,'Color', purple,  'DisplayName'," Movmean" ) % Before camera switch removed
hold on
yline(threshold)
hold off

[On_percentages, Off_percentages, ON_mean, OFF_mean] = slmg_threshold_percentage_conditions (signal, t, tmp_start, filter_opt,  threshold)

