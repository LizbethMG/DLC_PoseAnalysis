%%% Signal analysis, lowpass filter 
%%% Nan->0 and Nan removed look equal. Withcam filtered looks as
%%% withoutcam_filtered
clc 
close all

%% Parameters to modify
%Led ON start of the experiment in seconds
% #todo: only works for multipe of 40 ms

%tmp_start = 12.8; % M45_C4:
tmp_start = 14.04; % M567_C4

%%

Cc_array=table2array(Cc);
Sry_array=table2array(Sry);
nb_total_frame=size(Cc_array,1);

spd=zeros(nb_total_frame,4);
spd(1,:)=missing;
tilt=[];
%Colors 
yellow = '#EDB120';
orange = "#D95319";
purple = 	"#7E2F8E";
green = "#77AC30";
light_blue = "#4DBEEE";


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

plot(x_wcam,y_wcam, 'k',  'DisplayName'," Without camera switch removed" ) % Before camera switch removed
hold on
plot(x,y,'Color', yellow, 'DisplayName'," Camera switch removed") % After camera switch removed
hold off

title('Instant speed')
xlabel('Time (s)')
ylabel('Pixels/s')
legend show


%% 2b) Camera switch + Median Filter 
%Filtered, Nan not removed 
n= 21;%18; % Window 720 ms
y_medfilt=medfilt1(y,n);
y_medfilt_omitnan=medfilt1(y_wcam,n, 'omitnan','truncate'); 

figure()
plot(x,y,'Color', yellow, 'DisplayName'," Camera switch removed") % After camera switch removed
hold on
plot(x,y_medfilt_omitnan, 'Color', green, 'DisplayName'," MedFilt21 omitnan")
hold on 
plot(x,y_medfilt, 'Color', orange, 'DisplayName'," MedFilt21")
hold on 

title('Instant speed - Camera swtich removed')
xlabel('Time (s)')
ylabel('Pixels/s')

legend show


%% 2c) Camera swith + Outliers removal 

[y_wo_outliers] = movmean(y_medfilt_omitnan,21);

figure()

plot(x,y_medfilt_omitnan, 'Color', green, 'DisplayName'," MedFilt21 omitnan")
hold on
plot(x,y_wo_outliers,'Color', purple,  'DisplayName'," Movmean" ) % Before camera switch removed
hold off

title('Instant speed')
xlabel('Time (s)')
ylabel('Pixels/s')
legend show

%% For Editrice figure only 

A = y_wo_outliers;
B = fillmissing(A,'movmedian',50);
F = smoothdata(B,'gaussian',15);
plot(x,F,'Color', purple,  'DisplayName'," Movmean" ) % Before camera switch removed
hold off

title('Instant speed')
xlabel('Time (s)')
ylabel('Pixels/s')
legend show


%% 3) Extract experimental data 


duration = 30*60;
tmp_end = tmp_start + 30*60; % True for ON OFF stimulation with 30 min duration
shifted_x = x-tmp_start;
idx_start= find(x==tmp_start); % index start experiment
idx_end = find(x==tmp_end); % index end experiment

% Extract data
%x_interm = (x-tmp_start);
%x_new = x_interm(1:30*60*25);
%y_new = y_wo_outliers(idx_start:idx_end);
x_new = shifted_x(shifted_x>=0 & shifted_x<duration);
y_new = y(shifted_x>=0 & shifted_x<duration);
y_NaN = isnan(y_new);

OFF1lim= 1:3*60*25; % 0 - 3 min
ON1lim = 3*60*25: 2*3*60*25; % 3- 6 min
OFF2lim = 2*3*60*25: 3*3*60*25; % 6- 9min
ON2lim = 3*3*60*25: 4*3*60*25; % 9- 12min
OFF3lim = 4*3*60*25: 5*3*60*25; % 12- 15min
ON3lim = 5*3*60*25: 6*3*60*25; % 15- 18min
OFF4lim = 6*3*60*25: 7*3*60*25; % 18- 21min
ON4lim = 7*3*60*25: 8*3*60*25; % 21- 24min
OFF5lim = 8*3*60*25: 9*3*60*25; % 24- 27min
ON5lim = 9*3*60*25: 10*3*60*25; % 27- 30min
figure()

OFF1 = area(x_new(OFF1lim),y_new(OFF1lim), 'FaceColor',orange);
hold on
ON1 = area(x_new(ON1lim),y_new(ON1lim), 'FaceColor',light_blue);
hold on
OFF2 = area(x_new(OFF2lim),y_new(OFF2lim), 'FaceColor',orange);
hold on
ON2 = area(x_new(ON2lim),y_new(ON2lim), 'FaceColor',light_blue);
hold on
OFF3 = area(x_new(OFF3lim),y_new(OFF3lim), 'FaceColor',orange);
hold on
ON3 = area(x_new(ON3lim),y_new(ON3lim), 'FaceColor',light_blue);
hold on
OFF4 = area(x_new(OFF4lim),y_new(OFF4lim), 'FaceColor',orange);
hold on
ON4 = area(x_new(ON4lim),y_new(ON4lim), 'FaceColor',light_blue);
hold on
OFF5 = area(x_new(OFF5lim),y_new(OFF5lim), 'FaceColor',orange);
hold on
ON5 = area(x_new(ON5lim),y_new(ON5lim), 'FaceColor',light_blue);
hold off

nsamplePerTry = 3*60*25;
NanRateOFF1 = sum(y_NaN(OFF1lim))/nsamplePerTry;
NanRateON1 = sum(y_NaN(ON1lim))/nsamplePerTry;
NanRateOFF2 = sum(y_NaN(OFF2lim))/nsamplePerTry;
NanRateON2 = sum(y_NaN(ON2lim))/nsamplePerTry;
NanRateOFF3 = sum(y_NaN(OFF3lim))/nsamplePerTry;
NanRateON3 = sum(y_NaN(ON3lim))/nsamplePerTry;
NanRateOFF4 = sum(y_NaN(OFF4lim))/nsamplePerTry;
NanRateON4 = sum(y_NaN(ON4lim))/nsamplePerTry;
NanRateOFF5 = sum(y_NaN(OFF5lim))/nsamplePerTry;
NanRateON5 = sum(y_NaN(ON5lim))/nsamplePerTry;


y_0forNan = y_new;
y_0forNan(y_NaN)=0;
MeanSpeedOFF1 = sum(y_0forNan(OFF1lim))/sum(~y_NaN(OFF1lim));
MeanSpeedON1 = sum(y_0forNan(ON1lim))/sum(~y_NaN(ON1lim));
MeanSpeedOFF2 = sum(y_0forNan(OFF2lim))/sum(~y_NaN(OFF2lim));
MeanSpeedON2 = sum(y_0forNan(ON2lim))/sum(~y_NaN(ON2lim));
MeanSpeedOFF3 = sum(y_0forNan(OFF3lim))/sum(~y_NaN(OFF3lim));
MeanSpeedON3 = sum(y_0forNan(ON3lim))/sum(~y_NaN(ON3lim));
MeanSpeedOFF4 = sum(y_0forNan(OFF4lim))/sum(~y_NaN(OFF4lim));
MeanSpeedON4 = sum(y_0forNan(ON4lim))/sum(~y_NaN(ON4lim));
MeanSpeedOFF5 = sum(y_0forNan(OFF5lim))/sum(~y_NaN(OFF5lim));
MeanSpeedON5 = sum(y_0forNan(ON5lim))/sum(~y_NaN(ON5lim));

Mean_off = mean ([MeanSpeedOFF1,MeanSpeedOFF2,MeanSpeedOFF3,MeanSpeedOFF4,MeanSpeedOFF5] );
Mean_on = mean ([MeanSpeedON1,MeanSpeedON2,MeanSpeedON3,MeanSpeedON4,MeanSpeedON5]);

%% Standard deviation of 3 min blocks 

tsstd = std(y_new(OFF1lim),'omitnan') 
OFF1 = area(x_new(OFF1lim),y_new(OFF1lim), 'FaceColor',orange);
hold on
yline(tsstd)
