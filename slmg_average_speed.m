function [OFFs, ONs, Mean_on, Mean_off, NanRateONs, NanRateOFFs ] = slmg_average_speed(t, s, tmp_start, filter_opt, b_list, mouseID, plot_opt)
% Input:
%   t: time vector
%   s: speed vector 
%   tmp_start: start LED timestamp
%   filter_opt: 1 to remove behavioral events as in "b_list" or 0 to not
%   remove any window
%   mouseID: example: M567_C4
%   plot_opt: 1 to plot results 0 to plot nothing
% Output: 

% Colors for plots
yellow = '#EDB120';
orange = "#D95319";
purple = 	"#7E2F8E";
green = "#77AC30";
light_blue = "#4DBEEE";

x = t;
y= s;

duration = 30*60;
tmp_end = tmp_start + 30*60; % True for ON OFF stimulation with 30 min duration
shifted_x = x-tmp_start;
idx_start= find(x==tmp_start); % index start experiment
idx_end = find(x==tmp_end); % index end experiment
       
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

if filter_opt == 0 % Without removing specific behaviours 
   
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
    
    % Gives the percentage of nan values for each OFF or On period
    % This represent the bad DLC detection since a centroid can't be calculated
    % without less than 3 points
    NanRateONs = [NanRateON1 NanRateON2 NanRateON3 NanRateON4 NanRateON5 ].*100; 
    NanRateOFFs = [NanRateOFF1 NanRateOFF2 NanRateOFF3 NanRateOFF4 NanRateOFF5 ].*100; 
    
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
    
    OFFs = [MeanSpeedOFF1 MeanSpeedOFF2 MeanSpeedOFF3 MeanSpeedOFF4 MeanSpeedOFF5];
    ONs = [MeanSpeedON1 MeanSpeedON2 MeanSpeedON3 MeanSpeedON4 MeanSpeedON5];

    Mean_off = mean ([MeanSpeedOFF1,MeanSpeedOFF2,MeanSpeedOFF3,MeanSpeedOFF4,MeanSpeedOFF5] );
    Mean_on = mean ([MeanSpeedON1,MeanSpeedON2,MeanSpeedON3,MeanSpeedON4,MeanSpeedON5]);

     if plot_opt == 1 % plot data colored by condition (On blue and Off in red)
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

        hold on
        threshold = std(y_new,'omitnan') ;       
        yline(threshold)
        hold off

        title('Instant speed ', mouseID)
        xlabel('Time (s)')
        ylabel('Pixels/s')
    end
else
    display('todo remove behaviors')
end