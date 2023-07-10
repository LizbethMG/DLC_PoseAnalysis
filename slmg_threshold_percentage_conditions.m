function [On_percentages, Off_percentages, ON_mean, OFF_mean] = slmg_threshold_percentage_conditions (signal, t, tmp_start, filter_opt,  threshold)

% This is made exclusevely for the ONOFF experiments of 30 min duration
% with alternating On and Off conditions
x = t; 
y= signal;

duration = 30*60;
%tmp_end = tmp_start + 30*60; % True for ON OFF stimulation with 30 min duration
shifted_x = x-tmp_start;

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


    Off_1 = slmg_threshold_percentage (y_new(OFF1lim), threshold);
    On_1 = slmg_threshold_percentage (y_new(ON1lim), threshold);

    Off_2 = slmg_threshold_percentage (y_new(OFF2lim), threshold);
    On_2 = slmg_threshold_percentage (y_new(ON2lim), threshold);

    Off_3 = slmg_threshold_percentage (y_new(OFF3lim), threshold);
    On_3 = slmg_threshold_percentage (y_new(ON3lim), threshold);

    Off_4 = slmg_threshold_percentage (y_new(OFF4lim), threshold);
    On_4 = slmg_threshold_percentage (y_new(ON4lim), threshold);

    Off_5 = slmg_threshold_percentage (y_new(OFF5lim), threshold);
    On_5 = slmg_threshold_percentage (y_new(ON5lim), threshold);

    On_percentages = [On_1; On_2; On_3; On_4; On_5];
    Off_percentages = [Off_1; Off_2; Off_3; Off_4; Off_5];

    ON_mean_below = mean (On_percentages(:,1));
    ON_mean_above = mean (On_percentages(:,2));
    ON_mean_nan = mean (On_percentages(:,3));

    OFF_mean_below = mean (Off_percentages(:,1));
    OFF_mean_above = mean (Off_percentages(:,2));
    OFF_mean_nan = mean (Off_percentages(:,3));

    ON_mean = [ ON_mean_below; ON_mean_above; ON_mean_nan];
    OFF_mean = [ OFF_mean_below; OFF_mean_above; OFF_mean_nan];

else
    disp ("todo")
end 
end
