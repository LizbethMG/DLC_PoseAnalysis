%%% Signal analysis, lowpass filter 
%%% Nan->0 and Nan removed look equal. Withcam filtered looks as
%%% withoutcam_filtered
%Cc_array=Cc;
%Sry_array=Sry;
Cc_array=table2array(Cc);
Sry_array=table2array(Sry);
nb_total_frame=size(Cc_array,1);

spd=zeros(nb_total_frame,4);
spd(1,:)=missing;
tilt=[];
for frame=2:nb_total_frame
    spd(frame,3)=(Cc_array(frame,3)-Cc_array(frame-1,3))/(Cc_array(frame,2)-Cc_array(frame-1,2));
    spd(frame,4)=(Cc_array(frame,4)-Cc_array(frame-1,4))/(Cc_array(frame,2)-Cc_array(frame-1,2));
    spd(frame,2)=sqrt((spd(frame,3)^2)+(spd(frame,4)^2));
    spd(frame,1)=(frame-1)*40;
     if Sry_array(frame,5)~=Sry_array(frame-1,5)
         tilt(end+1)=(frame-1)*40;
     end

end
%%%with cams 
figure(1)
x_wcam=spd(:,1)/1000;
y_wcam=spd(:,2);
plot(spd(:,1)/1000,y_wcam)
y_wcam_filtered=medfilt1(y_wcam,18);
figure(2)
plot(x_wcam,y_wcam_filtered)

%%% without cam
spd_wo_cam=spd; 
spd_wo_cam(1+(tilt/40),:)=missing;

spd_max=max(spd_wo_cam(:,2));
x=spd_wo_cam(:,1)/1000;
y=spd_wo_cam(:,2)/spd_max;


%%% unfiltered, Nan replaced by 0
%x=fillmissing(x,"constant",0);
%y=fillmissing(y,"constant",0);
figure(3)
plot(x,y)

%%% unfiltered, Nan removed
wo_cam_missing=rmmissing(spd_wo_cam);
x_missing=wo_cam_missing(:,1)/1000;
y_missing=wo_cam_missing(:,2)/max(wo_cam_missing(:,2));
figure(4)
plot(x_missing,y_missing)

%%% filtered, Nan not removed 
figure(5)
y_filtered=medfilt1(y,18);
plot(x,y_filtered)

%%% filtered, Nan removed  
y_missing_filtered=medfilt1(y_missing,18);
figure(6)
plot(x_missing,y_missing_filtered)


%% Smoothing 

%%% movmean 
y_wo_movmean=smoothdata(y_filtered,'movmean',18);
figure(7)
plot(x,y_wo_movmean)


%%%Median
y_wcf_median=smoothdata(y_wcam_filtered,'movmedian');
figure(8)
plot(x_wcam,y_wcf_median)

%%%Lowess
y_wcf_lowess=smoothdata(y_wcam_filtered,'lowess');
figure(9)
plot(x_wcam,y_wcf_lowess)
%%%Loess
y_wcf_loess=smoothdata(y_wcam_filtered,'loess');
figure(10)
plot(x_wcam,y_wcf_loess)

%%%rlowess
y_wcf_rlowess=smoothdata(y_wcam_filtered,'rlowess');
figure(11)
plot(x_wcam,y_wcf_rlowess)

%%%rloess
y_wcf_rloess=smoothdata(y_wcam_filtered,'rloess');
figure(12)
plot(x_wcam,y_wcf_rloess)


%% Mean 
nb_of_periodes=input("nb of periodes :");
spd_mean=zeros(1,nb_of_periodes);
spd_mean(1)=mean(rmmissing(y_wo_movmean(1:round(length(y_wo_movmean)/nb_of_periodes))));
for i=2:nb_of_periodes-1
spd_mean(i)=mean(rmmissing(y_wo_movmean(1+round((i-1)*length(y_wo_movmean)/nb_of_periodes):round(i*length(y_wo_movmean)/nb_of_periodes))));
end 
spd_mean(end)=mean(rmmissing(y_wo_movmean(1+round(i*length(y_wo_movmean)/nb_of_periodes):end)));

%% Export data 
Analysis_table=table(x,y,y_filtered,y_wo_movmean);
spd_mean_table=array2table(spd_mean);
writetable(Analysis_table,fullfile(destination_path,'Analysis_table.txt'));
writetable(spd_mean_table,fullfile(destination_path,'Spd_mean.txt'));





