%%% Plotting the centroid as a video_ all bp
%%% you must have analyzed the video beforehand with Code_centroid_all_bp and kept the "Results"
%%% variable (that contains Sry and Cc) 

v_cam1=VideoReader('C:\Users\mona.nashashibi\Work\Matlab\Tests\Place_holder\ch01_20201013200922_001_avi_58AZQ_converted_1000_frames.avi');
v_cam2=VideoReader('C:\Users\mona.nashashibi\Work\Matlab\Tests\Place_holder\ch02_avi_58AZQDLC_converted_1000_frames.avi');
Input_videos={v_cam1,v_cam2};
Frame_v_cam1=Input_videos{1,1}.NumFrames;
Frame_v_cam2=Input_videos{1,2}.NumFrames;
numberOfFrames=999; %fill in the appropriate number
%numberOfFrames = min([Frame_v_cam1 Frame_v_cam2])-1;
%% Preallocation
allTheFrames = cell(numberOfFrames,1);
vidHeight = 576;
vidWidth = 704;
allTheFrames(:) = {zeros(vidHeight, vidWidth, 3, 'uint8')};
allTheColorMaps = cell(numberOfFrames,1);
allTheColorMaps(:) = {zeros(256, 3)};
F = struct('cdata', allTheFrames, 'colormap', allTheColorMaps);


%% Extract frames from vid


for frame = 1 : numberOfFrames  
    cam=Results{2,2}(frame,5);

    % Calculate polygon and coord_centroid (actually, this part is now
    % unnecessary because already calculated but I didn't adapt the code) 

    Centre=Results{3,2}(frame,3:4);%write Centre=table2array(Results{3,2}(7,3:4)) if necessary
    list_pts=rmmissing(table2array(Results{3,2}(frame,5:14)));
    Coordx=table2array(Results{2,2}(frame,(list_pts.*3)+4));
    Coordy=table2array(Results{2,2}(frame,(list_pts.*3)+5));
    polygon=polyshape(Coordx,Coordy);
    [xpg,ypg]=centroid(polygon);


    %create fig
    %this_frame = read(v, frame+2);
    this_frame=read(Input_videos{1,table2array(cam)},frame+2);
    thisfig = figure();
    set(thisfig,'MenuBar','none','Toolbar','none','Visible','off');
    thisax = axes('Parent', thisfig,'YDir','reverse','XLim',[0 704],'YLim',[0 576]);
    hold (thisax,'on');
    image(this_frame, 'Parent', thisax);
    plot1 = plot(polygon,'Parent',thisax,'FaceColor',[0 0.447 0.741],...
        'FaceAlpha',0.35);
    plot2=plot(xpg,ypg,'MarkerSize',5,'Marker','o','LineStyle','none','Color','w');
    %datatip(plot2);
    hold(thisax,'off');


    F(frame)=getframe(gcf);
end

%%
video = VideoWriter('test_plot_58AZQ_fusion_all_bp_voyons.avi');
video.FrameRate = 25;
open(video);
% write the frames to the video
for frame=1:length(F)
    % convert the image to a frame
    frame_2 = F(frame) ;
    writeVideo(video, frame_2);
end
close(video)