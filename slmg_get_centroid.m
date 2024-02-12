function [] = slmg_get_centroid(path_to_experiment)
%PATH_TO_EXPERIMENT Function based on Mona's Code_centroid_all_bp but
%wrapped in a function.
%   Analyze csv files from deeplabcut and return coordinates and likelihood per camera,choice of camera,
%   coordinates of bodyparts with best cam at every frame, and coordinates
%   of centroid.
%   All bp are taken into account
%   !  Make sure the files in the folder are ordered by name, and not by date
%    1st subsection and last subsection requires changing filepaths, if needed
%    This code removes pts with lk<0.7


%%% Analyze csv files from deeplabcut and return coordinates and likelihood per camera,choice of camera,
%%% coordinates of bodyparts with best cam at every frame, and coordinates
%%% of centroid. 
%%% All bp are taken into account
%%% ! Make sure the files in the folder are ordered by name, and not by date 
%%% 1st subsection and last subsection requires changing filepaths, if needed
%%% This code removes pts with lk<0.7

%% Change paths according to experiment 
% folder in which all csv files to analyze are stored (remove all other csv files)

% One folder per mouse
 if exist(path_to_experiment, 'dir')
    cd(path_to_experiment);
    destination_path =path_to_experiment;

else
    fprintf('The specified path does not exist or is not a directory.\n');
    return; % Exit the function
 end



%% Import csv files in a cell and shorten name 
% Sry : summary : choice of cam and coordinates given best camera 
% Cc : centroid coordinates 
% Cell : stores tables of coordinates obtained with each camera



csvfiles=dir('*csv');
headers_table=["Scorer","Frames_ms","Sn_x","Sn_y","Sn_lk","lear_x","lear_y",...
    "lear_lk","rear_x","rear_y","rear_lk","he_x","he_y","he_lk","bc_x","bc_y",...
    "bc_lk","tb_x","tb_y","tb_lk",...
    "lfrontp_x","lfrontp_y","lfrontp_lk",...
    "rfrontp_x","rfrontp_y","rfrontp_lk","lbackl_x","lbackl_y","lbackl_lk",...
    "rbackl_x","rbackl_y","rbackl_lk",...
    "lk_corr_sn","lk_corr_lear","lk_corr_rear",...
    "lk_corr_he","lk_corr_bc","lk_corr_tb","lk_corr_lfrontp","lk_corr_rfrontp","lk_corr_lbackl","lk_corr_rbackl","Sum_lk_corr"];

headers_Sry=["Scorer","Frames_ms","Ch?","Cam_ID_if_same_lk","Ch>3s","Same_lk","Sn_x","Sn_y","Sn_lk","lear_x","lear_y",...
    "lear_lk","rear_x","rear_y","rear_lk","he_x","he_y","he_lk","bc_x","bc_y",...
    "bc_lk","tb_x","tb_y","tb_lk","lfrontp_x","lfrontp_y","lfrontp_lk",...
    "rfrontp_x","rfrontp_y","rfrontp_lk","lbackl_x","lbackl_y","lbackl_lk",...
    "rbackl_x","rbackl_y","rbackl_lk"];

headers_Cc=["Scorer","Frames_ms","x_centroid","y_centroid","Pt_1","Pt_2","Pt_3","Pt_4","Pt_5","Pt_6","Pt_7","Pt_8","Pt_9","Pt_10"];
Cell=cell(2,length(csvfiles));
% Returns cell = name of csv files as columns and their corresponding tables as rows 
for i=1:length(csvfiles)
    T=readtable(csvfiles(i).name);
    ix=strfind(csvfiles(i).name,'_');
    vidname=[csvfiles(i).name(1:4) csvfiles(i).name(ix(4):ix(5)-1)];
    Cell{1,i}=vidname;
    Cell{2,i}=T;
end

%% Calculate corrected lk for each csvfile and add them at the end of the table 
for k=1:size(Cell,2)  %k=nb of csv files ; j= nb of lk_corr to calculate ; i=row (frame)
    T=Cell{2,k};
    Frames_ms=T{:,1}./0.025;
    T=[T(:,1) array2table(Frames_ms) T(:,2:end)];
    Tarray=table2array(T);
    index=5:3:32; % cols corrsponding to lk in table 
    lk_corr=zeros(size(Tarray,1),length(index));
    for j=1:length(index)

        if Tarray(1,index(j))>0.7 && Tarray(2,index(j))>0.7
            lk_corr(1,j)=1;
        else
            lk_corr(1,j)=0;
        end

        for i=2:size(Tarray,1)-1
            if (Tarray(i,index(j))>0.7 && (Tarray(i-1,index(j))>0.7||Tarray(i+1,index(j))>0.7))||(Tarray(i-1,index(j))>0.7 && Tarray(i+1,index(j))>0.7)
                lk_corr(i,j)=1;
            else
                lk_corr(i,j)=0;
            end
        end
        
        if Tarray(size(Tarray,1)-1,index(j))>0.7 && Tarray(size(Tarray,1),index(j))>0.7
            lk_corr(size(Tarray,1),j)=1;
        else
            lk_corr(size(Tarray,1),j)=2;
        end
    end
    Sum_lk_corr=sum(lk_corr,2);
    %Tarray=[Tarray lk_corr Sum_lk_corr];
    Tarray(:,end+1:end+10)=lk_corr;
    Tarray(:,end+1)=Sum_lk_corr;
    T=array2table(Tarray);
    T = renamevars(T,1:width(T),headers_table);
    Cell{2,k}=T;

end

%% Calculate choice of camera, coordinates given best cam and coordinates centroid
%Results : stores Sry as row 2 and centroid as row 3. Row 1 is for labels
%and each col is an experiment 
warning('off','all')
Results={};
Results{1,1}="Experiment";
Results{2,1}="Coord_bodyparts";
Results{3,1}="Coord_centroid";
Results{4,1}="Analysis_table";
Results{5,1}="Percent_bodypart_per_cam";
Results{6,1}="Nb_bp_per_time";

for nb_condition_expm=1:size(Cell,2)/2
    num_row=[size(Cell{2,nb_condition_expm},1) size(Cell{2,size(Cell,2)./2+nb_condition_expm},1)];
    nb_rows_summary=min(num_row);
    T1=table2array(Cell{2,nb_condition_expm});
    T2=table2array(Cell{2,size(Cell,2)./2+nb_condition_expm});
    T1(nb_rows_summary+1:end,:)=[];
    T2(nb_rows_summary+1:end,:)=[];
    Sry=[];
    Sry(:,1:2)=T1(:,1:2);
    % Sry : Scorer (numero of frame) -> col 2
    % Frame (ms) -> col 2 
    % Cam ? -> col 3

    for i=1:nb_rows_summary
        if T1(i,43)>T2(i,43)
            Sry(i,3)=1;
        else
            Sry(i,3)=2;
        end
    end

    %Same cam if lk equals -> col 4
    if T1(1,43)>T2(1,43)
        Sry(1,4)=1;
    else
        Sry(1,4)=2;
    end

    for i=2:nb_rows_summary
        if (T1(i,43)>T2(i,43))||((T1(i,43)==T2(i,43))&&Sry(i-1,4)==1)
            Sry(i,4)=1;
        else
            Sry(i,4)=2;
        end
    end

    
    % Best cam, at least >3s recording -> col 5

    if (Sry(1,4)==1 && sum(Sry(i+1:i+74,4)==74)) ||(Sry(1,4)==1 && sum(Sry(i+1:i+74,4))~=148)
        Sry(1,5)=1;
    else
        Sry(1,5)=2;
    end

    for i=2:nb_rows_summary-74
        if (Sry(i,4)==1 && sum(Sry(i+1:i+74,4))==74)||(Sry(i,4)==1 && Sry(i-1,5)==1)||(Sry(i,4)==2 && sum(Sry(i+1:i+74,4))~=148 && Sry(i-1,5)==1)
            Sry(i,5)=1;
        else
            Sry(i,5)=2;
        end
    end

    for i=nb_rows_summary-74:nb_rows_summary
        if Sry(i,4)==1
           Sry(i,5)=1;
        else
           Sry(i,5)=2;
        end
    end

    % Same lk between two cams ? -> col 6
    for i=1:nb_rows_summary
        if T1(i,43)==T2(i,43)
            Sry(i,6)=1;
        else
            Sry(i,6)=0;
        end
    end

    % Coordinates given best cam
    for j=3:32
        for i=1:nb_rows_summary
            if Sry(i,5) ==1
                Sry(i,j+4)=T1(i,j);
            else
                Sry(i,j+4)=T2(i,j);
            end
        end
    end

    % Centroid calculus


    % Calculus coord_x (cox) and coord_y (coy) while removing points for
    % which lk<0.7
    rowcol=zeros(size(Sry,1),10);
    Coxlk=zeros(size(Sry,1),10);
    Coylk=zeros(size(Sry,1),10);
    Cc=zeros(size(Sry,1),14);
    for num_frame=1:size(Sry,1)
        Lk_find=find(Sry(num_frame,9:3:36)>0.7);
        Lk_find(end+1:10)=missing;
        rowcol(num_frame,:)= Lk_find;

        Coxlk_find=Sry(num_frame,4+3*rmmissing(rowcol(num_frame,:)));
        Coxlk_find(:,end+1:10)=missing;
        Coxlk(num_frame,:)= Coxlk_find;
        Coylk_find=Sry(num_frame,5+3*rmmissing(rowcol(num_frame,:)));
        Coylk_find(:,end+1:10)=missing;
        Coylk(num_frame,:)= Coylk_find;
    

        % Calculate every combinaison of pts and area of their polyshape
        Coxlkrm=rmmissing(Coxlk(num_frame,:));
        Coylkrm=rmmissing(Coylk(num_frame,:));
        R=[];
        v=1:length(Coxlkrm);
        A=3:1:length(Coxlkrm);


        for nb_points_dans_combi=1:length(A)
            K=nchoosek(v,A(nb_points_dans_combi));
            L=K;
            L(:,end+1:10)=missing;

            for nb_combi_par_pts=1:size(K,1)
                x=Coxlkrm(:,K(nb_combi_par_pts,1:end));
                y=Coylkrm(:,K(nb_combi_par_pts,1:end));
                R=[R;polyarea(x,y) L(nb_combi_par_pts,1:end)];
            end
        end

        % Calculate Coord_centroid for best combinaison of pts (best
        % combinaison of pts = the one with the biggest area)  
        %Cc : col 1 : num frames
        %     col 2 : frames (ms)
        %     col 3 : x_centroid
        %     col 4 : y_centroid
        %     col 5 to 10 : points of the polygon 
        % 1-> Sn, 2->left ear, 3-> right_ear, 4-> head, 5-> body_center (hump),
        % 6-> tailbase
        if isempty(R)==1
            Cc(num_frame,:)=missing;
        else
            [rowmax,~]=find(R==max(R(1:end,1)));
            xmax=Coxlk(num_frame,rmmissing(R(rowmax,2:end)));
            ymax=Coylk(num_frame,rmmissing(R(rowmax,2:end)));
            %pgns=polyshape(Coxlk(7,:),Coylk(7,:),'Simplify',false,'KeepCollinearPoints',true); % pg non simplifié est meilleur ?
            %pgs=polyshape(Coxlk(7,:),Coylk(7,:),'Simplify',true,'KeepCollinearPoints',false);
            pg=polyshape(xmax,ymax);
            [xc,yc]=centroid(pg);
            %Cc=[Cc;(j-1) (j-1).*40 xc yc];
            pos=rowcol(num_frame,rmmissing(R(rowmax,2:end)));
            pos(:,end+1:10)=missing;
            Cc(num_frame,:)=[(num_frame-1) (num_frame-1).*40 xc yc pos];
        end
    end


% Freq analysis
Cc_array=Cc;
Sry_array=Sry;
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

%%% without cam
spd_wo_cam=spd; 
spd_wo_cam(1+(tilt/40),:)=missing;

spd_max=max(spd_wo_cam(:,2));
x=spd_wo_cam(:,1)/1000;
y=spd_wo_cam(:,2)/spd_max;

y_filtered=medfilt1(y,18);
y_wo_movmean=smoothdata(y_filtered,'movmean');

Analysis_table=table(x,y,y_filtered,y_wo_movmean);


% Camchoice and percent_bp
nb_total_frame=size(Cc_array,1);

    %Code plot switch cam 
    % percent bodypart per camera 
    percent_bodypart_per_cam=zeros(nb_total_frame,6);
    for frame=1:nb_total_frame
        cam=Sry_array(frame,5);
        percent_bodypart_per_cam(frame,1)=(frame-1)*40;
        percent_bodypart_per_cam(frame,2)=cam;
        if cam==1
            ind=nb_condition_expm;
        else 
            ind=size(Cell,2)./2+nb_condition_expm;
        end
        percent_bodypart_per_cam(frame,3)=table2array(Cell{2,ind}(frame,43));
        percent_bodypart_per_cam(frame,4)=(table2array(Cell{2,ind}(frame,43)))/0.1;
        percent_bodypart_per_cam(frame,5)=(table2array(Cell{2,nb_condition_expm}(frame,43)))/0.1;
        percent_bodypart_per_cam(frame,6)=(table2array(Cell{2,size(Cell,2)./2+nb_condition_expm}(frame,43)))/0.1;
    end

    %plot nb_pts per time 

    nb_bp_per_time=zeros(3,13);
        for nb_pts=0:1:10
         nb_bp_per_time(1,nb_pts+1)=sum(percent_bodypart_per_cam(:,3)==nb_pts);
         nb_bp_per_time(2,nb_pts+1)=sum(percent_bodypart_per_cam(:,3)==nb_pts)*40;
           nb_bp_per_time(3,nb_pts+1)=sum(percent_bodypart_per_cam(:,3)==nb_pts)*100/nb_total_frame;
        end
    nb_bp_per_time(1,12)=sum(nb_bp_per_time(1,1:6));
    nb_bp_per_time(2,12)=sum(nb_bp_per_time(2,1:6));
    nb_bp_per_time(3,12)=sum(nb_bp_per_time(3,1:6));
    nb_bp_per_time(1,13)=sum(nb_bp_per_time(1,7:11));
    nb_bp_per_time(2,13)=sum(nb_bp_per_time(2,7:11));
    nb_bp_per_time(3,13)=sum(nb_bp_per_time(3,7:11)); % j'aurais aussi pu faire nb(i,8)=sum..(i,1:4)

    %Create tables 
    headers_percent_bodypart=["Frames_ms","Cam","Nb_Visible_pts","%_Visible_points_best_cam","%_Visible_points_cam_1","%_Visible_points_cam_2"];
    headers_nb_bp_per_time=["0pt","1pt","2pts","3pts","4pts","5pts","6pts","7pts","8pts","9pts","10pts","Sum_<6pts","Sum_>=6pts"];
    percent_bodypart_per_cam=renamevars(array2table(percent_bodypart_per_cam),1:width(array2table(percent_bodypart_per_cam)),headers_percent_bodypart);
    nb_bp_per_time=renamevars(array2table(nb_bp_per_time),1:width(array2table(nb_bp_per_time)),headers_nb_bp_per_time);

%Save results in a cell 
Sry=renamevars(array2table(Sry),1:width(array2table(Sry)),headers_Sry);
Cc=renamevars(array2table(Cc),1:width(array2table(Cc)),headers_Cc);
Results{1,nb_condition_expm+1}=csvfiles(nb_condition_expm).name(ix(4):ix(5)-1);
Results{2,nb_condition_expm+1}=Sry;
Results{3,nb_condition_expm+1}=Cc;
Results{4,nb_condition_expm+1}=Analysis_table;
Results{5,nb_condition_expm+1}=percent_bodypart_per_cam;
Results{6,nb_condition_expm+1}=nb_bp_per_time;
end

%% Export Results (each expm condition is in the same excel file)
%Export Sry and Cc, analysis_table, percent_bp, nb_bp_per_time

for i=2:size(Results,2)
    writetable(Results{2,i},fullfile(destination_path,...
        sprintf('Results%s.xlsx',Results{1,i})),'Sheet',Results{2,1});
    writetable(Results{3,i},fullfile(destination_path,...
        sprintf('Results%s.xlsx',Results{1,i})),'Sheet',Results{3,1});
    writetable(Results{4,i},fullfile(destination_path,...
        sprintf('Results%s.xlsx',Results{1,i})),'Sheet',Results{4,1});
    writetable(Results{5,i},fullfile(destination_path,...
        sprintf('Results%s.xlsx',Results{1,i})),'Sheet',Results{5,1});
    writetable(Results{6,i},fullfile(destination_path,...
        sprintf('Results%s.xlsx',Results{1,i})),'Sheet',Results{6,1});
end

%Export coord and lk for each cam 
for i=1:size(Cell,2)
    if i<1+(size(Cell,2)/2)
        writetable(Cell{2,i},fullfile(destination_path,...
            sprintf('Results%s.xlsx',Results{1,i+1})),'Sheet',Cell{1,i});
    else
        writetable(Cell{2,i},fullfile(destination_path,...
            sprintf('Results%s.xlsx',Results{1,1+i-(size(Cell,2)/2)})),'Sheet',Cell{1,i});
    end
end


%Export Analysis table 
for i=2:size(Results,2)
    %writetable(Results{4,i},fullfile('C:\Users\mona.nashashibi\Work\Matlab\Tests\Results',...
     writetable(Results{4,i},fullfile(destination_path , sprintf('Analysis_table%s.txt',Results{1,i})));
end

sprintf('End')

end

