clc, clear variables, close all

%% User variables

Author_Name = 'Nathan Nahrung';
Author_Email = 'nnahrung@gmail.com';
Track_Name = 'lakeside_track';
File_Name = 'model';% used for .sdf & .config
Path_width = 3;%greater than (error aprox = resolution) - Note: FSG Rule Competition Handbook 2020 DE6.3 - minimum 3 m
Cone_sep = 3;%less than (error aprox = resolution) - Note: FSG Rule Competition Handbook 2020 DE6.3 - cone sep should be less than 5 m and reduced at corners
clockwise = true;%false;%
resolution = 0.5;%grid unit in metres
interpulation = 4;%used to close gaps in GPS data
track_offset = [0, 0];%meters xy offset of track (by default center of start is [0, 0])
track_rotation = 0;%Starting direction is positive x, change this bearing anticlockwise with degrees
acceleration_test = true;%false;%
straight_length = 75; % acceleration_test variable (acceleration length)
break_length = 100; % acceleration_test variable (breaking length)
skidpad_test = false;%true;%
entrance_length = 6; % skidpad_test variable (lengh from start to 
mid_length = 15; % skidpad_test variable (length 
exit_length = 13; % skidpad_test variable 
big_cone_sep = 0.6;
Cone_Model_blue = '//qut_description/meshes/cone_blue.dae';
Cone_Model_yelw = '//qut_description/meshes/cone_yellow.dae';
Cone_Model_big = '//qut_description/meshes/cone_big.dae';
Cone_Model_orng = '//qut_description/meshes/cone.dae';

% Cone_Model_blue = '//eufs_description/meshes/cone_blue.dae';
% Cone_Model_yelw = '//eufs_description/meshes/cone_yellow.dae';
% Cone_Model_big = '//eufs_description/meshes/big_cone.dae';
% Cone_Model_orng = '//eufs_description/meshes/cone.dae';

%% Track center GPS Data (User variable, copy, modify and comment out as required)

% GPS data should be '.csv', column 1 latitude, column 2 longitude, column 3 altitude (not used, can be 0)

Data_Table = readtable('lakeside GPS.csv');
Data_Array = table2array(Data_Table);
Data_trimmed = Data_Array(5650:24780,1:3);

% Data_Table = readtable('winton GPS.csv');
% Data_Array = table2array(Data_Table);
% Data_trimmed = Data_Array(12160:13640,1:3);

% Data_Table = readtable('QUT GPS3.csv');
% Data_Array = table2array(Data_Table);
% Data_trimmed = Data_Array;

% Data_trimmed = [linspace(0, 0, 101)', linspace(0, 100, 101)', linspace(0, 0, 101)']*0.00001;

figure(2)
axis equal
plot(Data_trimmed(:,2),Data_trimmed(:,1))

Data_mean = mean(Data_trimmed);
Data_flat = lla2flat(Data_trimmed, Data_mean(1:2), 0, 0, 'WGS84');
Data_flat = [interp(Data_flat(:,1),interpulation), interp(Data_flat(:,2),interpulation), interp(Data_flat(:,3),interpulation)];

%% Track edge or individual cones (User variable modify and comment out as required))

Data_flat_blue = [NaN,NaN,NaN];
Data_flat_yelw = [NaN,NaN,NaN];
Data_flat_big = [NaN,NaN,NaN];
Cone_blue3 = [NaN,NaN];
Cone_yelw3 = [NaN,NaN];
Cone_big3 = [NaN,NaN];
Cone_orng3 = [NaN,NaN];

% Data_Table_blue = readtable('QUT GPS3.csv');
% Data_Array_blue = table2array(Data_Table_blue );
% Data_trimmed_blue = Data_Array_blue ;
% Data_mean_blue = mean(Data_trimmed_blue);
% Data_flat_blue = [Data_flat_blue; round(lla2flat(Data_trimmed_blue, Data_mean_blue (1:2), 0, 0, 'WGS84'))/2];
% Data_flat_blue = unique(Data_flat_blue, 'rows','stable');
% % Data_flat_blue = [interp(Data_flat_blue(:,2),interpulation), interp(Data_flat_blue(:,1),interpulation), interp(Data_flat_blue(:,3),interpulation)];
% for k = 1:size(Data_flat_blue,1)
%     Cone_cost = (sum((Data_flat_blue-Data_flat_blue(k,:)).^2, 2)).^0.5;
%     Cone_logi = logical(Cone_cost < (Cone_sep*0.5));
%     Data_flat_blue(Cone_logi,:) = ones(sum(Cone_logi),1)*mean(Data_flat_blue(Cone_logi,:));
% end
% for k = 1:size(Data_flat_blue,1)
%     Cone_cost = (sum((Data_flat_blue-Data_flat_blue(k,:)).^2, 2)).^0.5;
%     Cone_logi = logical(Cone_cost < (Cone_sep*0.5));
%     Data_flat_blue(Cone_logi,:) = ones(sum(Cone_logi),1)*mean(Data_flat_blue(Cone_logi,:));
% end
% Data_flat_blue = unique(Data_flat_blue, 'rows','stable');
% 
% Data_Table_yelw = readtable('QUT GPS3.csv');
% Data_Array_yelw = table2array(Data_Table_yelw);
% Data_trimmed_yelw = Data_Array_yelw;
% Data_mean_yelw = mean(Data_trimmed_yelw);
% Data_flat_yelw = [Data_flat_yelw; round(lla2flat(Data_trimmed_yelw, Data_mean_yelw(1:2), 0, 0, 'WGS84'))/2];
% Data_flat_yelw = unique(Data_flat_yelw, 'rows','stable');
% % Data_flat_yelw = [interp(Data_flat_yelw(:,2),interpulation), interp(Data_flat_yelw(:,1),interpulation), interp(Data_flat_yelw(:,3),interpulation)];
% for k = 1:size(Data_flat_yelw,1)
%     Cone_cost = (sum((Data_flat_yelw-Data_flat_yelw(k,:)).^2, 2)).^0.5;
%     Cone_logi = logical(Cone_cost < (Cone_sep*0.5));
%     Data_flat_yelw(Cone_logi,:) = ones(sum(Cone_logi),1)*mean(Data_flat_yelw(Cone_logi,:));
% end
% for k = 1:size(Data_flat_yelw,1)
%     Cone_cost = (sum((Data_flat_yelw-Data_flat_yelw(k,:)).^2, 2)).^0.5;
%     Cone_logi = logical(Cone_cost < (Cone_sep*0.5));
%     Data_flat_yelw(Cone_logi,:) = ones(sum(Cone_logi),1)*mean(Data_flat_yelw(Cone_logi,:));
% end
% Data_flat_yelw = unique(Data_flat_yelw, 'rows','stable');
% 
% Data_Table_big = readtable('QUT GPS3.csv');
% Data_Array_big = table2array(Data_Table_big);
% Data_trimmed_big = Data_Array_big;
% Data_mean_big = mean(Data_trimmed_big);
% Data_flat_big = [Data_flat_big; round(lla2flat(Data_trimmed_big, Data_mean_big(1:2), 0, 0, 'WGS84'))/2];
% Data_flat_big = unique(Data_flat_big, 'rows','stable');
% Data_flat_big = [interp(Data_flat_big(:,2),interpulation), interp(Data_flat_big(:,1),interpulation), interp(Data_flat_big(:,3),interpulation)];
% for k = 1:size(Data_flat_big,1)
%     Cone_cost = (sum((Data_flat_big-Data_flat_big(k,:)).^2, 2)).^0.5;
%     Cone_logi = logical(Cone_cost < (Cone_sep*0.5));
%     Data_flat_big(Cone_logi,:) = ones(sum(Cone_logi),1)*mean(Data_flat_big(Cone_logi,:));
% end
% for k = 1:size(Data_flat_big,1)
%     Cone_cost = (sum((Data_flat_big-Data_flat_big(k,:)).^2, 2)).^0.5;
%     Cone_logi = logical(Cone_cost < (Cone_sep*0.5));
%     Data_flat_big(Cone_logi,:) = ones(sum(Cone_logi),1)*mean(Data_flat_big(Cone_logi,:));
% end
% Data_flat_big = unique(Data_flat_big, 'rows','stable');
% 
% Data_Table_orng = readtable('QUT GPS3.csv');
% Data_Array_orng = table2array(Data_Table_orng);
% Data_trimmed_orng = Data_Array_orng;
% Data_mean_orng = mean(Data_trimmed_orng);
% Data_flat_orng = [Data_flat_orng; round(lla2flat(Data_trimmed_orng, Data_mean_orng(1:2), 0, 0, 'WGS84'))/2];
% Data_flat_orng = unique(Data_flat_orng, 'rows','stable');
% Data_flat_orng = [interp(Data_flat_orng(:,2),interpulation), interp(Data_flat_orng(:,1),interpulation), interp(Data_flat_orng(:,3),interpulation)];
% for k = 1:size(Data_flat_orng,1)
%     Cone_cost = (sum((Data_flat_orng-Data_flat_orng(k,:)).^2, 2)).^0.5;
%     Cone_logi = logical(Cone_cost < (Cone_sep*0.5));
%     Data_flat_orng(Cone_logi,:) = ones(sum(Cone_logi),1)*mean(Data_flat_orng(Cone_logi,:));
% end
% for k = 1:size(Data_flat_orng,1)
%     Cone_cost = (sum((Data_flat_orng-Data_flat_orng(k,:)).^2, 2)).^0.5;
%     Cone_logi = logical(Cone_cost < (Cone_sep*0.5));
%     Data_flat_orng(Cone_logi,:) = ones(sum(Cone_logi),1)*mean(Data_flat_orng(Cone_logi,:));
% end
% Data_flat_orng = unique(Data_flat_orng, 'rows','stable');

%% Processing Geo to Flat

if acceleration_test == true
    Data_flat = [linspace(0, 0, 1001)', linspace(0, straight_length + break_length - 1, 1001)', linspace(0, 0, 1001)'];
end
if skidpad_test == true
    Data_flat = [linspace(0, 0, 1001)', linspace(0, 40, 1001)', linspace(0, 0, 1001)'];
end

figure(3)
axis equal
hold on
plot(Data_flat(:,2),Data_flat(:,1))

Data_map = (round(Data_flat/resolution))*resolution;
Data_map2 = unique(Data_map, 'rows','stable');

figure(4)
axis equal
hold on
plot(Data_map2(:,2),Data_map2(:,1), 'k.')

%% Convert Path to Occupancy Grid

Data_grid = zeros(round((range(Data_map2(:,1:2))+[50,50])/resolution));
Data_grid(end+(100/resolution),end+(100/resolution)) = 0;
offset_x = (50-min(Data_map2(1:end,1)))/resolution;
offset_y = (50-min(Data_map2(1:end,2)))/resolution;
for k = 1:size(Data_map2,1)
    Data_grid(round((Data_map2(k,1))/resolution+offset_x),round((Data_map2(k,2))/resolution+offset_y)) = 1;
end
Data_grid = bwmorph(Data_grid, 'skel');

%% Create Track 

if (acceleration_test || skidpad_test) == true 
    filter_shape1 = ones(2+round(Path_width/resolution));
    filter_shape2 = ones(round(Path_width/resolution),1+round(Path_width/resolution));
    clockwise = true;
else
    filter_shape1 = fspecial('disk', (0.5*Path_width/resolution)+1) > 0;
    filter_shape2 = fspecial('disk', (0.5*Path_width/resolution)) > 0;
end

Data_wide = imdilate(Data_grid, filter_shape1);
Data_road = imdilate(Data_grid, filter_shape2);
Data_cone = Data_wide - Data_road;
Data_cone = bwmorph(Data_cone, 'skel');

figure(5)
axis equal
imshow(flip(Data_cone))

%% Convert Occupancy Grid to Points (unordered as path)

[Cone_map(:,1),Cone_map(:,2)] = find(Data_cone);

Cone_map(:,1) = (Cone_map(:,1))-offset_x;
Cone_map(:,2) = (Cone_map(:,2))-offset_y;

Cone_map_x = Cone_map(1:end,2)*resolution;
Cone_map_y = Cone_map(1:end,1)*resolution;

%% Seperate inside of track from outside

index = inpolygon(Cone_map_x, Cone_map_y, [Data_map(1:end,2);Data_map(1,2)], [Data_map(1:end,1);Data_map(1,1)]);

if ~any(index)
    index = Cone_map_y < 0;
end

index = index == clockwise;

Cone_blue1 = [Cone_map_x(~index), Cone_map_y(~index)];
Cone_yelw1 = [Cone_map_x(index), Cone_map_y(index)];
Cone_orng1 = [NaN, NaN];

if acceleration_test == true
    Cone_orng1 = Cone_blue1(~logical(Cone_blue1(:,1) < straight_length), :);
    Cone_blue1 = Cone_blue1(logical(Cone_blue1(:,1) < straight_length), :);
    Cone_orng1 = [Cone_orng1; Cone_yelw1(~logical(Cone_yelw1(:,1) < straight_length), :)];
    Cone_yelw1 = Cone_yelw1(logical(Cone_yelw1(:,1) < straight_length), :);
end
if skidpad_test == true
    Cone_orng1 = Cone_blue1(logical((Cone_blue1(:,1) < entrance_length) + (Cone_blue1(:,1) > exit_length)), :);
    Cone_blue1 = Cone_blue1(~logical((Cone_blue1(:,1) < entrance_length) + (Cone_blue1(:,1) > exit_length)), :);
    Cone_orng1 = [Cone_orng1; Cone_yelw1(logical((Cone_yelw1(:,1) < entrance_length) + (Cone_yelw1(:,1) > exit_length)), :)];
    Cone_yelw1 = Cone_yelw1(~logical((Cone_yelw1(:,1) < entrance_length) + (Cone_yelw1(:,1) > exit_length)), :);
    Cone_blue1 = [Cone_blue1(1,:);Cone_blue1(end,:)];
    Cone_yelw1 = [Cone_yelw1(1,:);Cone_yelw1(end,:)];
end

%% Track edges to cone locations

% Start_point = [mean([Data_flat(1,2),Data_flat(end,2)],2), mean([Data_flat(1,1),Data_flat(end,1)],2)];%Data_flat(1,1:2);%
Start_point = Data_flat(1,1:2);
End_point = Data_flat(end,1:2);

Start_line = nan(2,2);
End_line = nan(2,2);

n = 1; % set the start big cone location
Cone_cost_old = Path_width + Cone_sep;
for k = 1:size(Cone_blue1,1)
    Cone_cost = abs((sum((Start_point-Cone_blue1(k,:)).^2, 2)).^0.5);
    if Cone_cost < Cone_cost_old
        n = k;
        Cone_cost_old = Cone_cost;
    end
end
Start_line(1,:) = Cone_blue1(n,:);

n = 1; % set the start big cone location
Cone_cost_old = Path_width + Cone_sep;
for k = 1:size(Cone_yelw1,1)
    Cone_cost = abs((sum((Start_point-Cone_yelw1(k,:)).^2, 2)).^0.5);
    if Cone_cost < Cone_cost_old
        n = k;
        Cone_cost_old = Cone_cost;
    end
end
Start_line(2,:) = Cone_yelw1(n,:);

n = 1; % set the end big cone location
Cone_cost_old = Path_width + Cone_sep;
for k = 1:size(Cone_blue1,1)
    Cone_cost = abs((sum((End_point-Cone_blue1(k,:)).^2, 2)).^0.5);
    if Cone_cost < Cone_cost_old
        n = k;
        Cone_cost_old = Cone_cost;
    end
end
End_line(1,:) = Cone_blue1(n,:);

n = 1; % set the end big cone location
Cone_cost_old = Path_width + Cone_sep;
for k = 1:size(Cone_yelw1,1)
    Cone_cost = abs((sum((End_point-Cone_yelw1(k,:)).^2, 2)).^0.5);
    if Cone_cost < Cone_cost_old
        n = k;
        Cone_cost_old = Cone_cost;
    end
end
End_line(2,:) = Cone_yelw1(n,:);

if skidpad_test == false
    for k = 1:size(Cone_blue1,1)
        Cone_cost = (sum((Cone_blue1-Cone_blue1(k,:)).^2, 2)).^0.5;
        Cone_logi = logical(Cone_cost < (Cone_sep*0.5));
        Cone_blue1(Cone_logi,:) = ones(sum(Cone_logi),1)*mean(Cone_blue1(Cone_logi,:));
    end
    Cone_blue2 = Cone_blue1;
    for k = 1:size(Cone_blue2,1)
        Cone_cost = (sum((Cone_blue2-Cone_blue2(k,:)).^2, 2)).^0.5;
        Cone_logi = logical(Cone_cost < (Cone_sep*0.5));
        Cone_blue2(Cone_logi,:) = ones(sum(Cone_logi),1)*mean(Cone_blue2(Cone_logi,:));
    end
    Cone_blue3 = [Cone_blue3; unique(Cone_blue2,'rows','stable')];
    for k = 1:size(Cone_blue3,1)
        if Cone_blue3(k,1) == Cone_blue3(k,2)
           Cone_blue3(k,:) = [NaN,NaN];
        end
    end
    Cone_blue3(any(isnan(Cone_blue3), 2), :) = [];

    n = 1;
    Cone_cost_old = Path_width + Cone_sep;
    for k = 1:size(Cone_blue3,1)
        Cone_cost = abs((sum((Start_point-Cone_blue3(k,:)).^2, 2)).^0.5);
        if Cone_cost < Cone_cost_old
            n = k;
            Cone_cost_old = Cone_cost;
        end
    end
    Cone_blue3(n,:) = [NaN,NaN];

    Cone_blue3 = [Cone_blue3;Data_flat_blue(:,1:2)];

    Cone_blue3(any(isnan(Cone_blue3), 2), :) = [];

    for k = 1:size(Cone_yelw1,1)
        Cone_cost = (sum((Cone_yelw1-Cone_yelw1(k,:)).^2, 2)).^0.5;
        Cone_logi = logical(Cone_cost < (Cone_sep*0.5));
        Cone_yelw1(Cone_logi,:) = ones(sum(Cone_logi),1)*mean(Cone_yelw1(Cone_logi,:));
    end
    Cone_yelw2 = Cone_yelw1;
    for k = 1:size(Cone_yelw2,1)
        Cone_cost = (sum((Cone_yelw2-Cone_yelw2(k,:)).^2, 2)).^0.5;
        Cone_logi = logical(Cone_cost < (Cone_sep*0.5));
        Cone_yelw2(Cone_logi,:) = ones(sum(Cone_logi),1)*mean(Cone_yelw2(Cone_logi,:));
    end
    Cone_yelw3 = [Cone_yelw3; unique(Cone_yelw2,'rows','stable')];
    for k = 1:size(Cone_yelw3,1)
        if Cone_yelw3(k,1) == Cone_yelw3(k,2)
           Cone_yelw3(k,:) = [NaN,NaN];
        end
    end
    Cone_yelw3(any(isnan(Cone_yelw3), 2), :) = [];

    n = 1;
    Cone_cost_old = Path_width + Cone_sep;
    for k = 1:size(Cone_yelw3,1)
        Cone_cost = abs((sum((Start_point-Cone_yelw3(k,:)).^2, 2)).^0.5);
        if Cone_cost < Cone_cost_old
            n = k;
            Cone_cost_old = Cone_cost;
        end
    end
    Cone_yelw3(n,:) = [NaN,NaN];
    Start_point = mean([Start_line(1,:); Start_line(2,:)], 1);
end

if skidpad_test == true
    Cone_blue3 = Cone_blue1;
    Cone_yelw3 = Cone_yelw1;
end
    
Cone_orng1 = Cone_orng1.*[1,Cone_sep];
for k = 1:size(Cone_orng1,1)
    Cone_cost = (sum((Cone_orng1-Cone_orng1(k,:)).^2, 2)).^0.5;
    Cone_logi = logical(Cone_cost < (Cone_sep*0.5));
    Cone_orng1(Cone_logi,:) = ones(sum(Cone_logi),1)*mean(Cone_orng1(Cone_logi,:));
end
Cone_orng2 = Cone_orng1;
for k = 1:size(Cone_orng2,1)
    Cone_cost = (sum((Cone_orng2-Cone_orng2(k,:)).^2, 2)).^0.5;
    Cone_logi = logical(Cone_cost < (Cone_sep*0.5));
    Cone_orng2(Cone_logi,:) = ones(sum(Cone_logi),1)*mean(Cone_orng2(Cone_logi,:));
end
Cone_orng3 = [Cone_orng3; unique(Cone_orng2,'rows','stable')];
for k = 1:size(Cone_orng3,1)
    if Cone_orng3(k,1) == Cone_orng3(k,2)
       Cone_orng3(k,:) = [NaN,NaN];
    end
end
Cone_orng3 = Cone_orng3./[1,Cone_sep];
Cone_orng3(any(isnan(Cone_orng3), 2), :) = [];

Cone_yelw3 = [Cone_yelw3;Data_flat_yelw(:,1:2)];
Cone_yelw3(any(isnan(Cone_yelw3), 2), :) = [];

Start_line = [Start_line;Data_flat_big(:,1:2)];
Start_line(any(isnan(Start_line), 2), :) = [];

%% Translate and Rotate

% Data_flat(:,1:2) =  Data_flat(:,1:2) - fliplr(Start_point);
% Cone_blue3 = Cone_blue3 - Start_point;
% Cone_yelw3 = Cone_yelw3 - Start_point;
% Cone_big3 = Cone_big3 - Start_point;
% Cone_orng3 = Cone_orng3 - Start_point;
% End_line = End_line - Start_point;
% End_point = End_point - Start_point;
% Start_line = Start_line - Start_point;
% Start_point = Start_point - Start_point;
% 
% track_rotation = 180 - track_rotation - atan2d(Start_line(2,1) - Start_line(1,1), Start_line(2,2) - Start_line(1,2));
% Rotation_matrix = [cosd(track_rotation), -sind(track_rotation); sind(track_rotation), cosd(track_rotation)];
% 
% Data_flat(:,1:2) = (Rotation_matrix * Data_flat(:,1:2)')';
% Cone_blue3 = fliplr((Rotation_matrix * flip(Cone_blue3'))');
% Cone_yelw3 = fliplr((Rotation_matrix * flip(Cone_yelw3'))');
% Cone_big3 = fliplr((Rotation_matrix * flip(Cone_big3'))');
% Cone_orng3 = fliplr((Rotation_matrix * flip(Cone_orng3'))');
% Start_line = fliplr((Rotation_matrix * flip(Start_line'))');
% Start_point = fliplr((Rotation_matrix * flip(Start_point'))');
% End_line = fliplr((Rotation_matrix * flip(End_line'))');
% End_point = fliplr((Rotation_matrix * flip(End_point'))');
% 
% Data_flat(:,1:2) =  fliplr(track_offset) + Data_flat(:,1:2);
% Cone_blue3 = track_offset + Cone_blue3;
% Cone_yelw3 = track_offset + Cone_yelw3;
% Cone_big3 = track_offset + Cone_big3;
% Cone_orng3 = track_offset + Cone_orng3;
% Start_line = track_offset + Start_line;
% Start_point = track_offset + Start_point;
% End_line = track_offset + End_line;
% End_point = track_offset + End_point;

%% Finalise

% if skidpad_test == false
%     if acceleration_test == true
%         End_line = Start_line + [straight_length, 0];
%         End_point = Start_point + [straight_length, 0];
%     else
%         End_line = Start_line;
%         End_point = Start_point;
%     end
% end
% End_line(any(isnan(End_line), 2), :) = [];
% 
% if skidpad_test == true
%     Start_line = (Start_line + End_line)/2;
%     Start_point = (Start_point + End_point)/2;
%     End_line = Start_line;
%     End_point = Start_point;
% end

Cone_big3 = [Cone_big3 ;[Start_line; End_line] + [big_cone_sep/2, 0]; [Start_line; End_line] - [big_cone_sep/2, 0]];
Cone_big3(any(isnan(Cone_big3), 2), :) = [];

figure(1)
axis equal
hold on
set(gca,'Color',[0.25, 0.25, 0.25])
plot(Data_flat(1:end,2),Data_flat(1:end,1))
plot(Cone_blue3(1:end,1), Cone_blue3(1:end,2), 'b.')
plot(Cone_yelw3(1:end,1), Cone_yelw3(1:end,2), 'y.')
plot(Cone_big3(1:end,1), Cone_big3(1:end,2), 'o', 'Color', [1,0.5,0])
plot(Cone_orng3(1:end,1), Cone_orng3(1:end,2), '.', 'Color', [1,0.5,0])
plot(Start_line(1:end,1), Start_line(1:end,2), 'r')
plot(End_line(1:end,1), End_line(1:end,2), 'r')
% plot(Start_point(1:end,1), Start_point(1:end,2), 'k.')
% plot(End_point(1:end,1), End_point(1:end,2), 'k.')

%% Make files

str = sprintf('%s.sdf', File_Name);
sdf = fopen(str, 'w');
fprintf( sdf, '<?xml version=''1.0''?>\n<sdf version=''1.6''>\n  <model name=''%s''>\n', Track_Name);
for k = 1:size(Cone_blue3,1)
    Cone_rot = pi * rand();
    fprintf(sdf, '    <link name=''link_%u''>\n', k);
    fprintf(sdf, '      <pose frame=''''>%f %f 0.070966 0 0 0</pose>\n', Cone_blue3(k,1), Cone_blue3(k,2));%, Cone_rot);
    fprintf(sdf, '      <collision name="collision">\n');
    fprintf(sdf, '        <geometry>\n');
    fprintf(sdf, '          <mesh>\n');
    fprintf(sdf, '            <uri>model:%s</uri>\n', Cone_Model_blue);
    fprintf(sdf, '          </mesh>\n');
    fprintf(sdf, '        </geometry>\n');
    fprintf(sdf, '      </collision>\n');
    fprintf(sdf, '      <visual name="visual">\n');
    fprintf(sdf, '        <geometry>\n');
    fprintf(sdf, '          <mesh>\n');
    fprintf(sdf, '            <uri>model:%s</uri>\n', Cone_Model_blue);
    fprintf(sdf, '          </mesh>\n');
    fprintf(sdf, '        </geometry>\n');
    fprintf(sdf, '      </visual>\n');
    fprintf(sdf, '    </link>\n');
end
for m = 1:size(Cone_yelw3,1)
    Cone_rot = pi * rand();
    fprintf(sdf, '    <link name=''link_%u''>\n', m+k);
    fprintf(sdf, '      <pose frame=''''>%f %f 0.070966 0 0 0</pose>\n', Cone_yelw3(m,1), Cone_yelw3(m,2));%, Cone_rot);
    fprintf(sdf, '      <collision name="collision">\n');
    fprintf(sdf, '        <geometry>\n');
    fprintf(sdf, '          <mesh>\n');
    fprintf(sdf, '            <uri>model:%s</uri>\n', Cone_Model_yelw);
    fprintf(sdf, '          </mesh>\n');
    fprintf(sdf, '        </geometry>\n');
    fprintf(sdf, '      </collision>\n');
    fprintf(sdf, '      <visual name="visual">\n');
    fprintf(sdf, '        <geometry>\n');
    fprintf(sdf, '          <mesh>\n');
    fprintf(sdf, '            <uri>model:%s</uri>\n', Cone_Model_yelw);
    fprintf(sdf, '          </mesh>\n');
    fprintf(sdf, '        </geometry>\n');
    fprintf(sdf, '      </visual>\n');
    fprintf(sdf, '    </link>\n');
end
for x = 1:size(Cone_big3,1)
    fprintf(sdf, '    <link name=''link_%u''>\n', x+m+k);
    fprintf(sdf, '      <pose frame=''''>%f %f 0.070966 0 0 0</pose>\n', Cone_big3(x,1), Cone_big3(x,2));
    fprintf(sdf, '      <collision name="collision">\n');
    fprintf(sdf, '        <geometry>\n');
    fprintf(sdf, '          <mesh>\n');
    fprintf(sdf, '            <uri>model:%s</uri>\n', Cone_Model_big);
    fprintf(sdf, '          </mesh>\n');
    fprintf(sdf, '        </geometry>\n');
    fprintf(sdf, '      </collision>\n');
    fprintf(sdf, '      <visual name="visual">\n');
    fprintf(sdf, '        <geometry>\n');
    fprintf(sdf, '          <mesh>\n');
    fprintf(sdf, '            <uri>model:%s</uri>\n', Cone_Model_big);
    fprintf(sdf, '          </mesh>\n');
    fprintf(sdf, '        </geometry>\n');
    fprintf(sdf, '      </visual>\n');
    fprintf(sdf, '    </link>\n');
end
for z = 1:size(Cone_orng3,1)
    fprintf(sdf, '    <link name=''link_%u''>\n', z+x+m+k);
    fprintf(sdf, '      <pose frame=''''>%f %f 0.070966 0 0 0</pose>\n', Cone_orng3(z,1), Cone_orng3(z,2));
    fprintf(sdf, '      <collision name="collision">\n');
    fprintf(sdf, '        <geometry>\n');
    fprintf(sdf, '          <mesh>\n');
    fprintf(sdf, '            <uri>model:%s</uri>\n', Cone_Model_orng);
    fprintf(sdf, '          </mesh>\n');
    fprintf(sdf, '        </geometry>\n');
    fprintf(sdf, '      </collision>\n');
    fprintf(sdf, '      <visual name="visual">\n');
    fprintf(sdf, '        <geometry>\n');
    fprintf(sdf, '          <mesh>\n');
    fprintf(sdf, '            <uri>model:%s</uri>\n', Cone_Model_orng);
    fprintf(sdf, '          </mesh>\n');
    fprintf(sdf, '        </geometry>\n');
    fprintf(sdf, '      </visual>\n');
    fprintf(sdf, '    </link>\n');
end
fprintf( sdf, '    <static>0</static>\n    <allow_auto_disable>1</allow_auto_disable>\n  </model>\n</sdf>');
fclose(sdf);

str = sprintf('%s.config', File_Name);
config = fopen(str, 'w');
fprintf(config, '<?xml version="1.0" ?>\n');
fprintf(config, '<model>\n');
fprintf(config, '    <name>%s</name>\n', Track_Name);
fprintf(config, '    <version>1.0</version>\n');
fprintf(config, '    <sdf version="1.6">%s.sdf</sdf>\n', File_Name);
fprintf(config, '    <author>\n');
fprintf(config, '        <name>%s</name>\n', Author_Name);% change for author
fprintf(config, '        <email>%s</email>\n', Author_Email);% change for author
fprintf(config, '    </author>\n');
fprintf(config, '    <description>This config file has been automaticly generated</description>\n');
fprintf(config, '</model>');
fclose(config);