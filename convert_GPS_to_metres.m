%% This code is for QUT Motorsport. It was written using Matlab R2018b durring the 2020 design camp by the Autonomus Driving Team.
clc, clear variables, close all

Author_Name = 'Nathan Nahrung';%Listed author should be point of contact for troubleshooting
Author_Email = 'nnahrung@gmail.com';%Listed authors contact email

%% User variables

Track_Name = 'lakeside track';%Name of track relivant to GPS data or standard track type (e.g. acceleration or skidpad)
File_Name = 'model';%used for .sdf & .config
track_width_x = 3;%greater than (error aprox = resolution) - Note: FSG Rule Competition Handbook 2020 DE6.3 - minimum 3 m
track_width_y = 3;%greater than (error aprox = resolution) - Note: FSG Rule Competition Handbook 2020 DE6.3 - minimum 3 m
Cone_sep = 5;%less than (error aprox = resolution) - Note: FSG Rule Competition Handbook 2020 DE6.3 - cone sep should be less than 5 m and reduced at corners
tall_cone_sep = 0.6;%distance between cone pairs at start and finish line
resolution = 0.1;%grid unit in metres
interpulation = 4;%used to close gaps in GPS data
track_offset = [0, 0];%meters xy offset of track (by default center of start is [0, 0])
track_rotation = 0;%Starting direction is positive x, change this bearing anticlockwise with degrees

% Below can be set true or false to activate or deactive features
track_direction_clockwise = false;%true;% 

% Below can be set true or false to activate or deactive features
acceleration_test = false;%true;%
acceleration_straight_length = 75;%acceleration_test variable (acceleration length)
acceleration_break_length = 100;%acceleration_test variable (breaking length)

% Below can be set true or false to activate or deactive features
skidpad_test = false;%true;%
skidpad_entrance_length = 7.4;%skidpad_test variable end of entrance from [0,0]
skidpad_mid_length = 15;%skidpad_test variable mid point of circles from [0,0] 
skidpad_exit_length = 22.6;%skidpad_test variable start of exit from [0,0]
skidpad_inner_diameter = 15.2;%Diameter of inside cone circle
skidpad_outer_diameter = 21.2;%Diameter of outside cone cicle
skidpad_circle_sep = 18.25;%Distance between turning circles circles

% Below can be set true or false to activate or deactive features
Add_extra_blue_cones = false;%true;%
Add_extra_yelw_cones = false;%true;%
Add_extra_tall_cones = false;%true;%
Add_extra_orng_cones = false;%true;%

%% Add used model mesh files with locations (as below)

Cone_Model_blue = '//qut_description/meshes/cone_blue.dae';
Cone_Model_yelw = '//qut_description/meshes/cone_yellow.dae';
Cone_Model_tall = '//qut_description/meshes/cone_tall.dae';
Cone_Model_orng = '//qut_description/meshes/cone.dae';

%% Track center GPS Data (User variable, copy, modify and comment out as required)
% NOTE: GPS data should be '.csv', column 1 latitude, column 2 longitude, column 3 altitude (not used, can be 0)
% NOTE: The App "Ultra GPS Logger" can be used to collect GPS data and export as '.csv'
% NOTE: Ensure the data is trimmed to a clean lap. Start point will be at the start of the GPS data.

Data_Table = readtable('lakeside GPS.csv');
Data_Array = table2array(Data_Table);
Data_trimmed = Data_Array(5650:24780,1:3);

% Data_Table = readtable('winton GPS.csv');
% Data_Array = table2array(Data_Table);
% Data_trimmed = Data_Array(12160:13640,1:3);

% Data_Table = readtable('QUT GPS3.csv');
% Data_Array = table2array(Data_Table);
% Data_trimmed = Data_Array;

figure(2)
axis equal
plot(Data_trimmed(:,2),Data_trimmed(:,1))

%% Track edge or individual cones (User variable modify and comment out as required))

Data_flat_blue = [NaN,NaN,NaN];
Data_flat_yelw = [NaN,NaN,NaN];
Data_flat_tall = [NaN,NaN,NaN];
Data_flat_orng = [NaN,NaN,NaN];
Cone_blue3 = [NaN,NaN];
Cone_yelw3 = [NaN,NaN];
Cone_tall3 = [NaN,NaN];
Cone_orng3 = [NaN,NaN];

if Add_extra_blue_cones == true
    Data_Table_blue = readtable('QUT GPS3.csv');
    Data_Array_blue = table2array(Data_Table_blue );
    Data_trimmed_blue = Data_Array_blue ;
    Data_mean_blue = mean(Data_trimmed_blue);
    Data_flat_blue = [Data_flat_blue; round(lla2flat(Data_trimmed_blue, Data_mean_blue (1:2), 0, 0, 'WGS84'))/2];
    Data_flat_blue = unique(Data_flat_blue, 'rows','stable');% Data_flat_blue = [interp(Data_flat_blue(:,2),interpulation), interp(Data_flat_blue(:,1),interpulation), interp(Data_flat_blue(:,3),interpulation)];
    for k = 1:size(Data_flat_blue,1)
        Cone_cost = (sum((Data_flat_blue-Data_flat_blue(k,:)).^2, 2)).^0.5;
        Cone_logi = logical(Cone_cost < (Cone_sep*0.5));
        Data_flat_blue(Cone_logi,:) = ones(sum(Cone_logi),1)*mean(Data_flat_blue(Cone_logi,:));
    end
    Data_flat_blue = unique(Data_flat_blue, 'rows','stable');
end

if Add_extra_yelw_cones == true
    Data_Table_yelw = readtable('QUT GPS3.csv');
    Data_Array_yelw = table2array(Data_Table_yelw);
    Data_trimmed_yelw = Data_Array_yelw;
    Data_mean_yelw = mean(Data_trimmed_yelw);
    Data_flat_yelw = [Data_flat_yelw; round(lla2flat(Data_trimmed_yelw, Data_mean_yelw(1:2), 0, 0, 'WGS84'))/2];
    Data_flat_yelw = unique(Data_flat_yelw, 'rows','stable');% Data_flat_yelw = [interp(Data_flat_yelw(:,2),interpulation), interp(Data_flat_yelw(:,1),interpulation), interp(Data_flat_yelw(:,3),interpulation)];
    for k = 1:size(Data_flat_yelw,1)
        Cone_cost = (sum((Data_flat_yelw-Data_flat_yelw(k,:)).^2, 2)).^0.5;
        Cone_logi = logical(Cone_cost < (Cone_sep*0.5));
        Data_flat_yelw(Cone_logi,:) = ones(sum(Cone_logi),1)*mean(Data_flat_yelw(Cone_logi,:));
    end
    Data_flat_yelw = unique(Data_flat_yelw, 'rows','stable');
end

if Add_extra_tall_cones == true
    Data_Table_tall = readtable('QUT GPS3.csv');
    Data_Array_tall = table2array(Data_Table_tall);
    Data_trimmed_tall = Data_Array_tall;
    Data_mean_tall = mean(Data_trimmed_tall);
    Data_flat_tall = [Data_flat_tall; round(lla2flat(Data_trimmed_tall, Data_mean_tall(1:2), 0, 0, 'WGS84'))/2];
    Data_flat_tall = unique(Data_flat_tall, 'rows','stable');% Data_flat_tall = [interp(Data_flat_tall(:,2),interpulation), interp(Data_flat_tall(:,1),interpulation), interp(Data_flat_tall(:,3),interpulation)];
    for k = 1:size(Data_flat_tall,1)
        Cone_cost = (sum((Data_flat_tall-Data_flat_tall(k,:)).^2, 2)).^0.5;
        Cone_logi = logical(Cone_cost < (Cone_sep*0.5));
        Data_flat_tall(Cone_logi,:) = ones(sum(Cone_logi),1)*mean(Data_flat_tall(Cone_logi,:));
    end
    Data_flat_tall = unique(Data_flat_tall, 'rows','stable');
end

if Add_extra_orng_cones == true
    Data_Table_orng = readtable('QUT GPS3.csv');
    Data_Array_orng = table2array(Data_Table_orng);
    Data_trimmed_orng = Data_Array_orng;
    Data_mean_orng = mean(Data_trimmed_orng);
    Data_flat_orng = [Data_flat_orng; round(lla2flat(Data_trimmed_orng, Data_mean_orng(1:2), 0, 0, 'WGS84'))/2];
    Data_flat_orng = unique(Data_flat_orng, 'rows','stable');% Data_flat_orng = [interp(Data_flat_orng(:,2),interpulation), interp(Data_flat_orng(:,1),interpulation), interp(Data_flat_orng(:,3),interpulation)];
    for k = 1:size(Data_flat_orng,1)
        Cone_cost = (sum((Data_flat_orng-Data_flat_orng(k,:)).^2, 2)).^0.5;
        Cone_logi = logical(Cone_cost < (Cone_sep*0.5));
        Data_flat_orng(Cone_logi,:) = ones(sum(Cone_logi),1)*mean(Data_flat_orng(Cone_logi,:));
    end
    Data_flat_orng = unique(Data_flat_orng, 'rows','stable');
end

%% Processing Geo to Flat (Convert GPS to Meters)

Data_mean = mean(Data_trimmed);%Find approximate center of center track (used for map projection)
Data_flat = lla2flat(Data_trimmed, Data_mean(1:2), 0, 0, 'WGS84');%Converts GPS to meters by planar projection around center of the track
Data_flat = [interp(Data_flat(:,1),interpulation), interp(Data_flat(:,2),interpulation), interp(Data_flat(:,3),interpulation)];%interpulation helps smothly fill GPS gaps

if acceleration_test == true%Overrides the loaded GPS map if true
    Data_flat = [linspace(0, 0, 1001)', linspace(0, acceleration_straight_length + acceleration_break_length - 1, 1001)', linspace(0, 0, 1001)'];
end
if skidpad_test == true%Overrides the loaded GPS map if true
    Data_flat = [linspace(0, 0, 1001)', linspace(0, 40, 1001)', linspace(0, 0, 1001)'];
end

figure(3)
axis equal
hold on
plot(Data_flat(:,2),Data_flat(:,1))

Data_map = (round(Data_flat/resolution))*resolution;%Reduces resolution of map
Data_map2 = unique(Data_map, 'rows','stable');%Reduces size of map data

figure(4)
axis equal
hold on
plot(Data_map2(:,2),Data_map2(:,1), 'k.')

%% Convert Path to Occupancy Grid

Data_grid = zeros(round((range(Data_map2(:,1:2))+[50,50])/resolution));%Make blank occupancy grid
Data_grid(end+(100/resolution),end+(100/resolution)) = 0;%Increase size of blank occupancy grid
offset_x = (50-min(Data_map2(1:end,1)))/resolution;%Define offset x
offset_y = (50-min(Data_map2(1:end,2)))/resolution;%Define offset x
for k = 1:size(Data_map2,1)
    Data_grid(round((Data_map2(k,1))/resolution+offset_x),round((Data_map2(k,2))/resolution+offset_y)) = 1;
end
Data_grid = bwmorph(Data_grid, 'skel');%Use image morph orperation to reduce path to thin line

%% Create Track 

if (acceleration_test || skidpad_test) == true%Rectanglular shape used to detirmine width of some tracks
    filter_shape1 = ones(2+round(track_width_x/resolution));
    filter_shape2 = ones(round(track_width_x/resolution),1+round(track_width_x/resolution));
    track_direction_clockwise = true;
else%Eliptical shape used to detirmine width of typical tracks
    filter_shape = fspecial('disk', ((track_width_y+track_width_x)/resolution)) > 0;
    filter_shape1 = imresize(filter_shape,[(track_width_y/resolution)+2 (track_width_x/resolution)+2]);
    filter_shape2 = imresize(filter_shape,[(track_width_y/resolution) (track_width_x/resolution)]);
end

Data_wide = imdilate(Data_grid, filter_shape1);%Outer width of cones
Data_road = imdilate(Data_grid, filter_shape2);%Inner width of cones
Data_cone = Data_wide - Data_road;%Subtract Inner from Outer width of cones to creat 2 parralel lines
Data_cone = bwmorph(Data_cone, 'skel');%Use image morph orperation to reduce cones to two thin lines

figure(5)
axis equal
imshow(flip(Data_cone))

%% Convert Occupancy Grid to Points (unordered as path)

[Cone_map(:,1),Cone_map(:,2)] = find(Data_cone);%Convert Occupancy grid to vector of points

Cone_map(:,1) = (Cone_map(:,1))-offset_x;%Remove offset x
Cone_map(:,2) = (Cone_map(:,2))-offset_y;%Remove offset y

Cone_map_x = Cone_map(1:end,2)*resolution;%Remove scalling in x
Cone_map_y = Cone_map(1:end,1)*resolution;%Remove scalling in y

%% Seperate inside of track from outside

index = inpolygon(Cone_map_x, Cone_map_y, [Data_map(1:end,2);Data_map(1,2)], [Data_map(1:end,1);Data_map(1,1)]);%For a closed loop track find cones on inside of track

if ~any(index)%This step finds if no cones are inside of track and if so assumes the track it straight so uses and alternative method to seperate cone types
    index = Cone_map_y < 0;
end

index = index == track_direction_clockwise;%If anticlockwise the cone list will invert

Cone_blue1 = [Cone_map_x(~index), Cone_map_y(~index)];%Assigns cones as blue
Cone_yelw1 = [Cone_map_x(index), Cone_map_y(index)];%Assigns cones as yellow
Cone_orng1 = [NaN, NaN];%Assigns cones as orange

% Below steps used to fix cone colours for special tracks
if acceleration_test == true
    Cone_orng1 = Cone_blue1(~logical(Cone_blue1(:,1) < acceleration_straight_length), :);
    Cone_blue1 = Cone_blue1(logical(Cone_blue1(:,1) < acceleration_straight_length), :);
    Cone_orng1 = [Cone_orng1; Cone_yelw1(~logical(Cone_yelw1(:,1) < acceleration_straight_length), :)];
    Cone_yelw1 = Cone_yelw1(logical(Cone_yelw1(:,1) < acceleration_straight_length), :);
end
if skidpad_test == true
    Cone_orng1 = Cone_blue1(logical((Cone_blue1(:,1) < skidpad_entrance_length) + (Cone_blue1(:,1) > skidpad_exit_length)), :);
    Cone_blue1 = Cone_blue1(~logical((Cone_blue1(:,1) < skidpad_entrance_length) + (Cone_blue1(:,1) > skidpad_exit_length)), :);
    Cone_orng1 = [Cone_orng1; Cone_yelw1(logical((Cone_yelw1(:,1) < skidpad_entrance_length) + (Cone_yelw1(:,1) > skidpad_exit_length)), :)];
    Cone_yelw1 = Cone_yelw1(~logical((Cone_yelw1(:,1) < skidpad_entrance_length) + (Cone_yelw1(:,1) > skidpad_exit_length)), :);
    Cone_blue1 = [Cone_blue1(1,:);Cone_blue1(end,:)];
    Cone_yelw1 = [Cone_yelw1(1,:);Cone_yelw1(end,:)];
end

%% Track edges to cone locations

Start_point = fliplr(Data_flat(1,1:2));%Define start point
End_point = fliplr(Data_flat(end,1:2));%Define end point

Start_line = nan(2,2);
End_line = nan(2,2);

[~,index] = sort(sum((Cone_blue1(:,:)-End_point).^2, 2).^0.5);
Cone_blue1 = Cone_blue1(index,:);
End_line(1,:) = Cone_blue1(1,:);%Find point of blue cone side closest the end

[~,index] = sort(sum((Cone_yelw1(:,:)-End_point).^2, 2).^0.5);
Cone_yelw1 = Cone_yelw1(index,:);
End_line(2,:) = Cone_yelw1(1,:);%Find point of yellow cone side closest the end

[~,index] = sort(sum((Cone_blue1(:,:)-Start_point).^2, 2).^0.5);
Cone_blue1 = Cone_blue1(index,:);
Start_line(1,:) = Cone_blue1(1,:);%Find point of blue cone side closest the start

[~,index] = sort(sum((Cone_yelw1(:,:)-Start_point).^2, 2).^0.5);
Cone_yelw1 = Cone_yelw1(index,:);
Start_line(2,:) = Cone_yelw1(1,:);%Find point of yellow cone side closest the start

if skidpad_test == false%Not used for skidpad
    Cone_cost = (sum((Cone_blue1-Start_line(1,:)).^2, 2)).^0.5;
    Cone_logi = logical(Cone_cost < (Cone_sep*0.1));
    Cone_blue1(Cone_logi,:) = [];
    for k = 1:size(Cone_blue1,1)%These steps place the blue cones at the maximum spaceing
        Cone_cost = (sum((Cone_blue1-Cone_blue1(k,:)).^2, 2)).^0.5;
        Cone_logi = logical(Cone_cost < (Cone_sep*0.5));
        Cone_blue1(Cone_logi,:) = ones(sum(Cone_logi),1)*mean(Cone_blue1(Cone_logi,:));
    end
    Cone_blue3 = [Cone_blue3; unique(Cone_blue1,'rows','stable')];
    for k = 1:size(Cone_blue3,1)
        if Cone_blue3(k,1) == Cone_blue3(k,2)
           Cone_blue3(k,:) = [NaN,NaN];
        end
    end
    Cone_blue3(any(isnan(Cone_blue3), 2), :) = [];
%     n = 1;
%     Cone_cost_old = track_width_x + Cone_sep;
%     for k = 1:size(Cone_blue3,1)
%         Cone_cost = abs((sum((Start_point-Cone_blue3(k,:)).^2, 2)).^0.5);
%         if Cone_cost < Cone_cost_old
%             n = k;
%             Cone_cost_old = Cone_cost;
%         end
%     end
%     Cone_blue3(n,:) = [NaN,NaN];
    Cone_blue3 = [Cone_blue3;Data_flat_blue(:,1:2)];
    Cone_blue3(any(isnan(Cone_blue3), 2), :) = [];

    Cone_cost = (sum((Cone_yelw1-Start_line(2,:)).^2, 2)).^0.5;
    Cone_logi = logical(Cone_cost < (Cone_sep*0.1));
    Cone_yelw1(Cone_logi,:) = [];
    for k = 1:size(Cone_yelw1,1)%These steps place the yellow cones at the maximum spaceing
        Cone_cost = (sum((Cone_yelw1-Cone_yelw1(k,:)).^2, 2)).^0.5;
        Cone_logi = logical(Cone_cost < (Cone_sep*0.5));
        Cone_yelw1(Cone_logi,:) = ones(sum(Cone_logi),1)*mean(Cone_yelw1(Cone_logi,:));
    end
    Cone_yelw3 = [Cone_yelw3; unique(Cone_yelw1,'rows','stable')];
    for k = 1:size(Cone_yelw3,1)
        if Cone_yelw3(k,1) == Cone_yelw3(k,2)
           Cone_yelw3(k,:) = [NaN,NaN];
        end
    end
    Cone_yelw3(any(isnan(Cone_yelw3), 2), :) = [];
%     n = 1;
%     Cone_cost_old = track_width_x + Cone_sep;
%     for k = 1:size(Cone_yelw3,1)
%         Cone_cost = abs((sum((Start_point-Cone_yelw3(k,:)).^2, 2)).^0.5);
%         if Cone_cost < Cone_cost_old
%             n = k;
%             Cone_cost_old = Cone_cost;
%         end
%     end
%     Cone_yelw3(n,:) = [NaN,NaN];
    Cone_yelw3 = [Cone_yelw3;Data_flat_yelw(:,1:2)];
    Cone_yelw3(any(isnan(Cone_yelw3), 2), :) = [];
    
    Start_point = mean([Start_line(1,:); Start_line(2,:)], 1);
end

if skidpad_test == true%inverts the blue and yellow cones
    Cone_blue3 = Cone_blue1;
    Cone_yelw3 = Cone_yelw1;
end

Cone_orng1 = Cone_orng1.*[1,Cone_sep];
for k = 1:size(Cone_orng1,1)%These steps place the orange cones at the maximum spaceing
    Cone_cost = (sum((Cone_orng1-Cone_orng1(k,:)).^2, 2)).^0.5;
    Cone_logi = logical(Cone_cost < (Cone_sep*0.5));
    Cone_orng1(Cone_logi,:) = ones(sum(Cone_logi),1)*mean(Cone_orng1(Cone_logi,:));
end
Cone_orng3 = [Cone_orng3; unique(Cone_orng1,'rows','stable')];
for k = 1:size(Cone_orng3,1)
    if Cone_orng3(k,1) == Cone_orng3(k,2)
       Cone_orng3(k,:) = [NaN,NaN];
    end
end
Cone_orng3 = [Cone_orng3./[1,Cone_sep];Data_flat_orng(:,1:2)];
Cone_orng3(any(isnan(Cone_orng3), 2), :) = [];

if skidpad_test == true%For the skidpad this adds the cones on the inner and outer circles
    blue_limit = min(Cone_blue3(:,2));
    yelw_limit = max(Cone_yelw3(:,2));
    mid_limit = mean([blue_limit, yelw_limit]);
    inner_cone_number = linspace(0, 1, ceil(pi*skidpad_inner_diameter/Cone_sep)+1);
    inner_circle = [(sin((inner_cone_number*pi*2)))', (cos((inner_cone_number*pi*2)))']*skidpad_inner_diameter*0.5;
    outer_cone_number = linspace(0, 1, ceil(pi*skidpad_outer_diameter/Cone_sep)+1);
    outer_circle = [(sin((outer_cone_number*pi*2)))', (cos((outer_cone_number*pi*2)))']*skidpad_outer_diameter*0.5;
    Cone_yelw3 = [Cone_yelw3; [skidpad_mid_length, mid_limit-skidpad_circle_sep/2] + [outer_circle(1:end-1,1),-outer_circle(1:end-1,2)]];
    Cone_blue3 = [Cone_blue3; [skidpad_mid_length, mid_limit+skidpad_circle_sep/2] + [outer_circle(1:end-1,1),outer_circle(1:end-1,2)]];
    Cone_yelw3 = Cone_yelw3(Cone_yelw3(:,2) < (yelw_limit+resolution), :);
    Cone_blue3 = Cone_blue3(Cone_blue3(:,2) > (blue_limit-resolution), :);
    Cone_yelw3 = [Cone_yelw3; [skidpad_mid_length, mid_limit+skidpad_circle_sep/2] + inner_circle(1:end-1,:)];
    Cone_blue3 = [Cone_blue3; [skidpad_mid_length, mid_limit-skidpad_circle_sep/2] + inner_circle(1:end-1,:)];
end

% Start_line = [Start_line;Data_flat_tall(:,1:2)];
Start_line(any(isnan(Start_line), 2), :) = [];
End_line(any(isnan(Start_line), 2), :) = [];

if acceleration_test == true
    Cone_tall3 = [Start_line;Start_line;End_line;End_line;Data_flat_tall(:,1:2);Cone_tall3];
else
    Cone_tall3 = [Start_line;Start_line;Data_flat_tall(:,1:2);Cone_tall3];
end

%% Translate and Rotate
%NOTE: The purpose of this section is to rotate and translate the map to orintate it to the start line

Data_flat(:,1:2) =  Data_flat(:,1:2) - fliplr(Start_point);
Cone_blue3 = Cone_blue3 - Start_point;
Cone_yelw3 = Cone_yelw3 - Start_point;
Cone_tall3 = Cone_tall3 - Start_point;
Cone_orng3 = Cone_orng3 - Start_point;
End_line = End_line - Start_point;
End_point = End_point - Start_point;
Start_line = Start_line - Start_point;
Start_point = Start_point - Start_point;

track_rotation = 180 - track_rotation - atan2d(Start_line(2,1) - Start_line(1,1), Start_line(2,2) - Start_line(1,2));
Rotation_matrix = [cosd(track_rotation), -sind(track_rotation); sind(track_rotation), cosd(track_rotation)];

Data_flat(:,1:2) = (Rotation_matrix * Data_flat(:,1:2)')';
Cone_blue3 = fliplr((Rotation_matrix * flip(Cone_blue3'))');
Cone_yelw3 = fliplr((Rotation_matrix * flip(Cone_yelw3'))');
Cone_tall3 = fliplr((Rotation_matrix * flip(Cone_tall3'))');
Cone_orng3 = fliplr((Rotation_matrix * flip(Cone_orng3'))');
Start_line = fliplr((Rotation_matrix * flip(Start_line'))');
Start_point = fliplr((Rotation_matrix * flip(Start_point'))');
End_line = fliplr((Rotation_matrix * flip(End_line'))');
End_point = fliplr((Rotation_matrix * flip(End_point'))');

Data_flat(:,1:2) =  fliplr(track_offset) + Data_flat(:,1:2);
Cone_blue3 = track_offset + Cone_blue3;
Cone_yelw3 = track_offset + Cone_yelw3;
Cone_tall3 = track_offset + Cone_tall3;
Cone_orng3 = track_offset + Cone_orng3;
Start_line = track_offset + Start_line;
Start_point = track_offset + Start_point;
End_line = track_offset + End_line;
End_point = track_offset + End_point;

%% Finalise

if skidpad_test == false
    if acceleration_test == true
        End_line = Start_line + [acceleration_straight_length, 0];
        End_point = Start_point + [acceleration_straight_length, 0];
    else
        End_line = Start_line;
        End_point = Start_point;
    end
end
End_line(any(isnan(End_line), 2), :) = [];

if skidpad_test == true
    Start_line = (Start_line + End_line)/2;
    Start_point = (Start_point + End_point)/2;
    End_line = Start_line;
    End_point = Start_point;
    Cone_blue2 = Cone_blue3;
    Cone_blue3 = Cone_yelw3;
    Cone_yelw3 = Cone_blue2;
    Cone_tall3(1:2,1:2) = Start_line(1:2,:);
    Cone_tall3(3:4,1:2) = Start_line(1:2,:);
end

Cone_tall3(1:2,1:2) = Cone_tall3(1:2,1:2) + [tall_cone_sep/2, 0];
Cone_tall3(3:4,1:2) = Cone_tall3(3:4,1:2) - [tall_cone_sep/2, 0];
if acceleration_test == true
    Cone_tall3(5:6,1:2) = End_line(:,1:2) + [tall_cone_sep/2, 0];
    Cone_tall3(7:8,1:2) = End_line(:,1:2) - [tall_cone_sep/2, 0];
end
Cone_tall3(any(isnan(Cone_tall3), 2), :) = [];

figure(1)
axis equal
hold on
set(gca,'Color',[0.25, 0.25, 0.25])
% plot(Data_flat(1:end,2),Data_flat(1:end,1))
plot(Cone_blue3(1:end,1), Cone_blue3(1:end,2), 'b.')
plot(Cone_yelw3(1:end,1), Cone_yelw3(1:end,2), 'y.')
plot(Cone_tall3(1:end,1), Cone_tall3(1:end,2), 'o', 'Color', [1,0.5,0])
plot(Cone_orng3(1:end,1), Cone_orng3(1:end,2), '.', 'Color', [1,0.5,0])
plot(Start_line(1:end,1), Start_line(1:end,2), 'r')
plot(End_line(1:end,1), End_line(1:end,2), 'r')
% plot(Start_point(1:end,1), Start_point(1:end,2), 'k.')
% plot(End_point(1:end,1), End_point(1:end,2), 'k.')

%% Make files
%Note: This section creates the files that are used by Gazebo

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
for x = 1:size(Cone_tall3,1)
    fprintf(sdf, '    <link name=''link_%u''>\n', x+m+k);
    fprintf(sdf, '      <pose frame=''''>%f %f 0.070966 0 0 0</pose>\n', Cone_tall3(x,1), Cone_tall3(x,2));
    fprintf(sdf, '      <collision name="collision">\n');
    fprintf(sdf, '        <geometry>\n');
    fprintf(sdf, '          <mesh>\n');
    fprintf(sdf, '            <uri>model:%s</uri>\n', Cone_Model_tall);
    fprintf(sdf, '          </mesh>\n');
    fprintf(sdf, '        </geometry>\n');
    fprintf(sdf, '      </collision>\n');
    fprintf(sdf, '      <visual name="visual">\n');
    fprintf(sdf, '        <geometry>\n');
    fprintf(sdf, '          <mesh>\n');
    fprintf(sdf, '            <uri>model:%s</uri>\n', Cone_Model_tall);
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