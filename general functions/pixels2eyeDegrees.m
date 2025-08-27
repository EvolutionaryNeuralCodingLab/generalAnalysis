%% Pixels to degreess
%%Dell screen specs

%%Acer active display area 587x330 mm
%%pixel_size = 33/(1080)

function [theta_x theta_y] = pixels2eyeDegrees(eye_to_monitor_distance,pixel_size,monitor_resolution)
% Parameters
% eye_to_monitor_distance = 21.5; % Distance from eye to monitor in cm
% pixel_size = 33/1080; % Size of one pixel in cm (e.g., 25 micrometers)
% monitor_resolution = [1920, 1080]; % Width and height in pixels
screen_center = monitor_resolution / 2; % Center pixel coordinates

% Create a grid of pixel coordinates
[x, y] = meshgrid(1:monitor_resolution(1), 1:monitor_resolution(2));

% Calculate the pixel distance from the center
delta_x = (x - screen_center(1)) * pixel_size;
delta_y = (y - screen_center(2)) * pixel_size;

% Convert pixel distances to degrees of visual angle
theta_x = 2 * atan(delta_x / (2 * eye_to_monitor_distance)) * (180 / pi);
theta_y = 2 * atan(delta_y / (2 * eye_to_monitor_distance)) * (180 / pi);

end

