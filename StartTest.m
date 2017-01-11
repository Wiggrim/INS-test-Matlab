clear;
clc;
close all;

start_position = [ 39.9, 116.3 ] / 180 * pi;
end_position = [ 41.3, 123.73 ] / 180 * pi;
velocity = 500;
height = 10000;
earth_flag = 1;
delta_l = 0.1;

% get IMU measurement and other profile
[ IMU_measurement, delta_t, start_info, gps_data ] = genIMUMeasurement( start_position, end_position, velocity, height, earth_flag, delta_l );

% extract data
accel_measurement = IMU_measurement(:,1:3);
gryo_measurement = IMU_measurement(:,4:6);
start_p = start_info(1,:);
start_v = start_info(2,:);
start_ati = start_info(3,:);

%
%   code below used to add different kinds of noise to gps_data and
%   IMU_measurement manualy
gps_data(:,3:end) = gps_data(:,3:end) + randn(size(gps_data(:,3:end)))*sqrt(0.005);

% index = randi(size(gps_data,1));
% tmp = 1:size(gps_data,1);
% gps_data(tmp>index,1:2) = gps_data(tmp>index,1:2) + 0.5;

accel_noise = sin( (1:length(accel_measurement))/length(accel_measurement)*4*pi )'*randn(1,3)*0.05;
accel_measurement = accel_measurement + accel_noise;
%

[ ins_position, ins_velocity, ins_attitude, recorder, RS_test ] = INSUpdate( gryo_measurement, accel_measurement, delta_t, gps_data, start_p, start_v, start_ati );

[judge, gx, rs] = RS_RAIM_check( RS_test );