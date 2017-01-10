clear;
clc;
close all;

start_position = [ 39.9, 116.3 ] / 180 * pi;
end_position = [ 41.3, 123.73 ] / 180 * pi;
velocity = 500;
height = 10000;
earth_flag = 1;
delta_l = 0.1;

[ accel_measurement, gryo_measurement, delta_t, start_p, start_v, start_ati, lat_prof, lon_prof, height_prof, v_prof, yaw_prof ] = genIMUMeasurement( start_position, end_position, velocity, height, earth_flag, delta_l );

gps_data = [lat_prof',lon_prof',height_prof', v_prof];

gps_data(:,3:end) = gps_data(:,3:end) + randn(size(gps_data(:,3:end)))*sqrt(0.005);

% index = randi(size(gps_data,1));
% tmp = 1:size(gps_data,1);
% gps_data(tmp>index,1:2) = gps_data(tmp>index,1:2) + 0.5;

accel_noise = sin( (1:length(accel_measurement))/length(accel_measurement)*4*pi )'*randn(1,3)*0.05;
accel_measurement = accel_measurement + accel_noise;

[ ins_position, ins_velocity, ins_attitude, rs ] = INSUpdate( gryo_measurement, accel_measurement, delta_t, gps_data, start_p, start_v, start_ati );