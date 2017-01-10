function [ accel_measurement, gryo_measurement, delta_t, start_p, start_v, start_ati, lat_prof, lon_prof, height_prof, velocity_prof, yaw_prof ] = genIMUMeasurement( start_position, end_position, velocity, height, wcg_flag, delta_l )

	% This matlab file is used to build the measurement of gryo and accel
	% The gryo and accel's measurement are both relative to the I-frame
	% Which is to say the measurement of gryo = w_b2l + w_l2e + w_e2i
	% And the measurement of accel = d(ve)/dt + cori_bias
	% One thing needed to be noticed is changing the reference frame is different from changing the axis

	HEIGHT = height;
	EARTH_FLAG = wcg_flag; % use the WCG
	VELOCITY = velocity; % nm/h
	START_POSITION = start_position;  %[ 39+54/60+27/3600, 116+23/60+17/3600 ] * pi / 180;
	END_POSITION = end_position; %[ 41.3, 123.73 ] * pi / 180;

	% Generate the footprint
	[ lat_prof, lon_prof, tc_prof, total_dist ] = greatcir( START_POSITION, END_POSITION, HEIGHT, EARTH_FLAG, delta_l );

	% Something needed in calculating the accel and gryo measurement
	prof_length = length( lat_prof );
	total_time = total_dist / VELOCITY * 3600;
	time_scale = ( 0:prof_length-1 )*( total_time/( prof_length-1 ) );
	velocity_norm_prof = VELOCITY*ones( 1, prof_length );	% nm/h
	height_prof = HEIGHT*ones( 1, prof_length );

	% Calculate the gryo measurement
	DCMnb_prof = dcmnbgen(tc_prof);	% The angle is expressed by transform matrix, from NED axis to body axis
	dthetbody = gendthet(DCMnb_prof);	% The w_b2l part
	DCMel_prof = dcmelgen(lat_prof, lon_prof, tc_prof);
	deltaer = earthrot(time_scale,DCMel_prof,DCMnb_prof);	% The w_e2i part
	deltacr = gendelcr(lat_prof,tc_prof,velocity_norm_prof, height_prof,time_scale,DCMnb_prof,DCMel_prof,EARTH_FLAG);	% The w_l2e part
	gryo_measurement = dthetbody; %dthetbody + deltacr + deltaer;
	gryo_error = zeros( size(gryo_measurement) );
	gryo_measurement = gryo_measurement + gryo_error;

	% Calculate the accel measurement
	velocity_prof = genvelpr(tc_prof,velocity_norm_prof);	% This velocity is in ENU axis
	deltav_b = gendv(velocity_prof,DCMnb_prof);	% d(ve)/dt
	dvcor = gendvcor(lat_prof,velocity_norm_prof,tc_prof,height_prof,time_scale,DCMnb_prof,DCMel_prof,EARTH_FLAG);	% cori_bias as well as gravity
	accel_measurement = deltav_b + dvcor;
	accel_error = zeros( size(accel_measurement) );
	accel_measurement = accel_measurement + accel_error;

	delta_t = (time_scale(2) - time_scale(1));
	v_ps = VELOCITY * 1.6878 * 0.3048;
	start_p = [ lat_prof(1), lon_prof(1), HEIGHT ];
	start_v = [ v_ps*cos(tc_prof(1)), v_ps*sin(tc_prof(1)), 0 ];	% NED
	start_ati = [ 0, 0, tc_prof(1) ];
    
    accel_measurement = accel_measurement / delta_t;
    gryo_measurement = gryo_measurement / delta_t;
    
    velocity_prof = [velocity_prof(:,2),velocity_prof(:,1),-velocity_prof(:,3)];
    
    yaw_prof = tc_prof;


end
















