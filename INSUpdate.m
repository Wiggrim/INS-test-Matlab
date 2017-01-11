function [ ins_position, ins_velocity, ins_attitude, recorder, RS_test ] = INSUpdate( gryo_measurement, accel_measurement, delta_t, gps_data, start_position, start_velocity, start_attitude )
	
	% This file is used to do INS update using the measurement got from simulation
	% The gryo_measurement is in body axis, and is calculated by delta_theta/dt
	% The accel_measurement is in body axis, and is calculated by delta_v/dt
	% The start_position is in lat, lon, height
	% The start_velocity is in NED axis, although the velocty in INS toolbox is in ENU axis
	% The start_attitude describe the relationship between the NED axis and body axis
	
	EARTH_A = 6378137.0;
	EARTH_E = 0.0818;
	EARTH_W = 7.2921155E-5;
	
	total_length = length( gryo_measurement ) + 1;
    
    recorder = zeros(3,total_length);
    
	time_scale = (0:total_length-1)*delta_t;
	ins_position = zeros( total_length, 3 );
	ins_velocity = zeros( total_length, 3 );
	ins_attitude = zeros( total_length, 3 );
	
	ins_position(1,:) = start_position;
	ins_velocity(1,:) = start_velocity;
	ins_attitude(1,:) = start_attitude;
	
	R_b2n = zeros(3,3); % transform body axis to NED axis
	w_n2i_pseu = zeros(3,1);    % rotation between n-frame and i-frame
	R = zeros(3,3); % transform rotation rate to eular angle rate
    
    % data space alloc for kalman-filter
    accel_accum_bias = zeros(3,1);
    F = zeros(9,9); % disturbation matrix of INS positioning equation
	noise_w = zeros(9,9);   % kalman-filter's state transformation noise variance
    noise_w(7,7) = 0.05;    % The accel_measurement's residue bias has in-coming gauss noise
    noise_w(8,8) = 0.05;
    noise_w(9,9) = 0.05;
    ins_variance = zeros(9,9,total_length); % the variance matrix of delta_p, delta_v, delta_a
    ins_variance(:,:,1) = eye(9,9)*0.005;
    ins_variance(1,1,1) = 1E-4; % position error is 100 meters
    ins_variance(2,2,1) = 1E-4;
    ins_variance_latch = zeros(9,9,total_length);
    observe_variance = eye(6,6)*0.005;  % the observe data variance matrix
    observe_variance(1,1) = 1E-4;
    observe_variance(2,2) = 1E-4;
    
    % data space alloc for RS-RAIM
    fai = eye(9,9);
    rs = zeros( 21, total_length );
    H = zeros( 9, 9, total_length );
    observe_vector = zeros( 6, total_length );
    gx = zeros(total_length,1);
    latch = zeros(9,6);
    judge = zeros(total_length,1);
    index = zeros(1,total_length);
	
	for k = 2 : total_length
	
		RM = EARTH_A*( 1-EARTH_E*EARTH_E ) / (1-(EARTH_E*sin(ins_position(k-1,1)))^2)^(3/2);
		RN = EARTH_A * ( 1-(EARTH_E*sin(ins_position(k-1,1)))^2 )^(-1/2);
		D = [ 1/(RM+ins_position(k-1,1)), 0, 0; 
			  0, 1/(RN+ins_position(k-1,1))/cos(ins_position(k-1,1)), 0;
			  0, 0, -1 ];
			  
		R_b2n(1,1) = cos(ins_attitude(k-1,3))*cos(ins_attitude(k-1,2));
		R_b2n(1,2) = cos(ins_attitude(k-1,3))*sin(ins_attitude(k-1,2))*sin(ins_attitude(k-1,1))-sin(ins_attitude(k-1,3))*cos(ins_attitude(k-1,1));
		R_b2n(1,3) = cos(ins_attitude(k-1,3))*sin(ins_attitude(k-1,2))*cos(ins_attitude(k-1,1)) + sin(ins_attitude(k-1,3))*sin(ins_attitude(k-1,1));
		R_b2n(2,1) = sin(ins_attitude(k-1,3))*cos(ins_attitude(k-1,2));
		R_b2n(2,2) = sin(ins_attitude(k-1,3))*sin(ins_attitude(k-1,2))*sin(ins_attitude(k-1,1)) + cos(ins_attitude(k-1,3))*cos(ins_attitude(k-1,1));
		R_b2n(2,3) = sin(ins_attitude(k-1,3))*sin(ins_attitude(k-1,2))*cos(ins_attitude(k-1,1)) - cos(ins_attitude(k-1,3))*sin(ins_attitude(k-1,1));
		R_b2n(3,1) = -sin(ins_attitude(k-1,2));
		R_b2n(3,2) = cos(ins_attitude(k-1,2))*sin(ins_attitude(k-1,1));
		R_b2n(3,3) = cos(ins_attitude(k-1,2))*cos(ins_attitude(k-1,1));
		
		w_n2i_pseu(1) = 2*EARTH_W*cos(ins_position(k-1,1)) + ins_velocity(k-1,2)/(RN+ins_position(k-1,3));
		w_n2i_pseu(2) = -ins_velocity(k-1,1)/(RM+ins_position(k-1,3));
		w_n2i_pseu(3) = -2*EARTH_W*sin(ins_position(k-1,1)) - ins_velocity(k-1,2)*tan(ins_position(k-1,1))/(RN+ins_position(k-1,3));
		
		R(1,1) = 1;
		R(1,2) = sin(ins_attitude(k-1,1))*tan(ins_attitude(k-1,2));
		R(1,3) = cos(ins_attitude(k-1,1))*tan(ins_attitude(k-1,2));
		R(2,2) = cos(ins_attitude(k-1,1));
		R(2,3) = -sin(ins_attitude(k-1,1));
		R(3,2) = cos(ins_attitude(k-1,2))*sin(ins_attitude(k-1,1));
		R(3,3) = cos(ins_attitude(k-1,2))*cos(ins_attitude(k-1,1));
        
        % formulate the disturb matrix F, assume the ins_position is close
        % to the reality
        F = F.*0;
        F(1,3) = -ins_velocity(k-1,1)/(RM+ins_position(k-1,3))^2;
        F(2,1) = ins_velocity(k-1,2)*sin(ins_position(k-1,1))/(RN+ins_position(k-1,3))/cos(ins_position(k-1,1))/cos(ins_position(k-1,1));
        F(2,3) = -ins_velocity(k-1,2)/((RN+ins_position(k-1,3))^2)/cos(ins_position(k-1,1));
        F(1:3,4:6) = D;
        F(4,1) = -2*ins_velocity(k-1,2)*EARTH_W*cos(ins_position(k-1,1));
        F(4,3) = -ins_velocity(k-1,1)*ins_velocity(k-1,3)/((RM+ins_position(k-1,3))^2) + ins_velocity(k-1,2)*ins_velocity(k-1,2)*tan(ins_position(k-1,1))/((RN+ins_position(k-1,3))^2);
        F(5,1) = 2*EARTH_W*(ins_velocity(k-1,1)*cos(ins_position(k-1,1))-ins_velocity(k-1,3)*sin(ins_position(k-1,1))) + ins_velocity(k-1,2)*ins_velocity(k-1,1)/(RN+ins_position(k-1,3))/cos(ins_position(k-1,1))/cos(ins_position(k-1,1));
        F(5,3) = -ins_velocity(k-1,2)*(ins_velocity(k-1,3)+ins_position(k-1,1)*tan(ins_position(k-1,1)))/((RN+ins_position(k-1,3))^2);
        F(6,1) = 2*ins_velocity(k-1,2)*EARTH_W*sin(ins_position(k-1,1));
        F(6,3) = ((ins_velocity(k-1,2))^2)/((RN+ins_position(k-1,3))^2) + ((ins_velocity(k-1,1))^2)/((RM+ins_position(k-1,3))^2);
        F(4,4) = ins_velocity(k-1,3)/(RM+ins_position(k-1,3));
        F(4,5) = -2*EARTH_W*sin(ins_position(k-1,1)) - 2*ins_velocity(k-1,2)*tan(ins_position(k-1,1))/(RN+ins_position(k-1,3));
        F(4,6) = ins_velocity(k-1,1)/(RM+ins_position(k-1,3));
        F(5,4) = 2*EARTH_W*sin(ins_position(k-1,1)) + ins_velocity(k-1,2)*tan(ins_position(k-1,1))/(RN+ins_position(k-1,3));
        F(5,5) = (ins_velocity(k-1,3)+ins_velocity(k-1,1)*tan(ins_position(k-1,1)))/(RN+ins_position(k-1,3));
        F(5,6) = 2*EARTH_W*cos(ins_position(k-1,1)) + ins_velocity(k-1,2)/(RN+ins_position(k-1,3));
        F(6,4) = -2*ins_velocity(k-1,1)/(RM+ins_position(k-1,3));
        F(6,5) = -2*EARTH_W*cos(ins_position(k-1,1)) - 2*ins_velocity(k-1,2)/(RN+ins_position(k-1,3));
        F(4:6,7:9) = R_b2n;
        F(7,7) = 0.005;
        F(8,8) = 0.005;
        F(9,9) = 0.005;
        
		% the ins positioning function
		ins_position(k,:) = ( ins_position(k-1,:)' + delta_t * ( D * ins_velocity(k-1,:)' ) )';
		ins_velocity(k,:) = ins_velocity(k-1,:) + delta_t*( R_b2n*((accel_measurement(k-1,:)-accel_accum_bias')') - cross(w_n2i_pseu, ins_velocity(k-1,:)') + [0;0;gravity(ins_position(k-1,1),ins_position(k-1,3))])';
		ins_attitude(k,:) = ins_attitude(k-1,:) + delta_t*( R*(gryo_measurement(k-1,:))' )';
        
        % the variance update funtion
        ins_variance(:,:,k) = (eye(9,9)+delta_t*F) * ins_variance(:,:,k-1) * (eye(9,9)+delta_t*F') + noise_w;
        ins_variance_latch(:,:,k) = ins_variance(:,:,k);
        fai = fai * (eye(9,9)+delta_t*F);
        
        % do feedback every 7 seconds
        if mod(k-1,10) == 0
        
            % the kalman-filter feedback
            observe_vector(:,k) = ([ins_position(k,:),ins_velocity(k,:)]-gps_data(k,1:6))';
            temp = ins_variance(:,:,k)*eye(9,6)*inv(eye(6,9)*ins_variance(:,:,k)*eye(9,6)+observe_variance);
            ins_variance(:,:,k) = ins_variance(:,:,k) - ins_variance(:,:,k)*eye(9,6)*inv(eye(6,9)*ins_variance(:,:,k)*eye(9,6)+observe_variance)*eye(6,9)*ins_variance(:,:,k);
            ins_position(k,:) = ins_position(k,:) - (temp(1:3,:)*observe_vector(:,k))';
            ins_velocity(k,:) = ins_velocity(k,:) - (temp(4:6,:)*observe_vector(:,k))';
            accel_accum_bias = accel_accum_bias + (temp(7:9,:)*observe_vector(:,k));

            % formulate the observe matrix H, its size is (21,18)
%             H(:,:,k) = zeros( 21, 18 );
%             H(1:6,1:6) = eye(6,6);
%             H(7:15,1:9, k) = fai*(eye(9,9)-latch*eye(6,9));
%             H(7:15,10:18, k) = -eye(9,9);
%             H(16:21,10:15, k) = eye(6,6);
            H(:,:,k) = fai * ( eye(9,9) - latch*eye(6,9) );

            %formulate the variance matrix V, its size is (21,21)
%             V = zeros(21,21);
%             V(1:6,1:6) = eye(6,6)*0.005;
%             V(7:15,7:15) = ins_variance(:,:,k);
%             V(16:21,16:21) = eye(6,6)*0.005;

            % formulate the observe vector, which consists of delta_x1,
            % delta_v1, zero*9, delta_x2, delta_v2
%             observe_vector_rs = zeros(21,1);
%             observe_vector_rs(1:6) = observe_vector(:,k-10);
%             observe_vector_rs(16:21) = observe_vector(:,k);

            % calculate the residue vector, assume the noise is of standard
            % gauss distribution
%             rs(:,k) = observe_vector_rs - H(:,:,k) * inv(H(:,:,k)'/V*H(:,:,k)) * H(:,:,k)'/(V)*observe_vector_rs;
%             judge(k) = rs(:,k)' / V * rs(:,k);

            %calculate the gx
%             S = inv( H(:,:,k)' * inv(V) * H(:,:,k) ) * H(:,:,k)' * inv(V);
%             A = zeros(21,12);
%             A(1:6,1:6) = eye(6,6);
%             A(16:21,7:12) = eye(6,6);
%             s_t = S(10,:);
%             Mx = s_t * A;
%             Ma = sqrt( A' * inv(V) * (eye(21,21)-H(:,:,k)*S) * A );
%             gx(k) = norm(Mx*Ma)^2;

            latch = temp;
            fai = eye(9,9);
            
            index(k) = 1;
        end
		
        recorder(:,k) = accel_accum_bias;
		
    end
	
    RS_test.trans_matrix = H;
    RS_test.ob_vector = observe_vector;
    RS_test.ob_variance = ins_variance_latch;
    RS_test.index = index;

end