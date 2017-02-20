function [ judge, judge_another, gx, rs, direction ] = RS_RAIM_check( RS_test_structure )
    
    % This function is used to do RS-RAIM test of data extracted from INS
    % positioning equation and kalman-filter
    
    transform_matrix_between_state = RS_test_structure.trans_matrix;
    index = RS_test_structure.index;
    observe_state_vector = RS_test_structure.ob_vector;
    state_variance = RS_test_structure.ob_variance;
    fai = RS_test_structure.fai;
    latch = RS_test_structure.latch;
    
    index(1) = 1;
    
    % extract useful data, which carrisponding to the KF renewing step
    transform_matrix_between_state = transform_matrix_between_state(:,:,index~=0);
    observed_state_vector = observe_state_vector(:,index~=0);
    state_variance = state_variance(:,:,index~=0);
    fai = fai(:,:,index~=0);
    latch = latch(:,:,index~=0);
    
    % this variable control how long the RS-test function is (ex. 3 KF states)
    used_state_num = 3;
    
    judge = zeros(size(observed_state_vector,2)-used_state_num,1);
    
    judge_another = zeros(size(observed_state_vector,2)-used_state_num,1);
    
    gx = zeros(size(observed_state_vector,2)-used_state_num,1);
    rs = zeros(used_state_num*15-9,size(observed_state_vector,2)-used_state_num);
    direction = zeros( used_state_num*15-9,size(observed_state_vector,2)-used_state_num );
    
    % do RS-test to all KF states
    for cnt = 1 : size(observed_state_vector,2)-used_state_num
    
        % formulate the observe matrix
        H = zeros( used_state_num*15-9, used_state_num*9 );
        for k = 1 : used_state_num-1
            H( (k-1)*15+1:(k-1)*15+6,(k-1)*9+1:(k-1)*9+6 ) = eye(6,6);
            H( (k-1)*15+7:(k-1)*15+15,(k-1)*9+1:(k-1)*9+18 ) = [transform_matrix_between_state(:,:,k+cnt),-eye(9,9)];
        end
        H((used_state_num-1)*15+1:(used_state_num-1)*15+6,(used_state_num-1)*9+1:(used_state_num-1)*9+6) = eye(6,6);

        % formulate the observe vector's variance matrix
        V = zeros( used_state_num*15-9, used_state_num*15-9 );
        for k = 1 : used_state_num-1
            V( (k-1)*15+1:(k-1)*15+6, (k-1)*15+1:(k-1)*15+6 ) = eye(6,6)*0.5 + state_variance(1:6,1:6,k+cnt-1);
            V( (k-1)*15+1:(k-1)*15+6, (k-1)*15+7:(k-1)*15+15 ) = (fai(:,:,k+cnt)*latch(:,:,k+cnt)*eye(6,6)*0.5)';
            V( (k-1)*15+7:(k-1)*15+15, (k-1)*15+7:(k-1)*15+15 ) = state_variance(:,:,k+cnt);
            V( (k-1)*15+7:(k-1)*15+15, (k-1)*15+1:(k-1)*15+6 ) = fai(:,:,k+cnt)*latch(:,:,k+cnt)*eye(6,6)*0.5;
        end
        V( (used_state_num-1)*15+1:(used_state_num-1)*15+6, (used_state_num-1)*15+1:(used_state_num-1)*15+6 ) = eye(6,6)*0.5+ state_variance(1:6,1:6,used_state_num+cnt-1);

        % formulate the observed vector consist of delta_x, delta_v
        ob_v = zeros( used_state_num*15-9, 1 );
        A = zeros( used_state_num*15-9, used_state_num*6 );
        idx = zeros( used_state_num*15-9, 1 );
        for k = 1 : used_state_num
            ob_v( (k-1)*15+1:(k-1)*15+6 ) = observed_state_vector(:,k+cnt-1);
            idx( (k-1)*15+1:(k-1)*15+6 ) = 1;
            A( (k-1)*15+1:(k-1)*15+6, (k-1)*6+1:(k-1)*6+6 ) = eye(6,6);
        end
        
        C = eye( used_state_num*15-9,used_state_num*15-9 );
        C = C(:,boolean(idx));
        
        [judge_another(cnt),~] = GLRT_Algorithm( H, C, V, ob_v );
        
        rs(:,cnt) = ob_v - H / (H'/ V * H) * H'/(V) * ob_v;
        judge(cnt) = rs(:,cnt)' / V * rs(:,cnt);

        S = inv( H' * inv(V) * H ) * H' * inv(V);
        s_t = S((used_state_num-1)*6+4,:);
        Mx = s_t * A;
        Ma = sqrt( A' * inv(V) * (eye(used_state_num*15-9,used_state_num*15-9)-H*S) * A );
        gx(cnt) = norm(Mx*Ma)^2;    % this gx gives the WMS algorithm's performance, but due to my KF uses the LMMSE algorithm, so this gx value are meanless
        
        direction(:,cnt) = A * Ma * Ma' * Mx';
    end
    
end

