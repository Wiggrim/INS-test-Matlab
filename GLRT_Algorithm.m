function [y, rank, kernel] = GLRT_Algorithm( H, C, variance, observe )
	
	% This function uses GLRT to test wether f equals zero in the model : x = Ha + Cf + w
	% Because of its linearity, result will be the same among Wald, Rao and GLRT
	
	H0 = H;
	H1 = [ H, C ];	% Here matrix H and matrix C will not cause linear appendence because the special design of H

	% Code below calculate the variable to be judged ---- 2ln(L(x))
	r0 = observe - H0*inv(H0'*inv(variance)*H0)*H0'*inv(variance)*observe;	% Using the weighted LS algorithm
	r1 = observe - H1*inv(H1'*inv(variance)*H1)*H1'*inv(variance)*observe;
	y = r0'*inv(variance)*r0 - r1'*inv(variance)*r1;
%     y = r0'*inv(variance)*r0;

	% Code below calculate the variable y's PDF ---- chi square
	rank = size( C, 2 );	% The unknown vector f's size
	inv_Fisher = inv(H1'*inv(variance)*H1);
	kernel = inv(inv_Fisher(end-rank+1:end,end-rank+1:end));	% The kernel should be of same property as Fisher matrix

end