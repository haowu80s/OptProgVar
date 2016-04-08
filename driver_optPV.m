%%
% This is a sample driver PROGram for generating optimal PROG-var for FPV
% model from FlameMaster produced flamelets files
%%
% seed random number generator
rng(0, 'twister');
%% prepare flamelet solutions
% convert the flamelets to a 3D table of (Z, lam, var)
fl_dir = './gri30_nonpremixed/';
fl_prefix = 'CH4_p01';
[ data_ZLV,  data_props] = fl2zlv( fl_dir, fl_prefix );
% compute the matrix for monotonicity constraits
[ C ] = build_mono_cons( data_ZLV,  data_props );
%% precomputed data for easy testing
% load('./sample.mat');
%% scale and center variables
[ X, data_props, C ] = scale_center_phi( data_ZLV, data_props, C, 25, 'none' );
%% specify parameters for optimizeAB
opts = struct('printEvery',0,'tol',1e-6,'maxIts',5000,'restart',100, 'continuation', true);
tol_AB = 1e-5;
lam_l1 = 1000.0;
mu_l2 = 1;
%% run optimizeAB
[ A, B, R2, objFun, tfocs_out, tfocs_optsOut ] = optimizeAB( X, data_props, C, lam_l1, mu_l2, opts, tol_AB);
%% post processing B to get PROG
[ PROG ] = post_proc_B( B, data_props, tol_AB);
%% testing 
figure(1)
PROG_data = zeros(data_props.nZ, data_props.nLambda);
for i=1:1:data_props.nZ
    for j=1:data_props.nLambda
        PROG_data(i,j) = PROG'*reshape(data_ZLV(i,j,:),[],1);
    end
end
PROG_data = PROG_data/max(max(PROG_data));
surf(data_props.ZList, data_props.lambda_list, PROG_data'); shading flat; colorbar;
view([0,0,90])
xlabel('Z')
ylabel('\Lambda')
title('PROG')
spec_name = 'Y_CO2';
figure(2);
hold on ;
for Z=[0.01, 0.055, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5]
    [~,iZ] = min((data_props.ZList-Z).^2);
    Ctest = PROG_data(iZ,:);
    Ytest = data_ZLV(iZ,:,data_props.mapIndex(spec_name));
    plot(Ctest, Ytest, '.-')
end
hold off;
xlabel('PROG');
ylabel(strrep(spec_name, 'Y_', ''));