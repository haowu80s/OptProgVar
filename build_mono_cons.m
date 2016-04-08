function [ C ] = build_mono_cons( data_ZLV,  data_props )
%build_mono_cons creats the matrix for monotonica constrains
%   input:
%       data_ZLV: 3D array (nz x nlam x nvar)
%       data_props: struct with info of data_ZLV
%   output:
%       C: C is a matrix of size n x p
%           n is the total number of points in the flamelet solutions, i.e.
%           nz*nlam. p is the total number of species, c_ij =
%           d(phi_j) / d(lambda) @ point i
%%
data_dlambda = zeros(size(data_ZLV));
lambda_list = data_props.lambda_list;
%%
h = waitbar(0,'Calculating the gradient');
steps = data_props.nPhi*data_props.nZ;
step = 0;
for iPhi=1:data_props.nPhi
    for iZ=1:data_props.nZ
        % waitbar
        step = step+1;
        waitbar(step / steps)
        % calculation
        y = squeeze(data_ZLV(iZ, :, iPhi));
        if mean(y)>sqrt(eps)
            sp = spaps(lambda_list, y/mean(y), 1e-5);
            sp2 = fnder(sp,1);
            y_d = fnval(sp2, lambda_list)*mean(y);
            data_dlambda(iZ, :, iPhi) = smooth(lambda_list, y_d,'lowess');
        else
            data_dlambda(iZ, :, iPhi) = 0;
        end
    end
end
close(h)
%%
C = zeros(data_props.nZ*data_props.nLambda, data_props.nPhi);
for iPhi=1:data_props.nPhi
    C(:,iPhi) = reshape(data_dlambda(:,:,iPhi),[],1);
end
