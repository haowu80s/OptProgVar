function [ X, data_props, C ] = scale_center_phi( data_ZLV, data_props, C, z_width, type )
%scale_center_phi scales and centers the mass fraction and produce
%contrains and data matrices accordingly
%   input:
%       data_ZLV: 3D array (nz x nlam x nvar)
%       data_props: struct with info of data_ZLV
%       C: constraints matrix
%       z_width: number of mixture fraction values within each group
%       type: scaling type
%   output:
%       C: constraints matrix after scaling
%       X: cells of data matrix
%       data_props: struct with info of X
%
%% centralize data by mean
data_props.yMean = zeros(size(data_ZLV));
for iPhi=1:data_props.nPhi
    for iZ=1:data_props.nZ
        data_props.yMean(iZ,:,iPhi) = mean(data_ZLV(iZ,:,iPhi));
    end
end
%% scale data by type
data_props.yScale = zeros(data_props.nPhi,1);
for iPhi=1:data_props.nPhi
    if mean(reshape(data_ZLV(:,:,iPhi), 1, []))<sqrt(eps) || quantile(reshape(data_ZLV(:,:,iPhi), 1, []), 0.95)<sqrt(eps)
        data_props.yScale(iPhi) = 0;
    else
        if strcmp(type, 'std')
            data_props.yScale(iPhi) = std(reshape(data_ZLV(:,:,iPhi), 1, []));
        elseif strcmp(type, 'mean')
            data_props.yScale(iPhi) = mean(reshape(data_ZLV(:,:,iPhi), 1, []));
        elseif strcmp(type, 'median')
            data_props.yScale(iPhi) = median(reshape(data_ZLV(:,:,iPhi), 1, []));
        elseif strcmp(type, 'pareto')
            data_props.yScale(iPhi) = sqrt(std(reshape(data_ZLV(:,:,iPhi), 1, [])));
        elseif strcmp(type, '95')
            data_props.yScale(iPhi) = quantile(reshape(data_ZLV(:,:,iPhi), 1, []), 0.95);
        elseif strcmp(type, '99')
            data_props.yScale(iPhi) = quantile(reshape(data_ZLV(:,:,iPhi), 1, []), 0.99);
        elseif strcmp(type, 'max')
            data_props.yScale(iPhi) = max([std(reshape(data_ZLV(:,:,iPhi), 1, [])), ...
                mean(reshape(data_ZLV(:,:,iPhi), 1, [])), ...
                median(reshape(data_ZLV(:,:,iPhi), 1, [])), ...
                sqrt(std(reshape(data_ZLV(:,:,iPhi), 1, []))), ...
                median(reshape(data_ZLV(:,:,iPhi), 1, [])), ...
                quantile(reshape(data_ZLV(:,:,iPhi), 1, []), normcdf(0.5)) - quantile(reshape(data_ZLV(:,:,iPhi), 1, []), normcdf(-0.5))]);
        elseif strcmp(type, 'none')
            data_props.yScale(iPhi) = 1.0;
        end
    end
end
yScaleND = reshape(data_props.yScale, 1,1,data_props.nPhi);
yScaleND(yScaleND==0) = inf;
yScaleND = repmat(yScaleND, data_props.nZ, data_props.nLambda, 1);
Data_ZLam_tmp = (data_ZLV-data_props.yMean)./yScaleND;
[n, ~] = size(C);
yScaleND = reshape(data_props.yScale, 1, []);
yScaleND(yScaleND==0) = inf;
C = C ./ repmat(yScaleND, n, 1);
[N,~] = size(C);
P = randperm(N,N);
C = C(P,:);
%% create data cells X
Xtmp = permute(Data_ZLam_tmp,[2,3,1]);
nZGrid = ceil(data_props.nZ/z_width);
X{nZGrid} = [];
iZStar = 1;
for i=1:nZGrid
    iZEnd = min(iZStar+z_width-1,data_props.nZ);
    for iZ = iZStar:iZEnd
        X{i} = [X{i};squeeze(Xtmp(:,:,iZ))];
    end
    if (iZEnd-iZStar<=5)
        warning('iZEnd-iZStar<=5');
    end
    iZStar = iZStar+z_width;
end
end

