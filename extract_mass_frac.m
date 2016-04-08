function [data_ZLV, data_props ] = extract_mass_frac( data, index_map, lambda_list, fl_order )
%extract_mass_frac removes variables that are not species mass fractions
%from the array of flamelet solutions and rebuilt data_props
%solutions
%   input:
%       data: cell of flamelet solutions
%       index_map: a map from variable name to variable index
%       lambda_list: array of calculated lambdas
%       fl_order: order of flamelets according to lambda
%   output:
%       data_ZLV: 3D array (nz x nlam x nvar)
%       data_props: struct with info of data_ZLV
%%
data_props.ZList = data{fl_order(1)}{index_map('Z')};
data_props.lambda_list = lambda_list;

data_props.nZ = length(data_props.ZList);
data_props.nLambda = length(lambda_list);
%% Filter out Phi's
data_props.nPhi = 0;
name_tmp = keys(index_map);
phi_list = {};
for i=1:length(index_map)
    if ~isempty(strfind(char(name_tmp(i)), 'Y_')) || strcmp(char(name_tmp(i)), 'T')
        data_props.nPhi=data_props.nPhi+1;
        phi_list = [phi_list,name_tmp(i)];
    end
end
%%
data_ZLV = zeros(data_props.nZ, data_props.nLambda, data_props.nPhi);
for i=1:data_props.nLambda
    for j=1:data_props.nPhi
        ytemp = data{fl_order(i)}{index_map(phi_list{j})};
        if mean(ytemp)>sqrt(eps)
            ztemp = data{fl_order(i)}{index_map('Z')};
            y_temp2 = interp1(ztemp, ytemp, data_props.ZList, 'spline');
            data_ZLV(:,i,j) = y_temp2;
        else
            data_ZLV(:,i,j) = 0;
        end
    end
end
%%
data_props.mapIndex = containers.Map(phi_list,1:data_props.nPhi);
%%
data_props.phiNames = phi_list';
end

