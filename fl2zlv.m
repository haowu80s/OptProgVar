function [ data_ZLV,  data_props] = fl2zlv( fl_dir, fl_prefix )
%fl2zlv convert the flamelets in a given directories to a 3D table of (Z, lam, var)
%   currenty supports FlameMaster output only, can replace readFMFiles with
%   other user defined function for reading files of other format
%   input:
%       fl_dir: dirctory of flamelet files
%       fl_prefix: prefix shared by flamelet files (to distinguish from
%       other files in fl_dir 
%   output:
%       data_ZLV: 3D array (nz x nlam x nvar)
%       data_props: struct with info of data_ZLV

%% read flamelets
[ data, props, index_map ] = readFMFiles( fl_dir, fl_prefix );
%% compute lambda from props
[ lambda_list,  fl_order] = props2lam( props );
%% filter out mass fractions
[data_ZLV, data_props ] = extract_mass_frac( data, index_map, lambda_list, fl_order );
