function [ lambda_list,  fl_order] = props2lam( props )
%props2lam converts the properties of flamelets (Tmax and chi_st) to
%lambda. lambda is defined as the arc-length along the s-curve. Thus, the
%combition of mixture fraction and lambda forms the coordinate chart for
%the FPV manifold, i.e phi = phi (Z, lambda) maps all the flamelet
%solutions
%   input:
%       props: output from readFMFiles
%   output:
%       lambda_list: array of calculated lambdas
%       fl_order: order of flamelets according to lambda

% intialize variables
nFlames = length(props);
lambda_list = zeros(1,nFlames);
chi_st_list = zeros(1,nFlames);
t_max_list = zeros(1,nFlames);
index_list = (1:nFlames);
% obtain list of Tmax and chi_st
for i=1:nFlames
    chi_st_list(i) = log10(props(i).chi_st);
    t_max_list(i) = props(i).Tmax;
end
% find the starting point for the s-curve
% starting point is defined as low chi_st and high t_max
t_max_ref = mean(t_max_list);
chi_st_ref = 1;
StartPointList = t_max_list/t_max_ref - chi_st_list/chi_st_ref;
[~, iStarFlame] = max(StartPointList);
% walk along the s-curve and build lambda
i_progress = 1;
ind = iStarFlame;
iFlame = index_list(ind);
lambda_list(iFlame) = max(StartPointList);
t_max_old = t_max_list(iFlame);
chi_st_old = chi_st_list(iFlame);
lambda_old = lambda_list(iFlame);
t_max_list(iFlame) = [];
chi_st_list(iFlame) = [];
index_list(iFlame) = [];
while true
    if isempty(t_max_list)
        break;
    end
    i_progress = i_progress+1;
    % find the next point that is closest to the previous one
    [d_lambda,ind] = min(sqrt((t_max_list/t_max_ref-t_max_old/t_max_ref).^2+(chi_st_list/chi_st_ref-chi_st_old/chi_st_ref).^2));
    iFlame = index_list(ind);
    lambda_list(iFlame) = lambda_old-d_lambda;
    t_max_old = t_max_list(ind);
    chi_st_old = chi_st_list(ind);
    lambda_old = lambda_list(iFlame);
    t_max_list(ind) = [];
    chi_st_list(ind) = [];
    index_list(ind) = [];
end
% clear temp variabels
clear chiStList TmaxList indexList Tmax_old chiSt_old Lambda_old;
% sourt according to lambda
[lambda_list, fl_order] = sort(lambda_list,'ascend');
end

