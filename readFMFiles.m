function [ data, props, index_map ] = readFMFiles( fl_dir, fl_prefix )
%readFMFiles reads FlameMaster flamelet files in fl_dir with prefix
%   fl_prefix. Support of orther fileformats can be achieved by modifying
%   readFM.
%   input:
%       fl_dir: dirctory of flamelet files
%       fl_prefix: prefix shared by flamelet files (to distinguish from
%       other files in fl_dir
%   output:
%       data: cell of flamelet solutions
%       props: cell of flamelet properties
%       index_map: a map from variable name to variable index

files = dir([fl_dir, fl_prefix, '*']);
nFlames = length(files);
props = [];
data = cell(1,nFlames);
names = cell(1,nFlames);
for i = 1:nFlames
    files(i).name
    [tmp, tmp1,  tmp2 ] = readFM([fl_dir,files(i).name]);
    props = [ props tmp2];
    data{i} = tmp;
    names{i} = tmp1;
end

clear tmp tmp1 tmp2
% build dictionary for species
index_map = dictName(names{1});
end

%% helper function: read flamelet file from filename
function [ data, names, opt] = readFM(filename)
fid = fopen(filename,'r');
tline = fgetl(fid);
if ( tline == 'header')
else
    error('Malformed Flamelets files');
end
notBody = 1;

opt = struct;
while ( notBody )
    tline = fgetl(fid);
    tmp = regexp(tline,'=','split');
    %tmp{1}
    switch strtrim(tmp{1})
        case {'body'}
            notBody = 0;
        case {'title'}
            opt.('title')=tmp{2};
        case {'mechanism'}
            opt.('mech') = tmp{2};
        case {'date'}
            opt.('date') = tmp{2};
        case {'fuel'}
            opt = setfield(opt,'fuel',tmp{2});
        case {'pressure'}
            val = sscanf(tmp{2},'%e');
            opt = setfield(opt,'P',val*1e5);
        case {'Z_st'}
            val = sscanf(tmp{2},'%e');
            opt = setfield(opt,'Zst',val);
        case {'chi_st'}
            val = sscanf(tmp{2},'%e');
            opt = setfield(opt,'chi_st',val);
        case {'FlameLoc'}
            val = sscanf(tmp{2},'%e');
            opt = setfield(opt,'floc',val);
        case {'Tmax'}
            val = sscanf(tmp{2},'%e');
            opt = setfield(opt,'Tmax',val);
        case {'numOfSpecies'}
            val = sscanf(tmp{2},'%d');
            opt = setfield(opt,'nsps',val);
        case {'gridPoints'}
            n = sscanf(tmp{2},'%d');
            opt = setfield(opt,'nx',n);
        case {'FuelSide'}
            tline = fgetl(fid);
            if ( ~strcmp(tline,'begin') )
                error('malformed FuelSide');
            end
            tline = fgetl(fid);
            while ( ~strcmp(tline,'end'))
                tmp1 = regexp(tline,'=','split');
                [token,remain] = strtok(strtrim(tmp1{1}),'-');
                remain = strrep(remain, '-', '_');
                switch token
                    case {'Temperature'}
                        val = sscanf(tmp1{2},'%e');
                        opt.('Fside').('Temp') =val;
                    case {'Massfraction'}
                        val = sscanf(tmp1{2},'%e');
                        opt.('Fside').(remain(2:end)) = val;
                end
                tline = fgetl(fid);
            end
        case {'OxidizerSide'}
            tline = fgetl(fid);
            if ( ~strcmp(tline,'begin') )
                error('malformed OxidizerSide');
            end
            tline = fgetl(fid);
            while ( ~strcmp(tline,'end'))
                tmp1 = regexp(tline,'=','split');
                [token,remain] = strtok(strtrim(tmp1{1}),'-');
                switch token
                    case {'Temperature'}
                        val = sscanf(tmp1{2},'%e');
                        opt.('Oside').('Temp') =val;
                    case {'Massfraction'}
                        val = sscanf(tmp1{2},'%e');
                        opt.('Oside').(remain(2:end)) = val;
                end
                tline = fgetl(fid);
            end
    end
end
%n
names= {};
data = {};
tline = fgetl(fid);
while ( ~strcmp(tline,'trailer') )
    names{end+1} = tline;
    data{end+1} = fscanf(fid,'%e',[ 1 n]);
    tline = fgetl(fid);
    tline = fgetl(fid);
end
fclose(fid);
end
%% helper function: convert list of variable names to the name-index map
function mapIndex = dictName(names)
tmp = nameTrans(names);
nvar = length(tmp);
mapIndex = containers.Map(tmp,1:nvar);
end
%% helper function: abreviate variable names
function out = nameTrans(names)
nvar = length(names);
out = cell(1,nvar);
for i = 1:nvar
    %remove whitespace
    tmp = strtrim(names{i});
    %remove units
    tmp = regexprep(tmp,'\s*\[.*\]\s*$','');
    %change massfraction
    tmp = regexprep(tmp,'massfraction-','Y_');
    tmp = regexprep(tmp,'massfraction','Y_');
    %change molefraction
    tmp = regexprep(tmp,'molefraction-','X_');
    tmp = regexprep(tmp,'molefraction','X_');
    %change ProdRate
    tmp = regexprep(tmp,'^ProdRate-','W_');
    tmp = regexprep(tmp,'^ProdRate','W_');
    %sanitize
    tmp = regexprep(tmp,'[\+,-,\*,\/]','_');
    switch(tmp)
        case {'temperature'}
            out{i} = 'T';
        case {'density'}
            out{i} = 'rho';
        otherwise
            out{i} = tmp;
    end
end
end
