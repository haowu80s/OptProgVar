function [ PROG ] = post_proc_B( B,  data_props, tol)
%POST_PROC_B post processes weight vector B to form the progress variable
%coefficient
B = B/norm(B,1);
B(abs(B)<= tol)=0;
phiNames = data_props.phiNames;
for i=1:data_props.nPhi
    if B(i)==0
        phiNames{i} = ' ';
    else
        phiNames{i} = strrep(phiNames{i}, 'Y_', '');
    end
end
%%
yScale = data_props.yScale;
yScale(yScale==0) = inf;
PROG = B./yScale;
PROG = PROG/norm(PROG,1);
figure(1)
bar(PROG);
set(gca, 'XTickLabel',phiNames, 'XTick',1:numel(phiNames))
ylabel('PROG Coeff');
end

