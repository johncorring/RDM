beta = [0:.25:2];
figure;

for k = 1:numel(beta)
    [nrg aff] = getCorrsTransformation(X, X, Y, beta(k));
    cParams = [aff; nrg];
    cConf.template = X;
    cConf.ctrl_pts = X;
    [transX transNu] = transformPointset(reshape(cParams', numel(cParams), 1)', cConf);
    cla;
    scatter3(X(:,1), X(:,2), X(:,3), 'rs');
    hold on;
    scatter3(Y(:,1), Y(:,2), Y(:,3), 'bx');
    hold on
    scatter3(transX(:,1), transX(:,2), transX(:,3), 'g*');
    drawnow;
    
    
end