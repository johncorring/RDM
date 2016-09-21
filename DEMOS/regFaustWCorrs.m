
F1 = '/3D/tr_reg_024.ply';
F2 = '/3D/tr_reg_034.ply';

[X triX] = read_ply(F1);
[Y triY] = read_ply(F2);

[mx dx] = size(triX);
[nx dx] = size(X);

[my dy] = size(triY);
[ny dy] = size(Y);

subRate = 20;
corrPtsInit = [5721 5211 5716 5437 2401 1591 2098 1561];

%get oriented points from the 3d pointsets
%first ptset
for i = 1:mx

    inds = triX(i,:);
    
    x1 = [X(inds(1),:)];
    x2 = [X(inds(2),:)];
    x3 = [X(inds(3),:)];
    
    xSource(i,:) = (x1 + x2 + x3)/3;
    nu = cross(x1-xSource(i,:), x2-xSource(i,:));
    
    xSourceNrms(i,:) = nu/norm(nu);
    [vx vy vz] = cart2sph(xSourceNrms(i,1), xSourceNrms(i,2), xSourceNrms(i,3));
    
    
end
 
nu1 = xSourceNrms(1:subRate:end,:);
sh1 = xSource(1:subRate:end,:);

%second ptset

for i = 1:my

    inds = triY(i,:);
    
    y1 = [Y(inds(1),:)];
    y2 = [Y(inds(2),:)];
    y3 = [Y(inds(3),:)];
    
    ySource(i,:) = (y1 + y2 + y3)/3;
    nu = cross(y1-ySource(i,:), y2-ySource(i,:));
    
    ySourceNrms(i,:) = nu/norm(nu);
    [vx vy vz] = cart2sph(ySourceNrms(i,1), ySourceNrms(i,2), ySourceNrms(i,3));
    
    
end
 
nu2 = ySourceNrms(1:subRate:end,:);
sh2 = ySource(1:subRate:end,:);

figure; 

subplot(2,2,1);
scatter3(sh1(:,1), sh1(:,2), sh1(:,3), 'ks');
hold on;
scatter3(sh2(:,1), sh2(:,2), sh2(:,3), 'rx');

title('Initial', 'FontSize', 20);
set(gca,'FontSize',20)

[xctrl0 yctrl0 zctrl0] = meshgrid(-.75:.25:.75, -2:.5:2,  -.5:.2:.5);
hold on; 
for k = 1:size(xctrl0,3)
    
    mesh(squeeze(xctrl0(:,:,k)), squeeze(yctrl0(:,:,k)), squeeze(zctrl0(:,:,k)), 'facecolor', 'none', 'edgecolor', 'k');
    hold on;
end


%sigma and lambda are the main operating parameters for RDM, sigma controls
%the spatial sensitivity (lower is more localized) and lambda controls the
%frequency sensitivity (again, lower is more localized). If the ptsets are
%far apart or are expected to differ by a large deformation then larger
%values are more appropriate to start
sigma = .08;
lambda = .3;
%these are the annealing params, see the bottom for more on them
sdecay = .6;
ldecay = .6;

tnu1 = nu1';
tnu2 = nu2';

config = rdmConfig(sh1, sh2, tnu1(:), tnu2(:), 'tps', 'free', [''], sigma, lambda, 1);
%beta is a parameter that penalizes the bending energy of the
%transformation, depending on the class (grbf or tps) and the control
%points (currently a 5 x 5 grid) beta has slightly different operating
%values
%TODO: explore different beta values/ see compute_kernel2.m and
%rdmL2Bregmann_Free.m
config.beta = .05;
config.max_iter = 500;
%this stores the whole optimization history so that we can track the
%movement of the template. useful for debugging
opthist = [];


for i = 1:3
    
    %call to the algorithm. should return a runtime and output a series of
    %optimization values
    [p x1p v1p history config2] = rdmRegister(config);
    
    %storing the new (transformed) template normals
    v1i(i,:) = v1p;
    %the updated configuation file is dumped back in technically could do
    %this in the LHS of call to rdmRegister
    config = config2;
    
    %when using 'free' the updated normal values are stored in the config
    %file. you can also fetch them from the end of p (the best parameters
    %for the objective, as found by fminunc)
    v1p = config.Vtarget';
    %if(1) since only 3 iters always viz//TODO: if you run this for more
    %iters then change the condition
    if(1)
        %create grid to viz the warp
        
        twarp0 = [xctrl0(:) yctrl0(:) zctrl0(:)];
        twarp = twarp0;
        pk = history.x(end,:);
        %transform grid
        [twarp] = transform_pointset(twarp, 'grbf-free', ... 
            [pk], config.ctrl_pts);
        [nctrl mctrl] = size(twarp);
        
        %(I made it square above)
        xctrl1 = reshape(twarp(:,1),  size(xctrl0));
        yctrl1 = reshape(twarp(:,2),  size(xctrl0));
        zctrl1 = reshape(twarp(:,3),  size(xctrl0));
        
        %plotting the grid and updated template, target normals
        subplot(2,2,i+1);
        scatter3(sh1(:,1), sh1(:,2), sh1(:,3), 'rs', 'LineWidth', 3);
        hold on; 
        scatter3(x1p(:,1), x1p(:,2), x1p(:,3), 'bx', 'LineWidth', 3);
        
        hold on; 
        for k = 1:size(xctrl1,3)

            mesh(squeeze(xctrl1(:,:,k)), squeeze(yctrl1(:,:,k)), squeeze(zctrl1(:,:,k)), 'facecolor', 'none', 'edgecolor', 'k');
            hold on;
        end
        title( strcat(strcat('After ',int2str(i)), ' iterations' ), 'FontSize', 20);
        set(gca,'FontSize',20);
        set(gca, 'XLim', [-1.5 1.5]);
        set(gca, 'YLim', [-1.5 1.5]);
        drawnow;
    end
    
    %this is the deterministic annealing part//TODO: tinker with the sdecay
    %and ldecay params. IDEA: when using 'free' lambda should start high 
    %and end roughly at lambda^2 <= sigma. Lambda should start high because
    %we don't want to be too sensitive to incorrect normal vals on the
    %target
    config.sigma = sigma * sdecay^i;
    config.lambda = lambda * ldecay^i;
end
