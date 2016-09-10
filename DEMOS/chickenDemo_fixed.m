%get target ptset//this ptset is ordered and comes from extracting the 
%boundary from the MPEG7 silhouette images
load('MPEG7/CHICKENS/curves/ckshapes.mat');

x1 = ckshapes{8}(1:5:end,:);
x = [];
%this ensures that the ptset doesn't have repeated locations (bad for
%normal estimation)
for i = 2:size(x1,1)
    if(x1(i,1) ~= x1(i-1,1) | x1(i,2) ~= x1(i-1,2))
        x = [x; x1(i,:)];
    end
    
end
x1 = x;
[n1 d] = size(x1);
%center and scale
x1 = x1 - repmat(mean(x1), n1, 1);
x1 = x1/(max(max(abs(x1))));

%TODO:add some random rotation here... 
theta = pi/4*(rand(1)-.5);
R1 = [cos(theta) sin(theta); -sin(theta) cos(theta)];
%note the reflection of the x-axis // some chickens face +, some face -
x1 = (x1*R1)*([-1 0; 0 1]) + 0*repmat(randn(1,2),n1,1);
v1 = estimateNormals(x1, 2, 1);
x0 = x1;
[n1 d] = size(x1);

v1 = cos(kron(v1+pi,ones(2,1)) - kron(ones(n1,1), pi/2*[0; 1]));


%load a template chicken
x2 = ckshapes{10}(5:5:end,:);

[n2 d] = size(x2);
x2 = x2 - repmat(mean(x2), n2, 1);
x2 = x2/max(max(abs(x2)));

v2 = estimateNormals(x2, 2, 1);
v2 = cos(kron(v2,ones(2,1)) - kron(ones(n2,1), pi/2*[0;1]));

%here I've zapped the normals, because we're going to estimate them
v1 = v1';
v2 = v2';

%sigma and lambda are the main operating parameters for RDM, sigma controls
%the spatial sensitivity (lower is more localized) and lambda controls the
%frequency sensitivity (again, lower is more localized). If the ptsets are
%far apart or are expected to differ by a large deformation then larger
%values are more appropriate to start
sigma = .15;
lambda = .45;
%these are the annealing params, see the bottom for more on them
sdecay = .6;
ldecay = .6;

%some vis stuff
figure;
subplot(2,2,1);
plot(x2(:,1), x2(:,2), '-bx', 'LineWidth', 3);
hold on; plot(x1(:,1), x1(:,2), '-rs', 'LineWidth', 3);
title('Initial', 'FontSize', 20);
set(gca,'FontSize',20)
[xctrl0 yctrl0] = meshgrid(-1.5:.15:1.5);
hold on; mesh(xctrl0, yctrl0, 0*xctrl0, 'facecolor', 'none', 'edgecolor', 'k');

set(gca, 'XLim', [-1.5 1.5]);
set(gca, 'YLim', [-1.5 1.5]);


%create a configuration struct. This contains all of the info that the
%algorithm needs to run
config = rdmConfig(x1, x2, v1', v2', 'grbf', 'fixed', [''], sigma, lambda, 1);
%beta is a parameter that penalizes the bending energy of the
%transformation, depending on the class (grbf or tps) and the control
%points (currently a 5 x 5 grid) beta has slightly different operating
%values
%TODO: explore different beta values/ see compute_kernel2.m and
%rdmL2Bregmann_Free.m
config.beta = .005;
config.max_iter = 500;
%this stores the whole optimization history so that we can track the
%movement of the template. useful for debugging
opthist = [];


for i = 1:3
    
    %call to the algorithm. should return a runtime and output a series of
    %optimization values
    [p x2p v2p history config2] = rdmRegister(config);
    
    %storing the new (transformed) template normals
    v2i(i,:) = v2p;
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
        [xctrl0 yctrl0] = meshgrid(-1.5:.15:1.5);
        twarp0 = [xctrl0(:) yctrl0(:)];
        twarp = twarp0;
        pk = history.x(end,:);
        %transform grid
        [twarp] = transform_pointset(twarp, 'grbf-free', ... 
            [pk], config.ctrl_pts);
        [nctrl mctrl] = size(twarp);
        
        %(I made it square above)
        xctrl1 = reshape(twarp(:,1),  [sqrt(nctrl) sqrt(nctrl)]);
        yctrl1 = reshape(twarp(:,2),  [sqrt(nctrl) sqrt(nctrl)]);
        
        %plotting the grid and updated template, target normals
        subplot(2,2,i+1);
        quiver(x1(:,1), x1(:,2), v1p(1:2:end)', v1p(2:2:end)', 'rs', 'LineWidth', 3);
        hold on; plot(x2p(:,1), x2p(:,2), '-bx', 'LineWidth', 3);
        hold on; mesh(xctrl1, yctrl1, 0*xctrl1, 'facecolor', 'none', 'edgecolor', 'k');
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