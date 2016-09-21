for j = 1:20
curve = alignedSh{j};
curve0 = alignedSh{1};
cla;
plot( curve0(:,1), curve0(:,2), 'b');
hold on;
plot( curve(:,1), curve(:,2), 'r');
title(int2str(j));
pause
end
 
figure;
for j = 10:20
 
    psihf = (1 + abs(batsdfs{j})/50).^(-2) .* exp( sqrt(-1)*batsdfs{j}/20);
    subplot(3,2,1);
    imagesc(abs(psihf).^2);
    colormap(gca, 'hot');
    title('|\Psi_{C}|^2', 'FontSize', 14);
    set(gca, 'YDir', 'normal');
    axis off;
    subplot(3,2,2);
    imagesc(20*angle(psihf));
    hold on; plot(batshapes{j}(:,1), batshapes{j}(:,2), 'b', 'LineWidth', 2);
    colormap(gca, 'hot');
    set(gca, 'YDir', 'normal');
set(gca, 'Clim', [-50 200]);
 
 
    axis off;
    title('\Theta_{C}', 'FontSize', 14);
     
     
    psig = exp(- batsdfs{j}.^2/350) .* exp( sqrt(-1)*batsdfs{j}/20);
    subplot(3,2,3);
    imagesc(abs(psig).^2);
    colormap(gca, 'hot');
    set(gca, 'YDir', 'normal');
    title('|\Psi_{G}|^2', 'FontSize', 14);
    axis off;
    subplot(3,2,4);
    imagesc(20*angle(psig));
    hold on; plot(batshapes{j}(:,1), batshapes{j}(:,2), 'b', 'LineWidth', 2);
    colormap(gca, 'hot');
    set(gca, 'YDir', 'normal');
set(gca, 'Clim', [-50 200]);
 
 
    title('\Theta_{G}', 'FontSize', 14);
    axis off;
     
    psig = exp(- max(batsdfs{j},0).^2/350) .* exp( sqrt(-1)*batsdfs{j}/20);
    subplot(3,2,5);
    imagesc(abs(psig).^2);
    colormap(gca, 'hot');
    set(gca, 'YDir', 'normal');
    title('|{\Psi^{\cap}_{g}}|^2', 'FontSize', 14);
    axis off;
    subplot(3,2,6);
    imagesc(20*angle(psig));
    hold on; plot(batshapes{j}(:,1), batshapes{j}(:,2), 'b', 'LineWidth', 2);
    colormap(gca, 'hot');
set(gca, 'Clim', [-50 200]);
 
 
    set(gca, 'YDir', 'normal');
    title('{\Theta^{\cap}_{g}}', 'FontSize', 14);
    axis off;
     
    pause;
end
 
figure;
colormap(gca, 'hot');
axis off
imagesc(batsdfs{10});
axis off
set(gca, 'YDir', 'normal');
title('SDF', 'FontSize', 14);
hold on; plot(batshapes{j}(:,1), batshapes{j}(:,2), 'b', 'LineWidth', 2);
set(gca, 'Clim', [-50 200]);
 
 
cdata = flipud(padarray(cdata, [2 2]));
cdata = imresize(cdata,[256 256])/max(max(cdata));
imX1 = zeros(size(cdata)+100);
imX1(51:end-50,51:end-50) = cdata;
cdata = flipud(padarray(cdata, [2 2]));
cdata = imresize(cdata,[256 256])/max(max(cdata));
imX2 = zeros(size(cdata)+100);
imX2(51:end-50,51:end-50) = cdata;
imX31 = [imX1 imX2-imX2];
imX32 = [imX1-imX1 imX2];
sdf1 = bwdist(imX31);
sdf1 = double(sdf1) - double(bwdist(1-imX31));
sdf2 = bwdist(imX32);
sdf2 = double(sdf2) - double(bwdist(1-imX32));
 
figure; subplot(3,3,1); imagesc(sdf1);
hold on; contour(sdf1, [0 .1], 'LineWidth', 2, 'Color', 'b');
colormap(gca, 'hot');
set(gca, 'CLim', [-200 300]);
set(gca, 'YDir', 'normal');
title('SDF of a Chicken', 'FontSize', 14);
axis off;
subplot(3,3,2); imagesc(sdf2);
hold on; contour(sdf2, [0 .1], 'LineWidth', 2, 'Color', 'b');
set(gca, 'YDir', 'normal');
set(gca, 'CLim', [-200 300]);
title('SDF of another Chicken', 'FontSize', 14);
axis off;
%for min//not adding
msdf = min(sdf1,sdf2);
subplot(3,3,3); imagesc(msdf);
hold on; contour(msdf, 'LineWidth', 2, 'Color', 'b');
set(gca, 'YDir', 'normal');
set(gca, 'CLim', [-200 300]);
title('Sum of SDFs', 'FontSize', 14);
axis off;
 
psi1 = exp( - sdf1.^2/250 + sqrt(-1)*sdf1/25);
 
psi2 = exp( - sdf2.^2/250 + sqrt(-1)*sdf2/25);
 
subplot(3,3,4); imagesc(25*angle(psi1));
hold on; contour(angle(psi1), [0 .1], 'LineWidth', 2, 'Color', 'b');
colormap(gca, 'hot');
set(gca, 'CLim', [-200 300]);
set(gca, 'YDir', 'normal');
title('G-SWF of aChicken', 'FontSize', 14);
axis off;
subplot(3,3,5); imagesc(25*angle(psi2));
hold on; contour(angle(psi2), [0 .1], 'LineWidth', 2, 'Color', 'b');
set(gca, 'YDir', 'normal');
set(gca, 'CLim', [-200 300]);
title('G-SWF of another Chicken', 'FontSize', 14);
axis off;
subplot(3,3,6); imagesc(25*angle(psi1+psi2));
hold on; contour(angle(psi1+psi2), [0 .01], 'LineWidth', 2, 'Color', 'b');
set(gca, 'YDir', 'normal');
set(gca, 'CLim', [-200 300]);
title('Superposition of G-SWF', 'FontSize', 14);
axis off;
 
 
psi11 = ( 2 + abs(sdf1)/80).^-2 .* exp(sqrt(-1)*sdf1/25);
 
psi22 = ( 2 + abs(sdf2)/80).^-2 .* exp(sqrt(-1)*sdf2/25);
 
 
subplot(3,3,7); imagesc(25*angle(psi11));
hold on; contour(angle(psi11), [0 .1], 'LineWidth', 2, 'Color', 'b');
colormap(gca, 'hot');
set(gca, 'CLim', [-200 300]);
set(gca, 'YDir', 'normal');
title('C-SWF of aChicken', 'FontSize', 14);
axis off;
subplot(3,3,8); imagesc(25*angle(psi22));
hold on; contour(angle(psi22), [0 .1], 'LineWidth', 2, 'Color', 'b');
set(gca, 'YDir', 'normal');
set(gca, 'CLim', [-200 300]);
title('C-SWF of another Chicken', 'FontSize', 14);
axis off;
subplot(3,3,9); imagesc(25*angle(psi11+psi22));
hold on; contour(angle(psi11+psi22), [0 .01], 'LineWidth', 2, 'Color', 'b');
set(gca, 'YDir', 'normal');
set(gca, 'CLim', [-200 300]);
title('Superposition of C-SWFs', 'FontSize', 14);
axis off;
 
%meshing
tSDF = batsdfs{10};
sh = batshapes{10};
sh2 = sh(1:2:end,:);
nu = estimate2dNormals(sh2,1,1);
nu = nu + pi;
figure; quiver(sh2(:,1), sh2(:,2), cos(nu), sin(nu));
%{
vert = [];
for i = 1:3:size(psi1,1);
for j = 1:3:size(psi1,2);
v1 = sub2ind([128 128], i,j);
v2 = sub2ind([128 128], i+1, j);
v3 = sub2ind([128 128], i+1, j+1);
vert = [[v1 v2 v3]; vert];
v1 = sub2ind([128 128], i,j);
v2 = sub2ind([128 128], i+1, j+1);
v3 = sub2ind([128 128], i, j+1);
vert = [[v1 v2 v3]; vert];
end
end
%}
[xpts ypts] = meshgrid(1:size(tSDF,1), 1:size(tSDF,2));
[psi1 phi1] = CWRXY(sh2, nu, xpts, ypts, 25, 30);
figure; imagesc(30*angle(psi1))
set(gca, 'YDir', 'normal');
%hold on; contour(50*angle(psi1), 0, 'Color', 'r');
%hold on; contour(50*angle(psi1), [0 .01], 'Color', 'r');
hold on; contour(30*angle(psi1), [-.5 0 .5], 'Color', 'b', 'LineWidth', 2);%, 'Color', 'k', 'LineWidth', 2);
hold on; contour(tSDF, [-.5 0 .5], 'Color', 'y', 'LineWidth', 2, 'LineStyle', ':');
 
axis off;
colormap(gca, 'hot');
set(gca, 'Position', [.1 .1 .8 .8]);
title('CWR vs. SDF', 'FontSize', 14);
set(gca, 'CLim', [-100 350]);
 
 
 
vert = [];
[m n] = size(psi1);
 
[Xmesh0 Ymesh0] = meshgrid(1:6:m, 1:6:n);
[p q] = size(Xmesh0);
for i = 1:p-1
    for j = 1:q-1
        v1 = sub2ind([p q], i,j);
        v2 = sub2ind([p q], i+1, j);
        v3 = sub2ind([p q], i+1, j+1);
        vert = [[v1 v2 v3]; vert];
        v1 = sub2ind([p q], i,j);
        v2 = sub2ind([p q], i+1, j+1);
        v3 = sub2ind([p q], i, j+1);
        vert = [[v1 v2 v3]; vert];
    end
end
 
pts0 = [Xmesh0(:) Ymesh0(:)];
[psi2 phi2] = meshCWRXY(sh2, nu, pts0, 55, 50);
figure; trisurf(vert, pts0(:,1), pts0(:,2), 50*angle(psi2), 'EdgeColor', 'none')
 
 
triZero0 = [];
lineZero0 = [];
for k = 1:size(vert)
    %for each triangle
   for j = 1:3
       %for each edge: check if it crosses 0; if so add the triangle to the
       %stack
       if(angle(psi2(vert(k,j))) < 0 & angle(psi2(vert(k,mod(j+1,3)+1))) > 0)
           triZero0 = [triZero0; k];
           %lineZero = [lineZero; j];
       end
   end
     
     
end
 
for i = 1:length(triZero0);
    %create line
    lnx = [pts0(vert(triZero0(i), 1),1) pts0(vert(triZero0(i), 2),1) pts0(vert(triZero0(i), 3),1) ...
        ];
    lny = [pts0(vert(triZero0(i), 1),2) pts0(vert(triZero0(i), 2),2) pts0(vert(triZero0(i), 3),2) ... 
        ];
    lnz = [angle(psi2(vert(triZero0(i), 1))) angle(psi2(vert(triZero0(i), 2))) ... 
        angle(psi2(vert(triZero0(i), 3)))];
     
    hold on; 
    trisurf([1 2 3], lnx, lny, lnz, 'FaceColor', 'b');
end
axis off;
colormap(gca, 'hot');
set(gca, 'Position', [.1 .1 .9 .8]);
title('CWR with ~ 900 mesh points', 'FontSize', 14);
set(gca, 'CLim', [-100 350]);
 
 
%try the cauchy
 
 
vert = [];
[m n] = size(tSDF);
 
[Xmesh0 Ymesh0] = meshgrid(1:3:m, 1:3:n);
[p q] = size(Xmesh0);
for i = 1:p-1
    for j = 1:q-1
        v1 = sub2ind([p q], i,j);
        v2 = sub2ind([p q], i+1, j);
        v3 = sub2ind([p q], i+1, j+1);
        vert = [[v1 v2 v3]; vert];
        v1 = sub2ind([p q], i,j);
        v2 = sub2ind([p q], i+1, j+1);
        v3 = sub2ind([p q], i, j+1);
        vert = [[v1 v2 v3]; vert];
    end
end
 
pts0 = [Xmesh0(:) Ymesh0(:)];
 
[psi4 phi4] = meshCauchyCWRXY(sh2, nu, pts0, 2, 50);
figure; trisurf(vert, pts0(:,1), pts0(:,2), 50*angle(psi4), 'EdgeColor', 'none')
 
triZero0 = [];
lineZero0 = [];
for k = 1:size(vert)
    %for each triangle
   for j = 1:3
       %for each edge: check if it crosses 0; if so add the triangle to the
       %stack
       if(angle(psi4(vert(k,j))) < 0 & angle(psi4(vert(k,mod(j+1,3)+1))) > 0)
           triZero0 = [triZero0; k];
           %lineZero = [lineZero; j];
       end
   end
     
     
end
 
for i = 1:length(triZero0);
    %create line
    lnx = [pts0(vert(triZero0(i), 1),1) pts0(vert(triZero0(i), 2),1) pts0(vert(triZero0(i), 3),1) ...
        ];
    lny = [pts0(vert(triZero0(i), 1),2) pts0(vert(triZero0(i), 2),2) pts0(vert(triZero0(i), 3),2) ... 
        ];
    lnz = [angle(psi4(vert(triZero0(i), 1))) angle(psi4(vert(triZero0(i), 2))) ... 
        angle(psi4(vert(triZero0(i), 3)))];
     
    hold on; 
    trisurf([1 2 3], lnx, lny, lnz, 'FaceColor', 'b');
end
axis off;
colormap(gca, 'hot');
set(gca, 'Position', [.1 .1 .9 .8]);
title('CWR with ~ 900 mesh points', 'FontSize', 14);
set(gca, 'CLim', [-100 350]);
%--end cauchy ex.
 
 
 
 
 
 
 
 
 
 
%create a random mesh
sh2 = sh(1:2:end,:);
nunew = estimate2dNormals(sh2, 1, 1);
nunew = nunew + pi;
[pts1] = rand(150,2)*355;
 
pts1 = [pts1; [1 1]; [1 355]; [355 1]; [355 355]; sh2; 
     sh+randn(size(sh))*2; sh+randn(size(sh))*2; sh+randn(size(sh))*2;];%... 
     %sh+randn(size(sh))*2; sh+randn(size(sh))*2; sh+randn(size(sh))*2;];
 
vert1 = delaunay(pts1);
tic
[psi3 phi3] = meshCWRXY(sh2, nunew, pts1, 25, 50);
toc
figure; trisurf(vert1, pts1(:,1), pts1(:,2), 50*angle(psi3), 'EdgeColor', 'none')
 
triZero = [];
lineZero = [];
for k = 1:size(vert1)
    %for each triangle
   for j = 1:3
       %for each edge: check if it crosses 0; if so add the triangle to the
       %stack
       if(angle(psi3(vert1(k,j))) < 0 & angle(psi3(vert1(k,mod(j+1,3)+1))) > 0)
           triZero = [triZero; k];
           %lineZero = [lineZero; j];
       end
   end
     
     
end
 
for i = 1:length(triZero);
    %create line
    lnx = [pts1(vert1(triZero(i), 1),1) pts1(vert1(triZero(i), 2),1) pts1(vert1(triZero(i), 3),1) ...
        ];
    lny = [pts1(vert1(triZero(i), 1),2) pts1(vert1(triZero(i), 2),2) pts1(vert1(triZero(i), 3),2) ... 
        ];
    lnz = [angle(psi3(vert1(triZero(i), 1))) angle(psi3(vert1(triZero(i), 2))) angle(psi3(vert1(triZero(i), 3))) ... 
         ];
     
    hold on; 
    trisurf([1 2 3], lnx, lny, lnz, 'FaceColor', 'b');
end
axis off;
colormap(gca, 'hot');
set(gca, 'Position', [.1 .1 .9 .8]);
title('CWR with ~ 1600 mesh, 700 source points', 'FontSize', 14);
set(gca, 'CLim', [-100 250]);
 
%build graph matrix for our algo
A = zeros(size(pts1,1), size(pts1,1));
for k = 1:size(vert1)
    %for each triangle
   for j = 1:3
       %for each edge: check if it crosses 0; if so add the triangle to the
       %stack
       if(50*angle(psi3(vert1(k,j))) < 1.5 & 50*angle(psi3(vert1(k,j))) > -1.5)
           %if this vertex is in a this strip
           %connect it with its neighbors, penalized by the change in the
           %angle
           A(vert1(k,j),vert1(k,mod(j+1,3)+1)) = (50*angle(psi3(vert1(k,j))) - ...
               50*angle(psi3(vert1(k,mod(j+1,3)+1))))^2 + (50*angle(psi3(vert1(k,j)))).^2;
       end
   end
     
     
end
A = sparse(A);
tic
 
 
dcirc = sqrt(sum((sh2 - circshift(sh2,+1)).^2,2));
 
%dsub = dcirc;
%truegDists = cumsum(dsub(end:-1:1));
%truegDists = truegDists(end:-1:1);
 
gDistEst = zeros(100,1);
truegDists = zeros(100,1);
 
figure;
firstinds = randperm(100);
secondinds = randperm(100);
B = diag(dcirc(2:end), -1);
B(1,end) = dcirc(end);
B = B + B';
B = sparse(B);
for k = 1:100
 
[Dists, path] = dijkstra(A, firstinds(k)+154, secondinds(k)+154);
vpath = pts1(path,:);
vpdist = sum(sqrt(sum((vpath - circshift(vpath,+1)).^2,2)));
gDistEst(k) = vpdist;
 
[truegd path2] = dijkstra(B, firstinds(k), secondinds(k));
truegDists(k) = truegd;
 
 
%{
cla;
lix = [pts1(path,1) pts1(path,1)];
liy = [pts1(path,2) pts1(path,2)];
plot(lix, liy, 'g', 'LineWidth', 3);
 
hold on; 
lix2 = [sh2(path2,1) sh2(path2,1)];
liy2 = [sh2(path2,2) sh2(path2,2)];
plot(lix2, liy2,  'r', 'LineWidth', 3);
pause;
%}
%toc;
%{
%}
end
figure; plot(truegDists);
hold on; plot(gDistEst, 'r')
 
%--end random mesh example
 
 
%3d mesh examples
sourceCat = 0;
sourceFileStr = strcat(strcat('tr_reg_0', strcat(int2str(sourceCat),int2str(0))),'.ply');
[X TRIV] = read_ply(sourceFileStr);%load(strcat(strcat('~/Dropbox/RDM/SHAPE_DATA/3d/nonrigid3d/gorilla', int2str(sourceCat)),'.mat'))
surface.X = X(:,1);
surface.Y = X(:,2);
surface.Z = X(:,3);
surface.TRIV = TRIV;
[m d] = size(surface.TRIV);
for i = 1:m
inds = surface.TRIV(i,:);
x1 = [surface.X(inds(1)) surface.Y(inds(1)) surface.Z(inds(1))];
x2 = [surface.X(inds(2)) surface.Y(inds(2)) surface.Z(inds(2))];
x3 = [surface.X(inds(3)) surface.Y(inds(3)) surface.Z(inds(3))];
xSource(i,:) = (x1 + x2 + x3)/3;
nu = cross(x1-xSource(i,:), x2-xSource(i,:));
xSourceNrms(i,:) = nu/norm(nu);
[vx vy vz] = cart2sph(xSourceNrms(i,1), xSourceNrms(i,2), xSourceNrms(i,3));
nuSource(i,:) = [vx vy];
end
 
nu1 = xSourceNrms(1:10:end,:);
sh1 = xSource(1:10:end,:);
 
vert3d = [];%(rand(1,3)-.5)*2;
vert3d = [vert3d; sh1; sh1+randn(size(sh1))*.05; sh1+randn(size(sh1))*.1;];% xSource-xSourceNrms*.1; xSource+xSourceNrms*.1;];%
[psi1 phi1] = meshCWRXY(sh1, nu1, vert3d, .003, 1);
angs = 1*angle(psi1);
% 
% [X Y Z] = meshgrid(-1: .05 :1);
% [newang] = TriScatteredInterp(vert3d(:,1), vert3d(:,2), vert3d(:,3), angle(psi1)); 
% 
% for i = 1:size(X,1); 
%     for j = 1:size(X,2);
%         for k = 1:size(X,3);
%             psi3(i,j,k) = newang(X(i,j,k), Y(i,j,k), Z(i,j,k));
%         end
%     end
% end
 
%figure;
%scatter3(vert3d(:,1), vert3d(:,2), vert3d(:,3), 1, angle(psi1), 'filled');
%axis off;
 
range = linspace((.001), (max(angs)-min(angs)),20);
range = range + min(angs);
figure;
for j = 1:19
   %get the vertices in this shell
   whj = find( range(j+1) > angs & angs > range(j));
   hold on;
   scatter3(vert3d(whj,1), vert3d(whj,2), vert3d(whj,3), 5,  angs(whj),'filled');
   set(gca, 'CLim', [range(1) range(end)]); 
   pause;
end
[dist path] = getGeodesicFromPsi(vert3d(1:size(sh1,1),:), angle(psi1(1:size(sh1,1))), 5, 300, .05, 10);
plotPath(vert3d(1:size(sh1),:), path);
%f1 = trisurf(tri, vert3d(:,1), vert3d(:,2), vert3d(:,3), angle(psi1));
 
%alpha(f1, .025);
 
 
 
%create his buddy
sourceCat = 0;
sourceFileStr = strcat(strcat('tr_reg_0', strcat(int2str(sourceCat),int2str(3))),'.ply');
[X TRIV] = read_ply(sourceFileStr);%load(strcat(strcat('~/Dropbox/RDM/SHAPE_DATA/3d/nonrigid3d/gorilla', int2str(sourceCat)),'.mat'))
surface.X = X(:,1);
surface.Y = X(:,2);
surface.Z = X(:,3);
surface.TRIV = TRIV;
[m d] = size(surface.TRIV);
for i = 1:m
inds = surface.TRIV(i,:);
x1 = [surface.X(inds(1)) surface.Y(inds(1)) surface.Z(inds(1))];
x2 = [surface.X(inds(2)) surface.Y(inds(2)) surface.Z(inds(2))];
x3 = [surface.X(inds(3)) surface.Y(inds(3)) surface.Z(inds(3))];
xSource2(i,:) = (x1 + x2 + x3)/3;
nu = cross(x1-xSource2(i,:), x2-xSource2(i,:));
xSourceNrms2(i,:) = nu/norm(nu);
%[vx vy vz] = cart2sph(xSourceNrms2(i,1), xSourceNrms2(i,2), xSourceNrms2(i,3));
%nuSource(i,:) = [vx vy];
end
sh2 = xSource2(1:10:end,:);
plotPath(sh2, path);
 
 
 
 
 
%try a 3d volume
sourceCat = 0;
sourceFileStr = strcat(strcat('tr_reg_0', strcat(int2str(sourceCat),int2str(3))),'.ply');
[X TRIV] = read_ply(sourceFileStr);%load(strcat(strcat('~/Dropbox/RDM/SHAPE_DATA/3d/nonrigid3d/gorilla', int2str(sourceCat)),'.mat'))
surface.X = X(:,1);
surface.Y = X(:,2);
surface.Z = X(:,3);
surface.TRIV = TRIV;
[m d] = size(surface.TRIV);
for i = 1:m
inds = surface.TRIV(i,:);
x1 = [surface.X(inds(1)) surface.Y(inds(1)) surface.Z(inds(1))];
x2 = [surface.X(inds(2)) surface.Y(inds(2)) surface.Z(inds(2))];
x3 = [surface.X(inds(3)) surface.Y(inds(3)) surface.Z(inds(3))];
xSource2(i,:) = (x1 + x2 + x3)/3;
nu = cross(x1-xSource2(i,:), x2-xSource2(i,:));
xSourceNrms2(i,:) = nu/norm(nu);
%[vx vy vz] = cart2sph(xSourceNrms2(i,1), xSourceNrms2(i,2), xSourceNrms2(i,3));
%nuSource(i,:) = [vx vy];
end
sh2 = xSource2(1:10:end,:);
plotPath(sh2, path);
 
%show how increasing the sampling rate improves geometry
 
%show how nonuniform sampling is ok
sh1 = xSource2(1:5:end,:);
nu2 = xSourceNrms2(1:5:end,:);
 
xmin = min(sh1(:,1));
ymin = min(sh1(:,2));
zmin = min(sh1(:,3));
xmax =max(sh1(:,1));
ymax =max(sh1(:,2));
zmax =max(sh1(:,3));
 
[X Y Z] = meshgrid(xmin:.01:xmax, ymin:.02:ymax, zmin:.01:zmax);
volCube = [X(:) Y(:) Z(:)]; [psi1 phi1] = meshCWRXY(sh1, nu2, volCube, .0007, 1);
psi2 = reshape(psi1, size(X));
[F1 V1] = MarchingCubes(X,Y,Z,angle(psi2),0);
%---end 3d mesh examples
 
%----build a solution to eikonal and make a movie
 
h1 = figure;
for k = 1:100
    
    acck = double(abs(tSDF) <= k);
    visk = double(abs(tSDF) <= k+2);
     
    subplot(1,2,1);
    imagesc( acck + 1 * visk);
    title('Visualization of FM');
    %set(gca, 'CLim', [0 2]);
    %caxis([0 1]);
    colormap(hot(3))
    colorbar( 'WestOutside','YTickLabels', {'Unvisited','Unvisited','Unvisited','Unvisited', ... 
        'Visited', 'Visited', 'Visited', 'Visited', ... 
        'Accepted', 'Accepted', 'Accepted', 'Accepted'}); %legend('Accepted', 'Visited');
    %set(gca, 'YDir', 'normal');
     
    subplot(1,2,2);
    imagesc(tSDF .* acck);
    title('Constructed SDF');
    %set(gca, 'YDir', 'normal');
    colormap(hot(256));
    pause(.25);
    F(k) = getframe(h1);
end
 
%---end movie for eikonal
 
 
%----build CWR and make a movie
 
h1 = figure;
pts = shapes{4};
nuangs = estimate2dNormals(pts,1,1);
[XX YY] = meshgrid(min(pts(:,1))-10:2:max(pts(:,1))+10,... 
min(pts(:,2))-10:2:max(pts(:,2))+10);
 
for k = 1:3:100
    
    ptsk = pts(1:100/k:end,:)
    nuangk = nuangs(1:100/k:end,:) + pi;
    psi1 = CWRXY(ptsk, nuangk, XX, YY, 1000/k, 50);
    subplot(1,2,1);
    %imagesc( acck + 1 * visk);
    quiver(ptsk(:,1), ptsk(:,2), cos(nuangk), sin(nuangk), 'LineWidth', 3);
    title('Visualization of Oriented Points');
    %set(gca, 'CLim', [0 2]);
    %caxis([0 1]);
    %colormap hot;
    %colorbar( 'WestOutside','YTickLabels', {'Unvisited','Unvisited','Unvisited','Unvisited', ... 
    %    'Visited', 'Visited', 'Visited', 'Visited', ... 
    %    'Accepted', 'Accepted', 'Accepted', 'Accepted'}); %legend('Accepted', 'Visited');
    %set(gca, 'YDir', 'normal');
     
    subplot(1,2,2);
    %imagesc(tSDF .* acck);
    imagesc(angle(psi1)*50);
    hold on; contour(angle(psi1)*50, [0 0], 'b');
    title('Constructed CWR');
    set(gca, 'YDir', 'normal');
    colormap hot
    pause(.25);
    F3(k) = getframe(h1);
end
 
%---end movie for eikonal
 
%%---build the cwr for an MPEG7 shape at varying levels of uniform sampling
%%and then compute a measure of accuracy of contour reconstruction by mean geodesic error 
%%from a curve recoverd via a level set to the sdf
 
load('SHAPES/BIRDS/curves/birdshapes.mat')
shapes = birdshapes;
 
 
 
numShapes = 20;
numSampleSpacing = 10;
numTrials = 50;
topRate = 50;
sigma = 200;
lamba = 50;
 
errMatMan = zeros(numShapes, numSampleSpacing, numTrials);
sampleRates = linspace(1,topRate,numSampleSpacing);
 
dists = pdist2(shapes{i});
 
ptDists = FastFloyd(dists);
 
for i = 1:numSamples
     
    ptsi = shapes{i};
    nuangs = estimate2dNormals(ptsi,1,1);
 
    %get the actual dists
    ptidistlist = sqrt(sum( (ptsi(2:end,:) - ptsi(1:end-1,:)).^2, 2));
    ptidists = diag(sqrt(sum( (ptsi(2:end,:) - ptsi(1:end-1,:)).^2, 2)), 1);
    ptidists = ptidists + ptidists';
    %ptDists = FastFloyd(dists);
     
    [X Y] = mesh(min(ptsi(:,1))-10:2:max(ptsi(:,1))+10, ... 
        min(ptsi(:,2))-10:2:max(ptsi(:,2))+10);
    vertsi = [X(:) Y(:)];
    for j = 1:numSampleSpacing
        sampleRatej = sampleRates(j);
        samplesj = [1:sampleRatej:size(ptsi,1)];
        ptsi_SampledAtRatej = ptsi(samplesj,:);
        nui_SampledAtRatej= nuangs(samplesj);
        psiij = meshCWRXY(ptsi_SampledAtRatej, nui_SampledAtRatej, vertsi, sigma, lambda);
         
         
        %TODO: curve comparison metrics
         
        for k = 1:numTrials
             
            %draw 2 random points
            rptk12 = randperm(samplesj)
             
            [estdistk path] = getGeodesicFromPsi(ptsi, ... 
                angle(psiij(1:size(sh1,1))), rptk12(1), rptk12(2), .05, 10);
 
            truedistk = sum(ptidistlist(min(rptk12(1:2)):max(rptsk12(1:2))));
            wrapdistk = sum(ptidistlist(max(rptk12(1:2)):end)) + ...
                sum(ptidistlist(1:min(rptk12(1:2))));
            truedistk = min(truedistk, wrapdistk);
             
        end
    end
end
 
 
%---