function [corrPtsX corrPtsY indices] = get3dCorrs(X, Y)


corrPtsX =  clickA3DPoint(X');
pause;

cla;

corrPtsY =  clickA3DPoint(Y');
pause;

indices = [];


end