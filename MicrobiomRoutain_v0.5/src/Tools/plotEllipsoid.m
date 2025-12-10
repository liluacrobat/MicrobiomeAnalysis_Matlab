function plotEllipsoid(X,Color)
%DEMOINERTIAELLIPSOID Demo program for the use of ellipsoids
%
%   Usage:
%   demoInertiaEllipsoid;
%
%   Example
%   demoInertiaEllipsoid
%
%   See also
%
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2011-06-21,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.
% Generate gaussian 3D data
if nargin<2
    Color='r';
end
points=X';
N=size(points,1);
Center=mean(points,1);
dis=pdist2(points,Center,'mahalanobis');

[~,idx]=sort(dis,'descend');
sel_idx=dis<dis(idx(round(N*0.33)));
points_filtered=points(sel_idx==1,:);

% display data

% axis equal;

% Fit a 3D inertia ellipsoid to data
elli = inertiaEllipsoid(points_filtered);

% draw the ellipsoid with transparency
drawEllipsoid(elli, 'FaceColor', Color, 'FaceAlpha', .3,'drawAxes', true,'drawEllipses', true, 'EllipseColor', 'k');
plot3(elli(1),elli(2),elli(3),'sk','MarkerSize',8);
for i=1:N
    plot3([points(i,1) elli(1)],[points(i,2) elli(2)],[points(i,3) elli(3)],'-k');
end
end