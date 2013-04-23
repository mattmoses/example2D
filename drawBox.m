function [ handle ] = drawBox( xx,r )
%DRAWBOX(x,r) draws a 3D box centered at x with size r
%   x and r are 3x1

% Version 1.0 Initial Release - April 2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
%  Copyright 2013 Matt Moses   
%  mmoses152 at gmail dot com
% 
%  This file is part of example2D. 
%
%  example2D is free software: you can redistribute it and/or modify
%  it under the terms of the GNU Lesser General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  example2D is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU Lesser General Public License for more details.
%
%  You should have received a copy of the GNU Lesser General Public License
%  along with example2D.  If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for face1
x1 = r(1)*[1 1 1 1]';
y1 = r(2)*[-1 1 1 -1]';
z1 = r(3)*[-1 -1 1 1]';

% for face 2 opposite face 1
x2 = -x1;
y2 = y1;
z2 = z1;

% for face 3
x3 = r(1)*[1 1 -1 -1]';
y3 = r(2)*[-1 1 1 -1]';
z3 = r(3)*[1 1 1 1]';

% for face 4
x4 = x3;
y4 = y3;
z4 = -z3;

% for face 5
x5 = r(1)*[1 1 -1 -1]';
y5 = r(2)*[-1 -1 -1 -1]';
z5 = r(3)*[-1 1 1 -1]';

% for face 6
x6 = x5;
y6 = -y5;
z6 = z5;

handle = patch([x1,x2,x3,x4,x5,x6]+xx(1),[y1,y2,y3,y4,y5,y6]+xx(2),[z1,z2,z3,z4,z5,z6]+xx(3),[0.2,0.6,0.8]);

end

