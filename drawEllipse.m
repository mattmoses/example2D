function drawEllipse( a,b,phi,x,y )
%DRAWELLIPSE 
%   example usage: drawEllipse(a,b,phi,x,y)

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

N = 80;

theta = linspace(0,2*pi * (N-1)/N,N);

coords = [x * ones(1,length(theta)); y * ones(1,length(theta))] + [cos(phi) , -sin(phi); sin(phi), cos(phi)] * [a*cos(theta); b*sin(theta)]; 

line(coords(1,:),coords(2,:),'Color','black');
line([coords(1,end) coords(1,1)],[coords(2,end) coords(2,1)],'Color','black'); 

end

