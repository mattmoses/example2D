function draw90s(x,y,phi,L,Linewidth,color_en)
% DRAW90S
% draw90s(x,y,phi,L,Linewidth)
% draw some 90 degree "diad" frames at x,y with angle phi
% L is length of lines
% if color_en ~= 0 then the lines are colored, otherwise black

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

for i = 1:length(x)
    hx = line([x(i) x(i)+L*cos(phi(i))],[y(i) y(i)+L*sin(phi(i))]);
    hy = line([x(i) x(i)-L*sin(phi(i))],[y(i) y(i)+L*cos(phi(i))]);
    set(hx, 'Linewidth', Linewidth, 'Color', 'k');
    set(hy, 'Linewidth', Linewidth, 'Color', 'k');
    if color_en
        set(hx, 'Color', 'r');
        set(hy, 'Color', 'g');
    end
end