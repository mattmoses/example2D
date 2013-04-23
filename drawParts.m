function drawParts(g_ta, g_gr, rp, rh, d, w_gr, w_ta, dch)
% DRAWPARTS
% Draw target and grasped parts for probabilistic assembly model.
%
% g_ta 3x3 2D transformation to target
% g_gr 3x3 2D transformation to gripper
% rp radius of pin on target
% rh radius of hole on gripper
% d distance between pins
% w_gr is side length of gripper
% w_ta is side length of target
% dch set to 1 to draw a crosshair symbol for debug

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

% here are coords for origin of frame at center of part
p(:,1) = [-d/2; -d/2; 1];
p(:,2) = [-d/2; d/2; 1];
p(:,3) = [d/2; d/2; 1];
p(:,4) = [d/2; -d/2; 1];

% here are coords for the corners of the base part
% and top part. Note: base part (target) has the pins, top part (gripped)
% has the holes.
out_ta(:,1) = [-w_ta/2; -w_ta/2; 1];
out_ta(:,2) = [-w_ta/2; w_ta/2; 1];
out_ta(:,3) = [w_ta/2; w_ta/2; 1];
out_ta(:,4) = [w_ta/2; -w_ta/2; 1];

out_gr(:,1) = [-w_gr/2; -w_gr/2; 1];
out_gr(:,2) = [-w_gr/2; w_gr/2; 1];
out_gr(:,3) = [w_gr/2; w_gr/2; 1];
out_gr(:,4) = [w_gr/2; -w_gr/2; 1];

% coordinates for crosshairs for debug
chc = g_gr * 1/2*[-w_gr w_gr 0 0; 0 0 -w_gr w_gr; 2 2 2 2];

% transformed points
np_ta = g_ta * p; % pins in the target
out_ta = g_ta * out_ta; % corners of the base
np_gr = g_gr * p; % holes in the gripped part
out_gr = g_gr * out_gr; % corners of the gripped part

for i = 1:4
    targPin = drawCircle(np_ta(1,i), np_ta(2,i), rp, 30);    
    targOutline = line([out_ta(1,:), out_ta(1,1)],[out_ta(2,:), out_ta(2,1)]);
    gripHol = drawCircle(np_gr(1,i), np_gr(2,i), rh, 30);
    gripOutline = line([out_gr(1,:), out_gr(1,1)],[out_gr(2,:), out_gr(2,1)]);
    set(gripOutline, 'Color', 'k', 'LineWidth', 1);
    set(gripHol, 'Color', 'k', 'LineWidth', 1);
    set(targOutline, 'Color', 'k', 'LineWidth', 2);
    set(targPin, 'Color', 'k', 'LineWidth', 2);
    if dch
        line(chc(1,1:2), chc(2,1:2));
        line(chc(1,3:4), chc(2,3:4));
    end
end

