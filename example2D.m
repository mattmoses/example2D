% example2D.m
%
% This code generates some plots for an example problem in probabilistic
% assembly

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

close all
clear all

% set the pseudo-random number generator seed
randn('seed',42);

% parameters for target and gripped parts
rpin = 0.2; % radius of pin on target
rhol = 0.375; % radius of hole in gripped part
d = 2.75; % distance between pins
w_gr = 4; % side length of gripped part
w_ta = 5; % side length of target part

L1 = 10; % length of link 1
L2 = 7.5; % length of link 2
usingrand = 0; % flag set to 0 if not using random joint angles
q1_0 = pi/4; % "commanded" joint angle 1
q2_0 = -pi/2; % "commanded" joint angle 2
q3_0 = pi/4; % "commanded" joint angle 3
axlimits2 = [12.1 12.7 1.3 2.3];

% % randomized joint angles for testing purposes
% q1_0 = 2*pi*rand(1);
% q2_0 = 2*pi*rand(1);
% q3_0 = 2*pi*rand(1);
% usingrand = 1;

var1 = 0.001; % variance of joint angle 1
var2 = 0.001; % variance of joint angle 2
var3 = 0.001; % variance of joint angle 3

% % ----------------------
% % As a sanity check, this part will show how the two pieces are mated 
% % together for one arbitrary transformation
xx = 0;
yy = 0;
phi = 0.065;
g = [cos(phi) -sin(phi) xx; sin(phi) cos(phi) yy; 0 0 1];
% alpha() will draw if draw enable is ~= 0
figure
alpha(g, rpin, rhol, d, w_gr, w_ta, 1);
axis equal

% ----------------------------
% This part draws the commanded frame and the ideal location of the parts.

% forward kinematics
c1 = cos(q1_0);
c2 = cos(q2_0);
c12 = cos(q1_0+q2_0);
s1 = sin(q1_0);
s2 = sin(q2_0);
s12 = sin(q1_0+q2_0);
x_com = L1*c1 + L2*c12;
y_com = L1*s1 + L2*s12;
phi_com = q1_0 + q2_0 + q3_0;

% draw lines to represent the SCARA robot arm
figure
hh1 = line([0 L1*c1],[0 L1*s1]);
hh2 = line([L1*c1 x_com],[L1*s1 y_com]);
set(hh1, 'Color', 'k');
set(hh2, 'Color', 'k');
hold on

% generate the gripper and target frames
g_gr = [cos(phi_com) -sin(phi_com) x_com; sin(phi_com) cos(phi_com) y_com; 0 0 1];
g_ta = g_gr;

% draw the target and gripped parts
drawCrosshair = 0; 
drawParts(g_ta, g_gr, rpin, rhol, d, w_gr, w_ta, drawCrosshair);
xlabel('x')
ylabel('y')
axis equal

% ------------------------------------------------------------------------
% this next part plots 2D contours of alpha with g_gr near g_ta
% the same g_ta as defined above is used

% Note that contour uses the same index convention as surf. From the
% help entry for SURF:
%     SURF(x,y,Z) and SURF(x,y,Z,C), with two vector arguments replacing
%     the first two matrix arguments, must have length(x) = n and
%     length(y) = m where [m,n] = size(Z).  In this case, the vertices
%     of the surface patches are the triples (x(j), y(i), Z(i,j)).
%     Note that x corresponds to the columns of Z and y corresponds to
%     the rows.

xylims = 0.3;
philims = 0.06;
nsa = 32; % number of samples along X and Y
nphi = 16; % number of angle samples
% Now make the sample grid. Coordinates are chosen for a range of g
% surrounding g_ta. 
xx = x_com + linspace(-xylims,xylims,nsa);
yy = y_com + linspace(-xylims,xylims,nsa);
[XX, YY] = meshgrid(xx,yy);
phiphi = phi_com + linspace(-philims, philims, nphi);
for k = 1:nphi
    k
    for i = 1:nsa
        for j = 1:nsa
            g = [cos(phiphi(k)) -sin(phiphi(k)) xx(i); sin(phiphi(k)) cos(phiphi(k)) yy(j); 0 0 1];
            asmCheck(i,j,k) = alpha(inv(g_ta)*g, rpin, rhol, d, w_gr, w_ta, 0);
            % uncomment this to directly view the nonzero samples of alpha
            % if (asmCheck(i,j,k)) plot(xx(i),yy(j),'k.'); end
        end
    end
    % if desired, plot the contour of alpha. note transpose.
    % [cc, hh] = contour(xx,yy,asmCheck(:,:,k)',0.5);
end


% ------------------------------------------------------------------------
% This part creates an ensemble of frames by randomly perturbing the
% commanded joint angles. A normal distribution is fit to the samples.
% Remember, for a random scalar x, var(a*x) = a^2 * var(x).
% For a random vector X, cov(A*X) = A * cov(X) * A'

% This is the commanded frame - where the robot would be if there was no
% perturbation.
frameL = 0.05;
figure
draw90s(x_com, y_com, phi_com, frameL, 3, 0);
hold on

% Now generate the random samples and fit the Gaussian.
nmansam = 1000; % number of samples for creating manipulator pdf
Cov = zeros(3,3); % initialize sample covariance matrix

% generate an ensemble of joint coordinates...
% Note: this generates joint angles that are centered around the 
% commanded position defined above, which coincides with g_ta, 
% loosely speaking: mean(g_gr) = g_ta. If you want to move the manipulator 
% slightly off so that mean(g_gr) ~= g_ta you can do that here by changing
% q1_0, q2_0, q3_0.
for i = 1:nmansam
    q1(i) = q1_0 + sqrt(var1) * randn(1);
    q2(i) = q2_0 + sqrt(var2) * randn(1);
    q3(i) = q3_0 + sqrt(var3) * randn(1);
end    

% do forward kinematics...
for i = 1:nmansam
    c1 = cos(q1(i));
    c2 = cos(q2(i));
    c12 = cos(q1(i)+q2(i));
    s1 = sin(q1(i));
    s2 = sin(q2(i));
    s12 = sin(q1(i)+q2(i));
    x(i) = L1*c1 + L2*c12;
    y(i) = L1*s1 + L2*s12;
    phi(i) = q1(i) + q2(i) + q3(i);
end

% draw the frames
draw90s3(x, y, phi, frameL, 1, 1);
xlabel('x')
ylabel('y')
axis equal

% find the sample means for each coordinate
x_sm = sum(x)/nmansam;
y_sm = sum(y)/nmansam;
phi_sm = sum(phi)/nmansam;

% draw a big frame at the sample mean for debug purposes
draw90s(x_sm, y_sm, phi_sm, frameL*5, 3, 1);

% find the sample covariance matrix. Note this is a vectorized version. Can
% do it with a for loop if you want...
Cov = Cov + [x - x_sm; y - y_sm; phi - phi_sm] * [x - x_sm; y - y_sm; phi - phi_sm]';
Cov = 1/(nmansam-1) * Cov; % note we use 1/(N-1) because we are using the sample mean

% for loop method for finding Cov
%Cov2 = 0
%for i = 1:nmansam
%    Cov2 = Cov2 + [x(i) - x_sm; y(i) - y_sm; phi(i) - phi_sm] * [x(i) - x_sm; y(i) - y_sm; phi(i) - phi_sm]';
%end
%Cov2 = 1/(nmansam-1) * Cov2;

% if we want to marginalize over phi, then we have 
Cov_xy = Cov(1:2,1:2);

% ok now we have a covariance matrix for the simulated end-effector data.
% we want to display this as an ellipsoid. diagonalizing Cov will give us
% the rotation matrix for rotating the ellipsoid
[Vcov, Dcov] = eigs(Cov);
[Vcov_xy, Dcov_xy] = eigs(Cov_xy);

% just in case Vcov is a reflection, this converts it to a plain rotation
Vcov = det(Vcov)*Vcov;
Vcov_xy = det(Vcov_xy)*Vcov_xy;

% generate an ellipsoid at the origin with radii equal to sqrt(eigenvalues) of
% Covariance matrix
Nep = 30; % dimension of X,Y,Z is Nep+1 x Nep+1
Sep = 1; % scale factor of ellipsoid
% This gives an ellipsoid centered on and aligned with the origin
[X,Y,Z]=ellipsoid(0,0,0, Sep*sqrt(Dcov(1,1)), Sep*sqrt(Dcov(2,2)), Sep*sqrt(Dcov(3,3)), Nep);

% now we reshape the X,Y,Z points, rotate and translate them, and reshape
% them back to a form surf and contour can plot.
temp = [reshape(X,(Nep+1)^2,1), reshape(Y,(Nep+1)^2,1), reshape(Z,(Nep+1)^2,1)];
temp = Vcov * temp';
X = reshape(temp(1,:)' + x_sm, Nep+1, Nep+1);
Y = reshape(temp(2,:)' + y_sm, Nep+1, Nep+1);
Z = reshape(temp(3,:)' + phi_sm, Nep+1, Nep+1);

% this shows a slice of the gaussian ellipsoid at phi_com
% for some reason contour doesn't draw well at phi_com = 0.
contour(X,Y,Z,phi_com+1E-06); 

% uncomment this to see the frame coordinates marginalized over phi
drawEllipse(Sep*sqrt(Dcov_xy(1,1)), Sep*sqrt(Dcov_xy(2,2)), acos(Vcov_xy(1,1)), x_sm, y_sm);

% ------------------------------------------------------------------------
% % This is just a double check of the covariance matrix.
% % This will plot a bunch of "virtual" sample points that can be
% % compared with the generated frames above.
% % first factor the covariance matrix as: Cov = Aco * Aco'
% Aco = Vcov * Dcov.^0.5;
% % now generate an ensemble using a standard normal distribution and
% % transform with Aco
% sampFrames = Aco * randn(3,nmansam) + repmat([x_com; y_com; phi_com],1,nmansam);
% draw90s3(sampFrames(1,:), sampFrames(2,:), sampFrames(3,:), frameL, 1, 0);
% plot(sampFrames(1,:),sampFrames(2,:),'k.');



% for i = 1:nphi
%     figure
%     %subplot(2,2,i);
%     [cc, hh] = contour(xx,yy,asmCheck(:,:,i)',0.5);
%     hold on
%     contour(X,Y,Z,phiphi(i))
%     xlabel('x')
%     ylabel('y')
%     if ~usingrand
%         axis equal
%         axis(axlimits2);
%         title(['\phi = ' num2str(phiphi(i))]);
%     end
% end

% ------------------------------------------------------------------------
% This part makes a 3D plot of rho and alpha. Alpha is shown approximated
% using the merge() function (see paper). Rho is shown as an ellipsoid
% representing the level set at a distance of Sep*(std) from the mean where
% std = sqrt(variance). (The radii of the ellipsoid are the square roots of
% the eigenvalues of the covariance, and the axes are aligned with the
% eigenvectors.)
figure
colormap gray
surf(X,Y,Z,0.5*ones(size(Z)))
hold on
rxx = xx(2)-xx(1);
ryy = yy(2)-yy(1);
rphi = phiphi(2)-phiphi(1);

% CC is the scale factor for the multivariate gaussian
CC = (2*pi)^(-3/2)*det(Cov)^(-1/2); 
% as a check, we integrate rho to make sure it comes out close to 1.0
% Note: In order for this summation to come out close to 1, the grid has to
% extend quite a ways out. If you are just integrating gamma, you only
% need to do the summation for alpha() that are nonzero.
% These are the parameters that yield a rho_int near 1:
% xylims = 2;
% philims = 0.2;
% nsa = 100; % number of samples along X and Y
% nphi = 32; % number of angle samples

rho_int = 0;
gamma = 0; % this is the number we have been looking for all along
iCov = inv(Cov); % inverse of covariance matrix

% this loop will draw alpha and calculate the summations
for i = 1:nsa
    i
    for j = 1:nsa
        for k = 1:nphi
            cpoint = [xx(i); yy(j); phiphi(k)] - [x_sm; y_sm; phi_sm];
            rho = CC * exp( -1/2 * cpoint' * iCov * cpoint );
            rho_int = rho_int + rho * rxx * ryy * rphi; 
            if asmCheck(i,j,k) 
                gamma = gamma + rho * rxx * ryy * rphi; 
                drawBox([xx(i), yy(j), phiphi(k)],[rxx/2.00001,ryy/2.00001,rphi/2.00001]);
            end
        end
    end
end

[x_com; y_com; phi_com]
[x_sm; y_sm; phi_sm]
Cov
rho_int
gamma

xlabel('x')
ylabel('y')
zlabel('\phi')
axis equal
camlight left

set(gca, 'DataAspectRatio', [2 2 1])
camlight(45, 1)
view(45,15)

% Note: here are the other values used to generate the figure rho_alpha_3D_1
% nsa = 32
% nphi = 16
% xx = x_com + linspace(-.3,.3,nsa);
% yy = y_com + linspace(-.3,.3,nsa);
% phiphi = phi_com + linspace(-0.06, 0.06, nphi);
% Sep = 1
% rpin = 0.2; % radius of pin on target
% rhol = 0.375; % radius of hole in gripped part
% d = 2.75; % distance between pins
% L1 = 10; % length of link 1
% L2 = 7.5; % length of link 2
% q1_0 = pi/4; % "commanded" joint angle 1
% q2_0 = -pi/2; % "commanded" joint angle 2
% q3_0 = pi/4; % "commanded" joint angle 3
% var1 = 0.001; % variance of joint angle 1
% var2 = 0.001; % variance of joint angle 2
% var3 = 0.001; % variance of joint angle 3
% x_com = 12.37436867076458
% y_com = 1.76776695296637
% phi_com = 0

% Here are the parameters used to generate the contour plots of alpha and
% rho.
% nsa = 200; % number of samples along X and Y
% nphi = 4; % number of angle samples
% xx = x_com + linspace(-.3,.3,nsa);
% yy = y_com + linspace(-.3,.3,nsa);
% phiphi = phi_com + linspace(-0.05, 0.05, nphi);

% The results will vary slightly from run to run because of the random
% samples. Cov, the sample means, rho_int, and gamma all have a little
% jitter from run to run.
% Typical values for gamma are 0.10 to 0.11


