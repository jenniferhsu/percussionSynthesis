% mesh2d.m
%
% A Matlab script implementing a rectilinear waveguide mesh.  This
% scheme is based on the STK Mesh2D class.  The calculations make use
% of two sets of wave variable matrices that are alternated to save
% state.  Note that the mesh boundaries do not appear to be perfectly
% fixed in the plots.  This occurs because the plotted values do not
% include the outer boundary junctions.
%
% by Gary Scavone
% MUMT614, McGill University, 2004.
% 
% retrieved from: http://www.music.mcgill.ca/~gary/618/matlab/mesh2d.m
% edited by Jennifer Hsu (Winter 2016)
%
% for the structures, these are the port-ins/port-outs for:
% vxp: the right-going wave
% vxm: the left-going wave
% vyp: the down-going wave
% vym: the up-going wave
%
%

% Length of signal to calculate.
N = 3*44100;

% Sample rate for playback.
fs = 44100;

% Plot settings.
plotOpt = 1;
nDisplay = 100;   % Ticks between plots.
pTime = 0.1;   % Seconds to pause

% Number of (x,y) grid junctions (assuming square shape).
NJ = 21;

% Initialize calculation matrices.
vxp = zeros(NJ,NJ);
vxm = zeros(NJ,NJ);
vyp = zeros(NJ,NJ);
vym = zeros(NJ,NJ);

vxp1 = zeros(NJ,NJ);
vxm1 = zeros(NJ,NJ);
vyp1 = zeros(NJ,NJ);
vym1 = zeros(NJ,NJ);

v = zeros(NJ-1,NJ-1);
y = zeros(1,N);

% Provide initial velocity profile (width = 20% of grid)
dim = floor(NJ/5);
start = floor((NJ-dim)/2);
vals = zeros(NJ-1,1);
vals(start:start+dim-1) = 0.25*sin(pi*[0:dim-1]/(dim));
valm = vals*transpose(vals);
vxp1(1:NJ-1,1:NJ-1) = valm;
vyp1(1:NJ-1,1:NJ-1) = valm;
vxm1(1:NJ-1,2:NJ) = valm;
vym1(2:NJ,1:NJ-1) = valm;

% set all inputs to one junction to be 1
%vxp1(round(NJ/2), round(NJ/2)) = 1;
%vxp1((NJ+1)/2, (NJ+1)/2) = 1;
%vxp1(NJ-1, NJ-1) = 0.25;
%vyp1(NJ-1, NJ-1) = 0.25;
%vxm1(NJ-1, NJ) = 0.25;
%vym1(NJ, NJ-1) = 0.25;

reflCoeff = -0.9;

% Do initial display.
v = 0.5 * (vxp1(1:NJ-1,1:NJ-1) + vxm1(1:NJ-1,2:NJ) + ...
           vyp1(1:NJ-1,1:NJ-1) + vym1(2:NJ,1:NJ-1));
surf(v);
colormap('pink')
axis([0, NJ, 0, NJ, -0.12 0.12])
%zlabel('Velocity');
%disp('paused ...');
%pause;

% for a movie
vidfile = VideoWriter('meshMovie.mp4','MPEG-4');
vidfile.FrameRate = 7;
open(vidfile);

for n = 1:90%N

  if ( mod(n,2) == 0 ), % tick0

    v = 0.5 * (vxp(1:NJ-1,1:NJ-1) + vxm(1:NJ-1,2:NJ) + ...
               vyp(1:NJ-1,1:NJ-1) + vym(2:NJ,1:NJ-1));

    % Update outgoing junction wave components.
    vxp1(1:NJ-1,2:NJ) = v - vxm(1:NJ-1,2:NJ);     % portOutsLower
    vyp1(2:NJ,1:NJ-1) = v - vym(2:NJ,1:NJ-1);     % portOutsLeft
    vxm1(1:NJ-1,1:NJ-1) = v - vxp(1:NJ-1,1:NJ-1); % portOutsUpper
    vym1(1:NJ-1,1:NJ-1) = v - vyp(1:NJ-1,1:NJ-1); % portOutsRight

    % Do boundary reflections.
    vxp1(1:NJ-1,1) = reflCoeff*vxm(1:NJ-1,1);  
    vxm1(1:NJ-1,NJ) = reflCoeff*vxp(1:NJ-1,NJ);
    vyp1(1,1:NJ-1) = reflCoeff*vym(1,1:NJ-1);
    vym1(NJ,1:NJ-1) = reflCoeff*vyp(NJ,1:NJ-1);
  
  else % tick1

    v = 0.5 * (vxp1(1:NJ-1,1:NJ-1) + vxm1(1:NJ-1,2:NJ) + ...
               vyp1(1:NJ-1,1:NJ-1) + vym1(2:NJ,1:NJ-1));

    % Update outgoing junction wave components.
    vxp(1:NJ-1,2:NJ) = v - vxm1(1:NJ-1,2:NJ);
    vyp(2:NJ,1:NJ-1) = v - vym1(2:NJ,1:NJ-1);
    vxm(1:NJ-1,1:NJ-1) = v - vxp1(1:NJ-1,1:NJ-1);
    vym(1:NJ-1,1:NJ-1) = v - vyp1(1:NJ-1,1:NJ-1);

    % Do boundary reflections.
    vxp(1:NJ-1,1) = reflCoeff*vxm1(1:NJ-1,1);
    vxm(1:NJ-1,NJ) = reflCoeff*vxp1(1:NJ-1,NJ);
    vyp(1,1:NJ-1) = reflCoeff*vym1(1,1:NJ-1);
    vym(NJ,1:NJ-1) = reflCoeff*vyp1(NJ,1:NJ-1);
  
  end

  if plotOpt == 1
    %if ( mod(n-1,nDisplay) == 0 ),
        surf(v);
        %axis([0, NJ, 0, NJ, -0.15 0.15])
        axis([0, NJ, 0, NJ, -0.12 0.12])
        colormap('pink')
        %zlabel('Velocity');
        disp('paused ...');
        %pause(pTime);
        %pause;
        F(n) = getframe(gcf); 
        writeVideo(vidfile, F(n));
    %end
  end
  
  % Do output calculation.
  y(n) = v((NJ+1)/2, (NJ+1)/2);
  %pause;
  %keyboard
  %y(n) = v(NJ/2,NJ/2);

end

close(vidfile)

% scale it for saving
y = y/(max([max(y), abs(min(y))]));

plot(y);
%sound(0.99*y/max(abs(y)), fs);