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

imgDir = 'meshImages/'

%% run the mesh and save jpeg images

% Sample rate for playback.
fs = 22025;
dur = 0.0015;

% Length of signal to calculate.
N = floor(dur*fs);

% Plot settings.
nDisplay = 1;   % Ticks between plots.
pTime = 0.01;   % Seconds to pause

% Number of (x,y) grid junctions (assuming square shape).
NJ = 31;

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

% Do initial display.
v = 0.5 * (vxp1(1:NJ-1,1:NJ-1) + vxm1(1:NJ-1,2:NJ) + ...
           vyp1(1:NJ-1,1:NJ-1) + vym1(2:NJ,1:NJ-1));
surf(v);
axis([0, NJ, 0, NJ, -0.13 0.13])
zlabel('velocity');
set(gcf, 'Position', [0, 0, 800, 800], 'PaperUnits', 'points')
saveas(gcf, [imgDir 'jpeg/0'], 'jpeg')

disp('paused ...');
pause;



for n = 1:N

  if ( mod(n,2) == 0 ) % tick0

    v = 0.5 * (vxp(1:NJ-1,1:NJ-1) + vxm(1:NJ-1,2:NJ) + ...
               vyp(1:NJ-1,1:NJ-1) + vym(2:NJ,1:NJ-1));

    % Update outgoing junction wave components.
    vxp1(1:NJ-1,2:NJ) = v - vxm(1:NJ-1,2:NJ);
    vyp1(2:NJ,1:NJ-1) = v - vym(2:NJ,1:NJ-1);
    vxm1(1:NJ-1,1:NJ-1) = v - vxp(1:NJ-1,1:NJ-1);
    vym1(1:NJ-1,1:NJ-1) = v - vyp(1:NJ-1,1:NJ-1);

    % Do boundary reflections.
    vxp1(1:NJ-1,1) = -0.99*vxm(1:NJ-1,1);
    vxm1(1:NJ-1,NJ) = -0.99*vxp(1:NJ-1,NJ);
    vyp1(1,1:NJ-1) = -0.99*vym(1,1:NJ-1);
    vym1(NJ,1:NJ-1) = -0.99*vyp(NJ,1:NJ-1);
  
  else % tick1

    v = 0.5 * (vxp1(1:NJ-1,1:NJ-1) + vxm1(1:NJ-1,2:NJ) + ...
               vyp1(1:NJ-1,1:NJ-1) + vym1(2:NJ,1:NJ-1));

    % Update outgoing junction wave components.
    vxp(1:NJ-1,2:NJ) = v - vxm1(1:NJ-1,2:NJ);
    vyp(2:NJ,1:NJ-1) = v - vym1(2:NJ,1:NJ-1);
    vxm(1:NJ-1,1:NJ-1) = v - vxp1(1:NJ-1,1:NJ-1);
    vym(1:NJ-1,1:NJ-1) = v - vyp1(1:NJ-1,1:NJ-1);

    % Do boundary reflections.
    vxp(1:NJ-1,1) = -0.99*vxm1(1:NJ-1,1);
    vxm(1:NJ-1,NJ) = -0.99*vxp1(1:NJ-1,NJ);
    vyp(1,1:NJ-1) = -0.99*vym1(1,1:NJ-1);
    vym(NJ,1:NJ-1) = -0.99*vyp1(NJ,1:NJ-1);
  
  end

  if ( mod(n,nDisplay) == 0 )
    surf(v);
    axis([0, NJ, 0, NJ, -0.13 0.13])
    zlabel('velocity');
    
    %disp('paused ...');
    %pause(pTime);

    set(gcf, 'Position', [0, 0, 800, 800], 'PaperUnits', 'points')
    saveas(gcf, [imgDir 'jpeg/' num2str(n)], 'jpeg')
    
  end
  
  % Do output calculation.
  y(n) = v(NJ-1,NJ-1);

end

%plot(y);
%sound(0.99*y/max(abs(y)), fs);

%% take the jpeg images and put them into a gif file
% see this answer for help:
% https://www.mathworks.com/matlabcentral/answers/94495-how-can-i-create-animated-gif-images-in-matlab

imgList = dir([imgDir 'jpeg/']);
imgList = imgList(3:end); % get rid of the '.' and '..'

gifPath = [imgDir 'meshGIF.gif'];

for i=0:length(imgList)-1
    
    [imgDir 'jpeg/' num2str(i) '.jpg']
    img = imread([imgDir 'jpeg/' num2str(i) '.jpg'], 'jpeg');
    [imind, cm] = rgb2ind(img, 256); 
    if i == 1 
        imwrite(imind, cm, gifPath, 'gif', 'Loopcount', inf, 'DelayTime', 0.25); 
    else 
        imwrite(imind, cm, gifPath, 'gif', 'WriteMode', 'append', 'DelayTime', 0.25); 
    end 
    
end