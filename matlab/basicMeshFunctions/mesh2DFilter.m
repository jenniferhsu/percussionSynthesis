function y = mesh2DFilter(ex, Nj, N, fs, inX, inY, outX, outY, B, A, plotOn)
%Y = mesh2DFilter(ex, Nj, N, fs, inX, inY, outX, outY, B, A, plotOn)
% this is a rectangular 2d mesh waveguide implementation that allows you
% to set the input excitation, specify the input/output points, and also
% change the filters at the boundaries. 
% inputs:
%   ex: the input excitation as an array
%   Nj: the number of junctions (x or y) assuming a square shape
%       this means that the square mesh will be of size NjxNj
%   N: the output signal length
%   fs: the sample rate
%   inX: input point x coordinate 
%   inY: input point y coordinate 
%   outX: listening point x coordinate 
%   outY: listening point y coordinate
%   B: filter coefficients (that multiply with the input)
%   A: filter coefficients (that multiply with the output)
%   plotOn: set to 1 to plot as an animation and 0 otherwise
% output:
%   y: the generated signal
%
% note: for inX, inY, outX and outY, the input signal is divided by 4 and
% input to the junction at (inX, inY) from all 4 directions. the output is 
% the velocity at junction (outX, outY)
%
% to call this function for a 10x10 square mesh:
%
% ex = 1;
% Nj = 10;
% N = 44100;
% fs = 44100;
% inX = 3;
% inY = 3;
% outX = 5;
% outY = 5;
% B = -[0.5 0.5];
% A = 1;
% plotOn = 0;
% y = mesh2DFilter(ex, Nj, N, fs, inX, inY, outX, outY, B, A, plotOn);
%

y = zeros(1,N);
Ns = Nj + 1;    % number of sections

% delay line structures
upper = zeros(Nj, Ns);
lower = zeros(Nj, Ns);
left = zeros(Ns, Nj);
right = zeros(Ns, Nj);

% scattering junction port structures
portInsUpper = zeros(Nj, Nj);
portInsLower = zeros(Nj, Nj);
portInsLeft = zeros(Nj, Nj);
portInsRight = zeros(Nj, Nj);
portOutsUpper = zeros(Nj, Nj);
portOutsLower = zeros(Nj, Nj);
portOutsLeft = zeros(Nj, Nj);
portOutsRight = zeros(Nj, Nj);

% velocity at the junctions
vj = zeros(Nj, Nj);

% for plotting, let's extend vj to have the boundaries
vjb = zeros(Nj+2, Nj+2);
            
           
% filters for the reflections - we'll use the simplest low pass filter
last_in = zeros(Nj,4); % 4 for each side that needs a reflection
%A = 1;


% gary scavone's initialization code makes it look much better
% if I use this, make sure to take away the input excitation?
% dim = floor(Nj/5);
% start = floor((Nj-dim)/2);
% vals = zeros(Nj-1,1);
% vals(start:start+dim-1) = 0.25*sin(pi*[0:dim-1]/(dim));
% valm = vals*transpose(vals);
% upper(1:Nj-1,1:Nj-1) = valm;
% lower(1:Nj-1,1:Nj-1) = valm;
% left(1:Nj-1,2:Nj) = valm;
% right(2:Nj,1:Nj-1) = valm;



% for each output sample
for n=1:N
    
    % get the input excitation value
    if length(ex) < n
        % do nothing
    else
        %upper(inX, inY) = upper(outX, outY) + ex(n); %<-- i think that's a
        %bug?
        upper(inX,inY) = upper(inX,inY) + (0.25*ex(n));
        lower(inX,inY+1) = lower(inX,inY+1) + (0.25*ex(n));
        left(inX+1,inY) = left(inX+1,inY) + (0.25*ex(n));
        right(inX,inY) = right(inX,inY) + (0.25*ex(n));
    end
    
    % save previous upper/lower/left/right rails
    upperOld = upper;
    lowerOld = lower;
    leftOld = left;
    rightOld = right;
    
    % get the port inputs
    portInsUpper = upperOld(:,1:Nj);
    portInsLower = lowerOld(:,2:Ns);
    portInsRight = rightOld(1:Nj,:);
    portInsLeft = leftOld(2:Ns,:);

    % calculate the velocity(?)
    vj = 0.5 * (portInsUpper + portInsLower + ...
        portInsLeft + portInsRight);

    % calculate the velocity at the port outputs
    portOutsUpper = vj - portInsLower;
    portOutsLower = vj - portInsUpper;
    portOutsRight = vj - portInsLeft; 
    portOutsLeft = vj - portInsRight; 
    
    % update upper/lower values
    upper(:,2:Ns) = portOutsUpper;
    lower(:,1:Nj) = portOutsLower;
    right(2:Ns,:) = portOutsRight;
    left(1:Nj,:) = portOutsLeft;
    
    % boundaries
    for j=1:Nj
        % left-most edge of mesh
        [upper(j,1), last_in(j,1)] = filter(B, A, lowerOld(j,1), last_in(j,1));  
        % right-most edge of mesh
        [lower(j,Ns), last_in(j,2)] = filter(B, A, upperOld(j,Ns), last_in(j,2)); 
        % highest edge of mesh
        [right(1,j), last_in(j,3)] = filter(B, A, leftOld(1,j), last_in(j,3));
        % lowest edge of mesh
        [left(Ns,j), last_in(j,4)] = filter(B, A, rightOld(Ns,j), last_in(j,4));
    end
    
    
    % output 
    y(n) = vj(outX, outY);
    %y(n) = vj(round(Nj/2), round(Nj/2));
    %y(n) = upper(outX, outY);
    
    if plotOn == 1
        if ( mod(n,1) == 0 ),
            vjb(2:end-1,2:end-1) = vj(1:end, :);
            surf(vjb);
            colormap Winter;
            %mesh(vj);
            axis([0, Nj+2, 0, Nj+2, -0.15 0.15])
            zlabel('Velocity');
            pause(0.1);
            %pause;
        end
    end
    
    %keyboard
    
end





