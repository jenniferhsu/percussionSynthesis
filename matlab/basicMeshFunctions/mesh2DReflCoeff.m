function y = mesh2DReflCoeff(ex, Nj, N, fs, inX, inY, outX, outY, refl, plotOn)
%y = mesh2DReflCoeff(ex, Nj, N, fs, inX, inY, outX, outY, refl, plotOn)
% this is a rectangular 2d mesh waveguide implementation that allows you
% to set the input excitation
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
%   refl: reflection coefficient (-1 for a perfect reflection,
%        -0.99 works well)
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
% ex = sin(2*pi*2004.6*(0:1/fs:1));
% Nj = 11;
% N = 44100;
% fs = 44100;
% inX = 5;
% inY = 5;
% outX = 5;
% outY = 5;
% refl = -0.95;
% plotOn = 1;
% y = mesh2DReflCoeff(ex, Nj, N, fs, inX, inY, outX, outY, refl, plotOn);
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

% initialize delay lines with some value?
%upper(round(Nj/2),round(Nj/2)) = 1.0;
%lower(Nj,Ns) = 0.0;
%left(Ns,Nj) = 0.0;
%right(Nj,Nj) = 0.0;

% velocity at the junctions
vj = zeros(Nj, Nj);

% for plotting
X = (1:Nj)';
XX = repmat(X, [1, Nj]);
Y = 1:Nj;
YY = repmat(Y, [Nj, 1]);

% for each output sample
for n=1:N
    
    % get the input excitation value
    if length(ex) < n
        % do nothing
    else
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
    
    %keyboard

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
    % left-most edge of mesh
    upper(:,1) = refl * lowerOld(:,1);
    % right-most edge of mesh
    lower(:,Ns) = refl * upperOld(:,Ns);
    % highest edge of mesh
    right(1,:) = refl * leftOld(1,:);
    % lowest edge of mesh
    left(Ns,:) = refl * rightOld(Ns,:);
    
    
    
    % output 
    %y(n) = vj(round(Nj/2), round(Nj/2));
    %y(n) = upper(outX, outY);
    %y(n) = sum(sum(vj)); % elliot said this is bassier
    
    y(n) = vj(outX, outY);
    
    %keyboard
    
    if plotOn == 1
        if ( mod(n,100) == 0 )
            %vj
            surf(XX, YY, vj);
            axis([0, Nj, 0, Nj, -1 1])
            zlabel('Velocity');
            pause(0.01);
            %pause
        end
    end
    
   
    
end





