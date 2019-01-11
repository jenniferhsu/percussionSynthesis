function y = mesh2DRectReflCoeff(ex, Njr, Njc, N, fs, inX, inY, outX, outY, refl, plotOn)
%y = mesh2DRectReflCoeff(ex, Njr, Njc, N, fs, inX, inY, outX, outY, refl, plotOn)
% this is a rectangular 2d mesh waveguide implementation that allows you
% to set the input excitation. THIS ONE IS RECTANGULAR NOT SQUARE
% inputs:
%   ex: the input excitation as an array
%   Njr: the number of row junctions
%   Njc: the number of column junctions
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
% to call this function for a 4x3 square mesh:
%
% ex = 1;
% Njr = 4;
% Njc = 3;
% N = 44100;
% fs = 44100;
% inX = 4;
% inY = 3;
% outX = 4;
% outY = 3;
% refl = -0.95;
% plotOn = 0;
% y = mesh2DRectReflCoeff(ex, Njr, Njc, N, fs, inX, inY, outX, outY, refl, plotOn)
%

y = zeros(1,N);
Nsr = Njr + 1;    % number of sections 
Nsc = Njc + 1;

% delay line structures
upper = zeros(Njr, Nsc);
lower = zeros(Njr, Nsc);
left = zeros(Nsr, Njc);
right = zeros(Nsr, Njc);

% scattering junction port structures 
portInsUpper = zeros(Njr, Njc);
portInsLower = zeros(Njr, Njc);
portInsLeft = zeros(Njr, Njc);
portInsRight = zeros(Njr, Njc);
portOutsUpper = zeros(Njr, Njc);
portOutsLower = zeros(Njr, Njc);
portOutsLeft = zeros(Njr, Njc);
portOutsRight = zeros(Njr, Njc);

% initialize delay lines with some value?
%upper(round(Nj/2),round(Nj/2)) = 1.0;
%lower(Nj,Ns) = 0.0;
%left(Ns,Nj) = 0.0;
%right(Nj,Nj) = 0.0;

% velocity at the junctions
vj = zeros(Njr, Njc);

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
    portInsUpper = upperOld(:,1:Njc);
    portInsLower = lowerOld(:,2:Nsc);
    portInsRight = rightOld(1:Njr,:);
    portInsLeft = leftOld(2:Nsr,:);

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
    upper(:,2:Nsc) = portOutsUpper;
    lower(:,1:Njc) = portOutsLower;
    right(2:Nsr,:) = portOutsRight;
    left(1:Njr,:) = portOutsLeft;
    
    % boundaries
    % left-most edge of mesh
    upper(:,1) = refl * lowerOld(:,1);
    % right-most edge of mesh
    lower(:,Nsc) = refl * upperOld(:,Nsc);
    % highest edge of mesh
    right(1,:) = refl * leftOld(1,:);
    % lowest edge of mesh
    left(Nsr,:) = refl * rightOld(Nsr,:);
    
    
    
    % output 
    %y(n) = vj(round(Nj/2), round(Nj/2));
    %y(n) = upper(outX, outY);
    %y(n) = sum(sum(vj)); % elliot said this is bassier
    
    y(n) = vj(outX, outY);
    
    %keyboard
    
    if plotOn == 1
        if ( mod(n,100) == 0 ),
            %vj
            surf(vj);
            axis([0, Nj, 0, Nj, -1 1]) % maybe change this
            zlabel('Velocity');
            pause(0.01);
            %pause;
        end
    end
    
   
    
end





