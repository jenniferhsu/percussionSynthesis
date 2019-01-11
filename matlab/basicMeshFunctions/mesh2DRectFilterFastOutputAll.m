function y = mesh2DRectFilterFastOutputAll(ex, Njr, Njc, N, fs, inX, inY, B, A, plotOn)
%Y = mesh2DRectFilterFastOutputAll(ex, Njr, Njc, N, fs, inX, inY, B, A, plotOn)
% this is a rectangular 2d mesh waveguide implementation that allows you
% to set the input excitation, specify the input point, and also
% change the filters at the boundaries. This is a fast version that only
% allows you to use first order filters. The resulting signal is the sum
% of all poins on the mesh (compare with mesh2DRectFilterFast() which 
% outputs from a single point
%
% inputs:
%   ex: the input excitation as an array
%   Njr: the number of row junctions
%   Njc: the number of column junctions
%   N: the output signal length
%   fs: the sample rate
%   inX: input point x coordinate 
%   inY: input point y coordinate 
%   B: filter coefficients (that multiply with the input)
%   A: filter coefficients (that multiply with the output)
%   plotOn: set to 1 to plot as an animation and 0 otherwise
% output:
%   y: the generated signal
%
% note: for inX, inY, the input signal is divided by 4 and
% input to the junction at (inX, inY) from all 4 directions. the output is
% the sum of velocities at every junction on the mesh at each timestep
%
% to call this function for a 11x11 square mesh:
%
% ex = 1;
% Njr = 31;
% Njc = 31;
% N = 44100;
% fs = 44100;
% inX = 31;
% inY = 31;
% B = -[0.9 0.1];
% A = 1;
% plotOn = 0;
% y = mesh2DRectFilterFastOutputAll(ex, Njr, Njc, N, fs, inX, inY, B, A, plotOn);
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

% velocity at the junctions
vj = zeros(Njr, Njc);

% for plotting, let's extend vj to have the boundaries
%vjb = zeros(Nj+2, Nj+2);
vjb = zeros(Njr+2, Njc+2);
            
           
% filters for the reflections - we'll use the simplest low pass filter
%last_in = zeros(Nj,4); % 4 for each side that needs a reflection
%%
%THE LAST IN STUFF NEEDS TO BE MODIFIED
%SOMETHING LIKE:
last_in_upper = zeros(Njr, 1);
last_in_lower = zeros(Njr, 1);
last_in_left = zeros(Njc, 1);
last_in_right = zeros(Njc, 1);
%%


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
    portInsUpper = upperOld(:,1:Njc);
    portInsLower = lowerOld(:,2:Nsc);
    portInsRight = rightOld(1:Njr,:);
    portInsLeft = leftOld(2:Nsr,:);

    % calculate the velocity(?)
    vj = 0.5 * (portInsUpper + portInsLower + ...
        portInsLeft + portInsRight);

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
    [upper(1:Njr,1), last_in_upper(1:Njr)] = simple_filter(B, A, lowerOld(1:Njr,1), last_in_upper(1:Njr)); 
    % right-most edge of mesh 
    [lower(1:Njr,Nsc), last_in_lower(1:Njr)] = simple_filter(B, A, upperOld(1:Njr,Nsc), last_in_lower(1:Njr)); 
    % highest edge of mesh
    [right(1,1:Njc), last_in_right(1:Njc)] = simple_filter(B, A, leftOld(1,1:Njc), last_in_right(1:Njc)');
    % lowest edge of mesh
    [left(Nsr,1:Njc), last_in_left(1:Njc)] = simple_filter(B, A, rightOld(Nsr,1:Njc)', last_in_left(1:Njc));
   

    % output 
    %y(n) = vj(outX, outY);
    %y(n) = vj(round(Nj/2), round(Nj/2));
    %y(n) = upper(outX, outY);
    y(n) = sum(sum(vj));
    
    if plotOn == 1
        if ( mod(n,1) == 0 )
            vjb(2:end-1,2:end-1) = vj(1:end, :);
            surf(vjb);
            colormap Winter;
            %mesh(vj);
            axis([0, Njr+2, 0, Njc+2, -0.15 0.15])
            zlabel('Velocity');
            pause(0.1);
            %pause;
        end
    end
    
    %keyboard
    
end





