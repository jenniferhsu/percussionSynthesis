x = 0:0.01:1;

% initial point
P0x = 0;
P0y = 1;

% end point
P2x = 1;
P2y = 0;

P1xVec = [zeros(1, 11) 0.1:0.1:1];
P1yVec = [0:0.1:1 ones(1, 10)];

% convex when P1x is 0
% concave when P1y is 1
% linear when P1x = 0 and P1y = 1

N = length(P1xVec);

figure
for i=1:N

    P1x = P1xVec(i)
    P1y = P1yVec(i)

    % Bezier curves
    Bx = (1-x) .* ((1-x)*P0x + x*P1x) + x .* ((1-x)*P1x + x*P2x);
    By = (1-x) .* ((1-x)*P0y + x*P1y) + x .* ((1-x)*P1y + x*P2y);

    plot(Bx, By)
    hold on
    pause
end
