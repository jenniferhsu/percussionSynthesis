% create a 2048-sample long sine table
N = 2048;
n = 0:1/N:(1-(1/N));
x = sin(2*pi*1*n);

% write to file
fileID = fopen('sinWavetable.txt', 'w');
for i=1:N
    fprintf(fileID, '%.12f ', x(i));
end
fclose(fileID);
