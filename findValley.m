%%=========================================================================
% Synopsis     : Finds valley between two mountains in a data stream
% Last updated : 2019-04-17
%%=========================================================================
clc;clear;close all;

%% generation of random data with two prominant peaks
hightMountain   = randi([10 20],1,1);
mountains       = [1:hightMountain hightMountain-1:-1:hightMountain/2 hightMountain/2+1:hightMountain hightMountain-1:-1:1];
inputData       = [randn(1,randi([30 100],1,1)) mountains randn(1,randi([30 100],1,1))];
inputData       = awgn(inputData,-2);

%% Smoothing the data
movAvgFiltOrder = 16;
smoothedData    = filter(ones(1,movAvgFiltOrder)/movAvgFiltOrder,1,inputData);

%% finding peaks from the smoothed curve
[vals,pos]      = findpeaks(smoothedData,"DoubleSided");
idx             = find(vals > mean([min(vals),max(vals)])); % removing peaks which are below the average,
                                                            % this is helpfull if some lower amplitude
                                                            % peaks are caught at sides of the prominent peaks.

peak_pos1       = pos(idx(1));    % farthest prominent peak on one side of valley
peak_pos2       = pos(idx(end));  % farthest prominent peak on other side of valley

%% correcting the peak positions for filter delay (Coarse tuning)
peak_pos1       = peak_pos1 -  movAvgFiltOrder/2;
peak_pos2       = peak_pos2 -  movAvgFiltOrder/2;

%% Looking for any other peaks nearby (Fine Tuning)
tuning_range    = min(movAvgFiltOrder,floor(abs((peak_pos2-peak_pos1-1)/2)));
p1              = peak_pos1-tuning_range:peak_pos1+tuning_range;
p2              = peak_pos2-tuning_range:peak_pos2+tuning_range;
[~,r1]          = max(inputData(p1));
[~,r2]          = max(inputData(p2));
peak_pos1       = p1(r1);
peak_pos2       = p2(r2);

%% Identifying valley
valley          = inputData(peak_pos1:peak_pos2);

%% findind minimum point in valley
[v,p]           = min(valley);
valley_amp      = v; % amplitude of minimum point
valley_pos      = peak_pos1+p-1; % finding the valley point in the entire array

%% Prints
fprintf('Valley point is %d\n',valley_pos);

%% Plots
figure('units','normalized','outerposition',[0 0 1 1])
plot(inputData,'-k');hold on;pause(1.5);
plot(peak_pos1:peak_pos1+length(valley)-1,valley,'--b','Linewidth',2);pause(1);
plot(valley_pos,valley_amp,'ro','Linewidth',1.5);
legend('Data','Valley','Valley Point');
