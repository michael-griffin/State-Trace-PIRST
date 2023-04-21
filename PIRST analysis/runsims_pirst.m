%Matlab Code used for running the simulations reported in 
%Benjamin, Griffin, and Douglas, "A nonparametric technique for analysis of state-trace functions:
%with an application to recognition memory"

%prepared by Michael Griffin
%last updated 6.25.2018

%This is the code used for running a series of simulations at once. After entering a set of noise
%values and latent interaction strengths, it will generate data (or load
%from already saved data), and then run the isotonic permutation test on it for n
%iterations. It then saves summary statistics such as the probability that swaps 
%led to a greater RSS, the mean RSS pre and post swap, and the percentage of points 
%in the overlap region.

%Two parameters that are set here are passed to simdata.m: a noise
%value (from the array noisevals) and an interaction value (from the array
%interactvals). These correspond to the gaussian noise that gets added to 
%each data point, and the latent interaction value that gets added (or 
%subtracted) as a constant to the y-axis. Both are taken from arrays that are 
%looped through, allowing runsims to generate all the results of each 
%combination of noise + interaction in one go.


starttime = tic;
savefile = 0; %save datasets into 'allsim.mat';
writefile = 0; %write results to excel file 'simulations_iso.xlsx'

runsim = 1; %If 0, load data from 'allsim.mat', if 1 runs simdata;
dataonly = 0; %use with savedata to generate and then save data without running the (time consuming) permutation analysis.

noisevals = [.1, .2, .4];
interactvals = [0, .1, .15, .2];
nrows = length(noisevals)*length(interactvals);

if runsim
    iterations = 10; 
    fulldata = cell(nrows,iterations);
else
    load('allsim.mat');
    iterations = size(fulldata,2);
end


isosummary = zeros(nrows,6); 


%Form of summary: rows are simulations are under different noise and interaction values
%Columns give different summary statistics.

%Row 1-4: No, weak, moderate, and strong interactions under low noise.
%Row 5-8: Ditto, moderate noise
%Row 9-12: Ditto, high noise
%Columns are:
%1,2: Base RSS, Swapped RSS
%3-4: P Swap led to greater, P swap led to less,
%5-6: P 'swap' didn't actually change order, percentage of points in overlap region.

for z = 1:length(noisevals)
    noise = noisevals(z); %.1, .2, .4
    
    for v = 1:length(interactvals)
        interactsize = interactvals(v);
        disp(['noise: ', num2str(noise), char(13), 'interact: ', num2str(interactsize)]);
        
        isos = zeros(iterations, 6); %1RSS no swap 2 Swap, 3 pgreater 4 pless/tie 5 psameorder 6 overlapavg
        
        row = v+(z-1)*length(interactvals); 
        for w = 1:iterations
            if runsim
                simdata;    
                fulldata{row,w} = data;
            else
                data = fulldata{row,w};
            end
            
            if dataonly
                continue;
            end
            
            %run the analysis
            isoperm;
            
            isos(w,1) = rssmeanbase;
            isos(w,2) = rssmeanswap;
            isos(w,3) = pgreater;
            isos(w,4) = 1-pgreater-psameorder; %This includes no change in rss as evidence not in favor of 2 process.
            isos(w,5) = psameorder;
            isos(w,6) = overlapavg;
        end
        
        if dataonly
            continue;
        end        
        
        for n = 1:6
            isosummary(row,n) = mean(isos(:,n)); 
        end
    end
end
if savefile
    save('allsim.mat', 'fulldata', 'spread', 'nhalf', 'concaveup');
end
if writefile
   xlswrite('simulations_iso.xlsx', isosummary, 1); 
end
elapsed = toc(starttime)/60