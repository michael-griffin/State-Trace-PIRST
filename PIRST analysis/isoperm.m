%Matlab Code used for running the simulations reported in 
%Benjamin, Griffin, and Douglas, "A nonparametric technique for analysis of state-trace functions:
%with an application to recognition memory"

%prepared by Michael Griffin
%last updated 6.25.2018

%This code runs the isotonic permutation test. It fits isotonic regression
%lines by calling calciso.m. It then determines the overlap region (points
%in the data set whose conditions labels are eligible to swapped), and then
%runs the permutation test for that subject n times: it swaps condition labels, 
%refits the isotonic regression lines to the new 'conditions', and records
%the new fit and whether it is better or worse than pre-swap.

%The permutations variable determines the number of swaps + refits done for
%each subject. The simulated data the test runs on will usually be generated/loaded 
%already by runsims.m, but the program can also be run on its own.





permutations = 100; %how many permutations are done for swapping condition labels.

if exist('runsim', 'var') && runsim == 0 %Data will already be loaded by runsims
else
    load('simulated.mat'); %contains the 'data' variable
end

allsubs = [1:size(data,1)]';
nsubs = size(data,1);

nhalf = size(data,2)/4;
npoints = nhalf*2; 
tol = 1e-9; %needed in pgreater down the line. Annoying.

summaryiso = zeros(nsubs, 6);
%Columns are:
%1-2: Base RSS, Swapped RSS
%3-4: P Swap led to greater, P swap led to less,
%5-6: P 'swap' didn't actually change order, percentage of points in overlap region.


for k = 1:nsubs
    crow = data(k,:);
    averages = [crow(1:npoints)' crow(npoints+1:end)' (1:npoints)' [ones(1,nhalf) 2*ones(1,nhalf)]'];
    sorted = sortrows(averages, 1);
    
    
    iso = zeros(npoints, 12);
    iso(:,1:size(sorted,2)) = sorted;
    %Col 1-2: Base data
    %Col 3-4: Condition 1 vs. Condition 2
    %Col 5: Overlap
    %Col 6: post swap condition
    %Col 7-9: Pre swap iso fits. (7: across cond, 8: Cond 1, 9: Cond 2)
    %Col 10-11: Post swap iso fits
    
    %Find overlap region. Points in this region will later be swapped for
    %the permutation test.
    %Find start:
    sorty = sortrows(averages,2);
    for n = 1:npoints-1
        index = iso(n,3);
        index = find(sorty(:,3) == index);
        %if the next point that has a higher x or y value is in the
        %other attention condition, or if there are any y values from
        %the other condition that are smaller than the current.
        if iso(n+1,4) ~= iso(n,4) || sum(sorty(1:index-1,4) ~= iso(n,4))
            break;
        end
        if index ~= npoints && sorty(index+1,4) ~= iso(n,4)
            break;
        end
    end
    start = n;
    
    
    indexone = find(iso(:,4) == 1);
    indextwo = find(iso(:,4) == 2);
    
    isoone = iso(indexone,:);
    isotwo = iso(indextwo,:);
    
    
    %Determines which are inside overlap region
    tocut = [];
    for n = 1:npoints
        switch iso(n,4)
            case 1
                if iso(n,1) > max(isotwo(:,1)) && iso(n,2) >= max(isotwo(:,2)) ...
                        || iso(n,1) < min(isotwo(:,1)) && iso(n,2) <= min(isotwo(:,2))
                    tocut = [tocut, n];
                end
            case 2
                if iso(n,1) > max(isoone(:,1)) && iso(n,2) >= max(isoone(:,2)) ...
                        || iso(n,1) < min(isoone(:,1)) && iso(n,2) <= min(isoone(:,2))
                    tocut = [tocut, n];
                end
        end
    end
    iso(:,5) = 1;
    iso(tocut,5) = 0;
    overlap = find(iso(:,5));
    
    %calciso is the program that does most of the work. Line below doesn't
    %actually get used for much, more a check to see if the fit's working.
    iso(:,7) = calciso(iso(:,1:2));
    %Base fits, to be compared to the permutations.
    iso(indexone,8) = calciso(iso(indexone,1:2));
    iso(indextwo,9) = calciso(iso(indextwo,1:2));
    
    
    %Calculate base RSS for both attention conditions
    rssbaseone = 0;
    rssbasetwo = 0;
    for n = 1:length(indexone)
        index = indexone(n);
        rssbaseone = rssbaseone + (iso(index,2)-iso(index,8))^2;
    end
    for n = 1:length(indextwo)
        index = indextwo(n);
        rssbasetwo = rssbasetwo + (iso(index,2)-iso(index,9))^2;
    end
    
    
    permtest = 0;
    if sum(iso(:,5)) > 0 && sum(iso(overlap,4) == 2) && sum(iso(overlap,4) == 1)%If there's any overlap, do permutation test.
        permtest = 1;
    end
    
    
    permdata = NaN(permutations,2);    
    
    %Swap conditions
    if permtest
        for p = 1:permutations
            swaps = overlap(randperm(length(overlap)));
            iso(:,6) = iso(:,4);
            iso(overlap,6) = iso(swaps,4);

            
            indexnewone = find(iso(:,6) == 1);
            indexnewtwo = find(iso(:,6) == 2);
            
            iso(:,10) = 0; %clears array for readability.
            iso(:,11) = 0;
            iso(indexnewone,10) = calciso(iso(indexnewone,1:2));
            iso(indexnewtwo,11) = calciso(iso(indexnewtwo,1:2));
            
            rss = 0;
            for n = 1:length(indexnewone)
                index = indexnewone(n);
                rss = rss + (iso(index,2)-iso(index,10))^2;
            end
            rssfull = rss;
            
            rss = 0;
            for n = 1:length(indexnewtwo)
                index = indexnewtwo(n);
                rss = rss + (iso(index,2)-iso(index,11))^2;
            end
            rssdiv = rss;
            
            permdata(p,1) = rssbaseone+rssbasetwo;
            permdata(p,2) = rssfull+rssdiv;
        end
        summaryiso(k,1) = rssbaseone+rssbasetwo; %RSS pre swap
        summaryiso(k,2) = mean(permdata(:,2)); %average swapped RSS
        summaryiso(k,3) = sum((permdata(:,2) - permdata(:,1)) > tol)/length(permdata);
        summaryiso(k,4) = sum((permdata(:,2) - permdata(:,1)) < -tol)/length(permdata);
        if length(overlap) < 50
            summaryiso(k,5) = 1/nchoosek(length(overlap), sum(iso(overlap,4) == 1)); 
        else
            summaryiso(k,5) = 0; %ridiculously small fraction. If nchoosek is used, gives warning message.
        end
    else
        summaryiso(k,1) = rssbaseone+rssbasetwo; %RSS pre swap
        summaryiso(k,2) = NaN;
        summaryiso(k,3) = NaN;
        summaryiso(k,4) = NaN;
        summaryiso(k,5) = NaN;
    end
    summaryiso(k,6) = length(overlap)/npoints;
end

relevant = 1:length(summaryiso);
toremove = [];
for n = 1:nsubs
    if sum(isnan(summaryiso(n,:)));
        toremove = [toremove n];
    end
end
relevant(toremove) = [];

%Measures recorded by runsims
rssmeanbase = mean(summaryiso(relevant,1));
rssmeanswap = mean(summaryiso(relevant,2));
pgreater = mean(summaryiso(relevant,3)); %Each subject has a p greater, this averages.
pless = mean(summaryiso(relevant,4));
psameorder = mean(summaryiso(relevant,5));
overlapavg = mean(summaryiso(relevant,6));

