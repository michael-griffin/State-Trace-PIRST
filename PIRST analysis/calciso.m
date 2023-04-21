function isofit = calciso(dat)
%Matlab Code used for running the simulations reported in 
%Benjamin, Griffin, and Douglas, "A nonparametric technique for analysis of state-trace functions:
%with an application to recognition memory"

%prepared by Michael Griffin
%last updated 6.25.2018


%This function fits an isotonic regression curve to data 'dat', via back
%averaging.


dat = sortrows(dat, 1);
dat = [dat dat(:,2)]; %last column will eventually become isofit

%Average across ties first.
uniq = unique(dat(:,1));
for n = 1:length(uniq)
    if length(find(dat(:,1) == uniq(n))) > 1 %if there's a tie.
        indexes = find(dat(:,1) == uniq(n));
        dat(indexes,3) = mean(dat(indexes,3));
    end
end

for n = 2:length(dat)
    if dat(n,3) < dat(n-1,3)
        done = 0; 
        index = n-1;
        newmean = mean(dat(index:n,3));
        
        %Back averages. If the new average is still smaller than the previous
        %point, that point is included, and then the average is recalculated.
        %Continues until the result is monotonic.
        while ~done
            problem = find(dat(1:index-1,3) > newmean, 1, 'last');
            if isempty(problem)
                done = 1;
            else
                index = problem;
                newmean = mean(dat(index:n,3));
            end
        end
        dat(index:n,3) = newmean;
    end
end
isofit = dat(:,3);
end

