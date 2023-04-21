%Matlab Code used for running the simulations reported in 
%Benjamin, Griffin, and Douglas, "A nonparametric technique for analysis of state-trace functions:
%with an application to recognition memory"

%prepared by Michael Griffin
%last updated 6.25.2018

%This program generates the simulated data that will eventually be tested
%by the isoperm program. It is called by runsims.m, but can also be run on
%its own. It creates and saves the 'data' variable to 'simulated.mat'. 

%Variables that were updated as we ran our simulations: 
%concaveup: determines shape of the underlying function 
%spreadsize: Determines how separated the two conditions are, with higher spreadsize leading to a smaller overlap region.
%nhalf: number of points per condition.

%noise and interactsize are also used, but normally are received from runsims.m.



concaveup = 0;          %0 = linear, 1 = concave up, -1 concave down
nsubs = 100;
spreadsize = .05;       %.25 or .05
nhalf = 4;              %50, 10, or 4
npoints = nhalf*2;

Strbase = .15; %Added to all Str values. 
duradj = .7/nhalf; 

spread = (0:nhalf) * duradj + duradj*.5;
[~,index] = min(abs(spread - spreadsize)); %Find smallest deviation from desired size where points still fall in between each other.
spread = spread(index);

if ~exist('runsim', 'var')
    noise = 0;          %.1, .2, .4
    interactsize = 0;   % 0, .1, .15, .2
end
interactsigns = Shuffle([-ones(nsubs/2,1); ones(nsubs/2,1)]);

data = zeros(nsubs, npoints*2);

for k = 1:nsubs
    averages = zeros(npoints,2);
    for yaxis = 0:1 
        for cond = 1:2
            if cond == 1
                cspread = spread/2;
                adj = 0;
            else
                cspread = -spread/2;
                adj = nhalf;
            end
            
            if cond == 1 && yaxis
                cinteract = interactsize*interactsigns(k);
            else
                cinteract = 0;
            end
            
            for n = 1:nhalf
                cdur = n*duradj;                
                Str = Strbase + cdur + cspread;
                StrRad = Str*(pi/2) + (3*pi)/2; %Get concave up section of circle.
                switch concaveup
                    case 0
                        observed = Str;
                    case 1
                        if yaxis == 0
                            observed = cos(StrRad);
                        else
                            observed = sin(StrRad)+1; %+1 so range is 0 to 1 rather than -1 to 0.
                        end
                    case -1   %Flip x and y to get concave down.
                        if yaxis == 0
                            observed = sin(StrRad)+1;
                        else
                            observed = cos(StrRad);
                        end
                end
                
                observed = observed + cinteract;
                observed = observed+randn*noise;
                averages(n+adj,1+yaxis) = observed;
            end
        end
    end
    data(k,1:npoints) = averages(:,1);
    data(k,npoints+1:end) = averages(:,2);
end

save('simulated.mat', 'data');