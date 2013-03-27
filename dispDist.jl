function [RsqPMF, RsqCount] = dispDist(T,binPos,EQUAL, VERBOSE, mode)
# [RsqPMF, RsqCount] = dispDist(T,binPos, [EQUAL], [VERBOSE], [mode])
# T       - Trajectory array, it could be in any dimension. T{traj}(frame,xyz)
# binPos  - Bin centers for the histograms of the logarithmic displacements. 
#           Must be a row vector of logarithmically spaced bin centers.
# EQUAL   - True if all traj have same length, else false. Default false
# VERBOSE - Verbose output. Default true 
# mode    - 'shift' or 'cut'. Default 'shift'
#
# RsqCount ~ RsqCountPMFance squared for every dt, and for every trajectory
#            RsqCount = {traj}[bin,dt], if EQUAL == true then it returns
#            RsqCount = [bin,dt] average of all traj. It starts at dt=1 
# RsqPMF   ~ RsqCount RsqCountPMFribution for every dt, normalized so that 
#            the sum(RsqPMF(:,i)) = 1. Every trajectory has the same weight, 
#            regardless of length. If EQUAL == true, same as RsqCount
#
# The mean of that is the MSD. (tested)
#
# gP 7/15/2011

if nargin < 3 || isempty(EQUAL)
    EQUAL = false;
end

if nargin < 4 || isempty(VERBOSE)
    VERBOSE = true;
end

if nargin < 5 || isempty(mode)
    mode = 'shift';
end

Nt = length(T);                             # Number of trajectories
Nbins = length(binPos);
RsqCount = cell(1,Nt);
RsqPMF = cell(1,Nt);
                                            # Convert from log bin center to edges
halfbin = (log10(binPos(2))-log10(binPos(1)))/2;
binEdge = [log10(binPos)-halfbin, log10(binPos(end))+halfbin];
binEdge = 10.^binEdge;


if strcmp(mode,'shift')
    
    if ~EQUAL
        
        if VERBOSE
            progressText(0,'Calculating dispDist')
        end
        
        for i = 1:Nt                        # Go through trajectories
            Np = size(T{i},1);              # Number of points in trajectory
            
            RsqCount{i} = zeros(Nbins, Np-1);
            RsqPMF{i} = zeros(Nbins, Np-1);
            
            for dt = 1:(Np-1)               # Time interval
                
                lag = 1:(Np-dt);
                
                Rsq = sum((T{i}(lag+dt,:) - T{i}(lag,:)).^2, 2);
                
                count = histc(Rsq,binEdge,1); count(end)=[];    # Remove last bin. It only contains the value for the edge of the bin
                RsqCount{i}(:,dt) = count;
                
                RsqSum = sum(RsqCount{i}(:,dt));
                if RsqSum == 0              # Prevents division by zero
                    RsqPMF{i}(:,dt) = zeros(Nbins, 1);
                else
                    RsqPMF{i}(:,dt) = count./RsqSum;
                end
            end
            
            if VERBOSE
                progressText(i/Nt)
            end
        end
        
    else                                    # For equal length 
        Np = size(T{1},1);
        RsqCount = zeros([length(binPos), Np-1]);
        RsqPMF = zeros([length(binPos), Np-1]);
        
        if VERBOSE
            progressText(0,'Calculating dispDist')
        end
        
        for dt = 1:(Np-1)                   # Go over time interval
            
            for i=1:Nt                      # Go over the traj 
                
                lag = 1:(Np-dt);
                                            
                Rsq = sum((T{i}(lag+dt,:) - T{i}(lag,:)).^2, 2);

#                 count = hist(Rsq,binPos)';                
                count = histc(Rsq,binEdge,1); count(end)=[];   # Remove last bin. It only contains the value for the edge of the bin
                RsqCount(:,dt) = RsqCount(:,dt) + count;
            end
            
            RsqSum = sum(RsqCount(:,dt));
            if RsqSum == 0
                RsqPMF(:,dt) = zeros(Nbins,1);
            else
                RsqPMF(:,dt) = RsqCount(:,dt)./ sum(RsqCount(:,dt));
            end
            
            if VERBOSE
                progressText(dt/(Np-1))
            end
        end
        
    end