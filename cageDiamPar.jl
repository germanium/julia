 
# CD = cageDiamPar(file, window)
# file    - File containing the trajectories. First column is time, the rest
#	    are x1,y1,z1...xN,yN,zN where N is the number of trajectories.
# window  - Sliding time window length in frames. Must be odd. Default is 7
#
# CD      ~ Returns caging diameter cell array for all the trajectories between 
#           every point and a radius of +/-window/2 time. For the points below 
#           window/2 and above (end - window/2) only a partial window is 
#           considered, hence they are not to be trusted. 
#           It's returned in the units of the trajectory.
#
# This is 1.5x times slower than the matlab code that uses the c function
# gP 3/26/2013

using Distance

function cageDiamPar(file::String, window::Int)

    T = readdlm(file)
    T = T[:, 2:end]       		# Remove first column containing time

    Nt = convert(Int64,size(T,2)/3)
    Np = size(T,1)
    CD = zeros(Np,Nt)

    halfWin = floor(window/2)

    iT2 = 1
    @parallel for iT=1:3:size(T,2)      # Go over x-dimension of traj
	
	TT = T[:, iT:(iT+2)]		# Single trajectory
	
	for t=1:Np                      # Go over time
            
            if t-halfWin < 1
		winPts = 1:(t+halfWin)
            elseif t+halfWin > Np
		winPts = (t-halfWin):Np
            else
		winPts = (t-halfWin):(t+halfWin)
            end
            
            Twin = TT[winPts,:]         # Points in window
            
            R = pairwise(Euclidean(), Twin', Twin')
            				
            CD[t, iT2] = max(R)       	# R will have duplicate displacements
	    				# but it doesn't matter in this case
	end				
	iT2 = iT2+1

    end

    return CD

end