

using Distributed

addprocs(3)

@everywhere using SharedArrays

S = SharedArray{Int,2}((3,4), init = 0.0f0)

@everywhere begin
	function myrange(q::SharedArray)
		idx = indexpids(q)
		if idx == 0 # This worker is not assigned a piece
			return 1:0, 1:0
		end
		nchunks = length(procs(q))
		splits = [round(Int, s) for s in range(0, stop=size(q,2), length=nchunks+1)]
1:size(q,1), splits[idx]+1:splits[idx+1]
	end
end
              
@everywhere begin
	function advection_chunk!(q, u, irange, jrange, trange)
		@show (irange, jrange, trange)  # display so we can see what's happening
		for t in trange, j in jrange, i in irange
			q[i,j,t+1] = q[i,j,t] + u[i,j,t]
		end
	q
	end
end
              
@everywhere advection_shared_chunk!(q, u)=advection_chunk!(q, u, myrange(q)..., 1:size(q,3)-1)

advection_serial!(q, u)=advection_chunk!(q, u, 1:size(q,1), 1:size(q,2), 1:size(q,3)-1);

function advection_parallel!(q, u)
	for t = 1:size(q,3)-1
		@sync @distributed for j = 1:size(q,2)
			for i = 1:size(q,1)
				q[i,j,t+1]= q[i,j,t] + u[i,j,t]
			end
		end
	end
	q
end


function advection_shared!(q, u)
	@sync begin
		for p in procs(q)
			@async remotecall_wait(advection_shared_chunk!, p, q, u)
		end
	end
	q
end


q = SharedArray{Float64,3}((500,500,500));

u = SharedArray{Float64,3}((500,500,500));

@time advection_serial!(q, u);

@time advection_parallel!(q, u);

@time advection_shared!(q, u);


######################################################################
#
#
#			Results
#
#
#
#	julia> @time advection_serial!(q, u);
#	(irange, jrange, trange) = (1:500, 1:500, 1:499)
#	  3.743086 seconds (362.55 k allocations: 18.751 MiB, 5.60% gc time)
#
#	julia> @time advection_parallel!(q, u);
#	  4.472376 seconds (1.13 M allocations: 57.873 MiB, 12.16% gc time)
#
#	julia> @time advection_shared!(q, u);
#	      From worker 3:    (irange, jrange, trange) = (1:500, 168:333, 1:499)
#	      From worker 4:    (irange, jrange, trange) = (1:500, 334:500, 1:499)
#	      From worker 2:    (irange, jrange, trange) = (1:500, 1:167, 1:499)
#	  1.985823 seconds (180.81 k allocations: 9.370 MiB)
###




