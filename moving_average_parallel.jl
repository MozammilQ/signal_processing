
using Distributed
addprocs(3)

@everywhere procs() using SharedArrays

function mov_avg_db_par(arr::SharedArray{Float32,1},M::Int)
	
	@everywhere workers() begin

		M=$M

		function mov_avg_worker(arr::SharedArray{Float32,1},M::Int,start_::Int,end_::Int)
			@info "worker " * string(myid()) * " at work"
			for i in 1:1:end_-start_-M+1
				temp_var=0.0f0
		                for j in i:1:i+M-1
		                        @inbounds temp_var+=arr[start_+j-1]
		                end
				temp_var/=M
				temp_var+=1.0f0
				temp_var=10.0f0*log10(temp_var)
				@inbounds arr[i]=temp_var
		        end
		end
	
	end


	@info "Main node mov_avg_db_par() invoked"
	chunk=round(Int,length(arr)/nworkers())
	
	wrkrs=workers()

	@sync for i in wrkrs
	
		start_= round(Int, chunk * ( i - wrkrs[1] )) + 1
		end_=start_+chunk+M-2
	
	
		if i==wrkrs[end]
			@fetchfrom i mov_avg_worker(arr,M,start_,length(arr))
		else
			@async @fetchfrom i mov_avg_worker(arr,M,start_,end_)
		end
	
	end
end


##########################################################################################################



function abs2_par!(abs2_arr::SharedArray{Float32,1},fft_reduced_arr::SharedArray{Complex{Float32},1})

	@everyone workers() begin
		function abs2_worker(abs2_arr::SharedArray{Float32,1},fft_reduced_arr::SharedArray{Complex{Float32},1},start_::Int,end_::Int)
		for i in start_:1:end_
			abs2_arr[i]=abs2(fft_arr[i])
		end
	end

	chunk=round(Int,length(abs2_arr)/nworkers())
end



function get_psd(fft_arr)
	abs2_arr=SharedArray{Float32,1}(round(Int,length(fft_arr)/decim);pids=workers())
	return psd_to_master
end






























ele=10_000_000;
using SharedArrays

arr=SharedArray{Float32,1}(ele;pids=workers());

for i in 1:ele
	arr[i]=rand()
end

M=30
mov_avg_db_par(arr,M);


