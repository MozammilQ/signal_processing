


####################################################################################################
################################     Declerations on Main node    ##################################
####################################################################################################

using Distributed
worker_ip="192.168.1.2"
worker_user="root"
num_cores=4
node_working_dir="/home/red/Desktop/Julia_Workspaces/sig_proc"
file_name="simple_fm_91.9MHz"

node_file=2
node_fft=3
node_psd=4
node_lpf=5

samp_rate=Int32(6000000)
center_freq=Int32(91900000)
frame_rate=Int8(25)
num_of_samps=Int64(samp_rate * (1/frame_rate))
decimation_ratio=Int8(16)
M=Int8(5)
num_of_ticks=Int8(50)
workers_specs=("$worker_user@$worker_ip",num_cores)
machines=Vector{typeof(workers_specs)}(undef,1)
machines[1]=workers_specs
workers_pids=addprocs(machines,max_parallel=num_cores,dir=node_working_dir)

###################################################################################################



##################################################################################################
#####################        Initilizations of worker nodes        ###############################
##################################################################################################

###############################        Shared Arrays        #####################################
@everywhere procs() using SharedArrays
shared_IQ=SharedArray{Complex{Float32},1}(num_of_samps;init=false,pids=[node_file,node_fft])
shared_fft=SharedArray{Complex{Float32},1}(num_of_samps;init=false,pids=[node_fft,node_psd,node_lpf])
#################################################################################################

#################################################################################################
#   node_file
#   This node will read file, compose IQ and update shared IQ buffer 
#################################################################################################

@everywhere [node_file] num_of_samps=$num_of_samps

@everywhere [node_file] begin

	function iq_file_open(file_name::String)
		global file_handle=open(file_name,"r")
		@info "File opened!"
	end

	function update_IQ_data!()
		raw_iq_buffer=Array{Float32}(undef,2,num_of_samps)
		read!(file_handle,raw_iq_buffer)
		@inbounds global shared_IQ=Complex.(raw_iq_buffer[1,1:end],raw_iq_buffer[2,1:end])
	end
	
	function iq_file_close()
		close(file_handle)
		@info "File closed!"
	end
end

###################################################################################################


###################################################################################################
# node_fft
# The node will update shared_fft with fft of data from node_file
###################################################################################################

@everywhere [node_fft] using FFTW
@everywhere [node_fft] update_fft!()=global shared_fft=fft(shared_IQ)

####################################################################################################


####################################################################################################
# node_psd
# This node will return moving average of signal passed
####################################################################################################
@everywhere [node_psd] num_of_samps=$num_of_samps
@everywhere [node_psd] decimation_ratio=$decimation_ratio
@everywhere [node_psd] M=$M
@everywhere [node_psd] begin

	function get_psd(shared_data::SharedArray{Complex{Float32},1})

		@inbounds psd_iq=(abs2.(shared_data) ./ num_of_samps)[1:decimation_ratio:end]

		len_psd_iq=length(psd_iq)

		@inbounds psd_iq[1:1:Int(len_psd_iq/2)]=psd_iq[Int(len_psd_iq/2):-1:1]
		@inbounds psd_iq[Int(len_psd_iq/2)+1:1:end]=psd_iq[end:-1:Int(len_psd_iq/2)+1]

		db_psd_len=len_psd_iq-M+1
		db_psd=Vector{Float32}(undef,db_psd_len)

		for i=1:1:db_psd_len
			avrg_psd=0.0f0
			for j=i:1:i+M-1
				@inbounds avrg_psd+=psd_iq[i]
			end
			avrg_psd/=M
			avrg_psd+=1.0f0
			@inbounds db_psd[i]=10*log10(avrg_psd)
		end
		
		@inbounds append!(db_psd,[0.0f0 for i in 1:1:Int((M-1)/2)])
		@inbounds prepend!(db_psd,[0.0f0 for i in 1:1:Int((M-1)/2)])

		return db_psd
	end
end
####################################################################################################





###################################################################################################
##########################      Main node executions       ########################################
###################################################################################################

tickrange_x=Array{Float32}(undef,1,num_of_ticks)
ticklabels_x=Array{String}(undef,1,num_of_ticks)
ticklabels_x.=""
x_label=""

ticklabels_y=Array{String}(undef,1,num_of_ticks)
ticklabels_y.=""
tickrange_y=Array{Float32}(undef,1,num_of_ticks)

freq_rec = range(0,num_of_samps-1,step=1) .* (samp_rate / num_of_samps)
freq_rec = freq_rec .- (samp_rate/2) .+ center_freq
freq_rec=freq_rec[1:decimation_ratio:end]
freq_rec=vec(freq_rec)

tickrange_x=range(minimum(freq_rec),maximum(freq_rec),length=num_of_ticks)
for i=1:num_of_ticks
	global ticklabels_x[i]=string(round((tickrange_x[i]/1E+06), digits=1))
end
if tickrange_x[1] >= 1E+06
	global x_label="frequency(MHz)"
end
ticklabels_x=vec(ticklabels_x)






@spawnat node_file iq_file_open(file_name)
@spawnat node_file update_IQ_data!()
@spawnat node_fft update_fft!()
psd_fft=@fetchfrom node_psd get_psd(shared_fft)

println(psd_fft)


@fetchfrom node_file iq_file_close()

