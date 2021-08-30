

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


using SharedArrays
shared_IQ=SharedArray{Complex{Float32},1}(num_of_samps)
shared_fft=SharedArray{Complex{Float32},1}(num_of_samps)


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


using FFTW
update_fft!()=global shared_fft=fft(shared_IQ)

function get_psd(shared_data)
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


