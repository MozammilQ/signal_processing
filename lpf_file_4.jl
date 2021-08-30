

function decimate(data_::Array{Float64,1},decim_ratio::Int64)
	return data_[1:decim_ratio:end]
end

function interpolate(data_::Array{Float64,1},inter_ratio::Int64)
	dummy_=zeros(inter_ratio-1)
	buffer=Float64[]
	for i in 1:1:length(data_)
		append!(buffer,data_[i])
		append!(buffer,dummy_)
	end
	return buffer
end

using FFTW
file_f="simple_fm_91.9MHz_FFT_120"

samp_rate=6_000_000
num_of_samps=samp_rate

f_handle=open(file_f,"r")
buffer=Array{Complex{Float64},1}(undef,num_of_samps)
max_freq=200_000
inter_r=147
decim_r=100

time_sec=10.0

fm_data=Float64[]

for i in 1:time_sec
	read!(f_handle,buffer)

	data_fft=buffer[1:1:max_freq]
	neg_f=buffer[end-max_freq+1:1:end]
	append!(data_fft,neg_f)
	
	data_t=ifft(data_fft)
	
	for i in 2:1:length(data_fft)
		dummy_=conj(data_t[i-1]) * data_t[i]
		append!(fm_data,atan(imag(dummy_)/real(dummy_)))
	end
end


volume=1.0E-1
fm_data = fm_data .* volume

cd("./PortAudio.jl/")
using Pkg
Pkg.activate(pwd())
using PortAudio
audio_stream=PortAudioStream("HDA Intel: ALC662 rev3 Analog (hw:0,0)")
write(audio_stream,fm_data)





