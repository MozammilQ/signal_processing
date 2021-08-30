
path="/home/red/Desktop/Julia_Workspaces/sig_proc"
include("./lpf.jl")
include("./resampler.jl")

file_f="NOAA_15_1376200_FFT_120"

sig_samp_rate=6_000_000
f_handle=open(file_f,"r")
time_dur=1
buffer=Array{Complex{Float64},1}(undef,sig_samp_rate*time_dur);

pre_demod_lpf_f=34_000
pre_demod_lpf_trans_bw=100
post_demod_lpf_f=15_000
post_demod_lpf_trans_bw=100
audio_samp_rate=44100

read!(f_handle,buffer)

var_1=lpf(buffer,pre_demod_lpf_f,pre_demod_lpf_trans_bw;input_domain=:frequency,output_domain=:time,output_samp_rate=:nyquist)
fm_data=polar_descr(var_1)

fm_data_fft=fft(fm_data)
#fm_data_fft[1]=Complex(0.0,0.0)

var_2=lpf(fm_data_fft,post_demod_lpf_f,post_demod_lpf_trans_bw;input_domain=:frequency,output_domain=:time,output_samp_rate=:nyquist)
a_data=angle.(var_2)



audio_data=Float64[]
time_sec=20.0

for i in 1:time_sec
	read!(f_handle,buffer)
	var_1=lpf(buffer,pre_demod_lpf_f,pre_demod_lpf_trans_bw;input_domain=:frequency,output_domain=:time,output_samp_rate=:nyquist)
	fm_data=polar_descr(var_1)
	append!(audio_data,fm_data)

end

path="/home/red/Desktop/Julia_Workspaces/sig_proc"
cd("./PortAudio.jl/")
using Pkg
Pkg.activate(pwd())
using PortAudio
audio_stream=PortAudioStream("HDA Intel: ALC662 rev3 Analog (hw:0,0)")

length(audio_data)
audio_data=decimate(interpolate(audio_data,2),3)
length(audio_data)
write(audio_stream,audio_data)


