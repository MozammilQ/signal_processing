

path="/home/red/Desktop/Julia_Workspaces/sig_proc"
include("./lpf.jl")
include("./resampler.jl")

file_f="simple_fm_91.9MHz_FFT_120"

sig_samp_rate=6_000_000
f_handle=open(file_f,"r")
time_dur=1
buffer=Array{Complex{Float64},1}(undef,sig_samp_rate*time_dur);

pre_demod_lpf_f=75_000
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
	fm_data=de_emph(fm_data)

	fm_data_fft=fft(fm_data)

	var_2=lpf(fm_data_fft,post_demod_lpf_f,post_demod_lpf_trans_bw;input_domain=:frequency,output_domain=:time,output_samp_rate=:nyquist)
	a_data=angle.(var_2)
	append!(audio_data,resample(a_data,audio_samp_rate))
end




path="/home/red/Desktop/Julia_Workspaces/sig_proc"
cd("./PortAudio.jl/")
using Pkg
Pkg.activate(pwd())
using PortAudio
audio_stream=PortAudioStream("HDA Intel: ALC662 rev3 Analog (hw:0,0)")





length(audio_data)
#audio_data=interpolate_gradual(audio_data,2)
write(audio_stream,audio_data)















path="/home/red/Desktop/Julia_Workspaces/sig_proc"
include("./lpf.jl")
include("./resampler.jl")

file_f="simple_fm_91.9MHz_FFT_120"
sig_samp_rate=6_000_000
f_handle=open(file_f,"r")
time_dur=1
buffer=Array{Complex{Float64},1}(undef,sig_samp_rate*time_dur)

pre_demod_lpf_f=75_000
pre_demod_lpf_trans_bw=100
post_demod_lpf_f=15_000
post_demod_lpf_trans_bw=100
audio_samp_rate=44100



include("/home/red/Desktop/Julia_Workspaces/sig_proc/lpf.jl")

buffer=Array{Float32,1}(undef,48000)
f_h=open("../raw_audio_float32","r")
include("$path/up_down_sampler.jl")

for i in 1:10
	read!(f_h,buffer)
	buffer_64=Float64.(buffer)
	write(audio_stream,angle.(ifft(lpf_t(fft(buffer_64),length(buffer_64),15_000,100))))
end





close(f_handle)



freq=200
audio_samp_rate=44100
time_sec=2
dt=1/audio_samp_rate
t_range=dt:dt:time_sec
audio_data_sn=sin.(2pi*freq.*t_range)

cd("./PortAudio.jl/")
using Pkg
Pkg.activate(pwd())
using PortAudio
audio_stream=PortAudioStream("HDA Intel: ALC662 rev3 Analog (hw:0,0)")
write(audio_stream,audio_data_sn)








scene_fft,scene_lpf=draw_psd_init(data_fft,data_lpf,length(data_fft),length(data_lpf),91900000,0)


while isopen(scene_fft)
	read!(f_handle,buffer)
	data_fft=buffer[1:1:max_freq]
	neg_f=buffer[end-max_freq+1:1:end]
	append!(data_fft,neg_f)
	data_lpf=lpf(data_fft,length(data_fft),f_p,trans_bw)
	draw_psd(scene_fft,scene_lpf,data_fft,data_lpf,length(data_fft),length(data_lpf))
end


freqs=Float64[12_500, 16_700, 20_800, 25_000, 29_200, 33_300]
samp_rate=100_000
dt=1/samp_rate
time_sec=1.0
t_range=dt:dt:time_sec
data_=Array{Float64,1}(undef,samp_rate)
data_.=0.0
for i in freqs
	data_=data_ .+ sin.(2pi*i .* t_range)
end
using FFTW
data_fft=fft(data_)
#scene_fft,scene_lpf=draw_psd_init(data_fft,data_fft,length(data_fft),length(data_fft),0,0)


filter_freq=18_000
include("./lpf.jl")
filter_data_=lpf_fft(data_fft,filter_freq)


include("psd_drawer.jl")
scene_fft,scene_lpf=draw_psd_init(data_fft,filter_data_,length(data_fft),length(filter_data_),0,0)





