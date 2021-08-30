
cd("/home/red/Desktop/Julia_Workspaces/sig_proc")


file_t="/run/media/red/21ccb3c5-d656-4b5b-9580-1dc35f680433/202011242128/202011161849_NOAA_15_raw_1376200.sdr_data"
file_f="NOAA_15_1376200_FFT_120_s"

#file_t="simple_fm_91.9MHz"
#file_f="simple_fm_91.9MHz_FFT_120"


file_t_handle=open(file_t,"r")
file_f_handle=open(file_f,"w")
using FFTW

samp_rate=6000000
time_sec=120
start_sec=60*5

seek(file_t_handle,(samp_rate*start_sec))



num_of_samps=samp_rate

buffer=Array{Float32,2}(undef,2,num_of_samps)

for i in 1:1:time_sec
	read!(file_t_handle,buffer)
	IQ_Complex=Complex.(Float64.(buffer[1,1:end]),Float64.(buffer[2,1:end]))
	write(file_f_handle,fft(IQ_Complex))
end
close(file_t_handle)
close(file_f_handle)
@info "Done converting to fft"


