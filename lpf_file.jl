

file_f="simple_fm_91.9MHz_FFT_120"

samp_rate=6_000_000
f_p=75_000
trans_bw=10_000


f_samp=samp_rate
num_of_samps=samp_rate
f_s=f_p+trans_bw
w_p=2pi*(f_p/f_samp)
w_s=2pi*(f_s/f_samp)
len_ham=round(Int,8pi/(w_s-w_p))
w_c=(w_p+w_s)/2

h_lpf=Array{Float64,1}(undef,len_ham)

for n in 1:1:len_ham
	h_lpf[n]= ( 0.54 - 0.46*cos( 2*(n-1)*pi/len_ham ) ) * ( (w_c/pi) * sinc(  ( w_c*((n-1)-(len_ham/2)) )/pi   ) ) 
end

using FFTW

h_fft=fft(h_lpf)


center_freq_fft=91900000
center_freq_lpf=0
decim=800
M=5
num_of_ticks=25

tickrange_x_fft=Array{Float64}(undef,1,num_of_ticks)
ticklabels_x_fft=Array{String}(undef,1,num_of_ticks)
ticklabels_x_fft.=""
x_label_fft=""

tickrange_x_lpf=Array{Float64}(undef,1,num_of_ticks)
ticklabels_x_lpf=Array{String}(undef,1,num_of_ticks)
ticklabels_x_lpf.=""
x_label_lpf=""

tickrange_y_fft=Array{Float64}(undef,1,num_of_ticks)	
ticklabels_y_fft=Array{String}(undef,1,num_of_ticks)
ticklabels_y_fft.=""

tickrange_y_lpf=Array{Float64}(undef,1,num_of_ticks)	
ticklabels_y_lpf=Array{String}(undef,1,num_of_ticks)
ticklabels_y_lpf.=""


f_handle=open(file_f,"r")
buffer=Array{Complex{Float64},1}(undef,num_of_samps)
read!(f_handle,buffer)


reduced_fft=buffer[1:1:round(Int,num_of_samps-(num_of_samps%len_ham))]
lpf_fft=Array{Complex{Float64},1}(undef,length(reduced_fft))

for i in 1:1:round(Int,(length(reduced_fft)/len_ham))
	start_=(len_ham*(i-1))+1
	end_=start_+len_ham-1
	lpf_fft[start_:1:end_]=h_fft .* reduced_fft[start_:1:end_] 
end
psd_fft=fftshift(abs2.(reduced_fft[1:decim:end])./num_of_samps^2)
psd_lpf=fftshift(abs2.(lpf_fft[1:decim:end])./num_of_samps^2)

for i=1:1:length(psd_fft)-M+1
	avrg_psd=0.0
	for j=i:1:i+M-1
		avrg_psd+=psd_fft[j]
	end
	avrg_psd/=M
	avrg_psd+=1.0
	psd_fft[i]=10*log10(avrg_psd)			
end
psd_fft[end-M+2:end].=psd_fft[end-M+1]
psd_fft=vec(psd_fft)

for i=1:1:length(psd_lpf)-M+1
	avrg_psd=0.0
	for j=i:1:i+M-1
		avrg_psd+=psd_lpf[j]
	end
	avrg_psd/=M
	avrg_psd+=1.0
	psd_lpf[i]=10*log10(avrg_psd)			
end
psd_lpf[end-M+2:end].=psd_lpf[end-M+1]
psd_lpf=vec(psd_lpf)

tickrange_y_fft=range(minimum(psd_fft),maximum(psd_fft),length=num_of_ticks)
y_label="Relative Power(dB)"
for i=1:num_of_ticks
	ticklabels_y_fft[i]=string(round(tickrange_y_fft[i],digits=2))
end
ticklabels_y_fft=vec(ticklabels_y_fft)

tickrange_y_lpf=range(minimum(psd_lpf),maximum(psd_lpf),length=num_of_ticks)
y_label="Relative Power(dB)"
for i=1:num_of_ticks
	ticklabels_y_lpf[i]=string(round(tickrange_y_lpf[i],digits=2))
end
ticklabels_y_lpf=vec(ticklabels_y_lpf)


fft_freq=fftshift((fftfreq(num_of_samps,samp_rate)[1:decim:end] .+ center_freq_fft))
lpf_freq=fftshift((fftfreq(num_of_samps,samp_rate)[1:decim:end] .+ center_freq_lpf))

tickrange_x_fft=range(minimum(fft_freq),maximum(fft_freq),length=num_of_ticks)
for i=1:num_of_ticks
	 ticklabels_x_fft[i]=string(round((tickrange_x_fft[i]/1E+06), digits=1))
end
if tickrange_x_fft[end] >= 1E+06
	 x_label_fft="frequency(MHz)"
elseif tickrange_x_fft[end] >= 1E+03
	x_label_fft="frequency(KHz)"
end
ticklabels_x_fft=vec(ticklabels_x_fft)


tickrange_x_lpf=range(minimum(lpf_freq),maximum(lpf_freq),length=num_of_ticks)
for i=1:num_of_ticks
	 ticklabels_x_lpf[i]=string(round((tickrange_x_lpf[i]/1E+06), digits=1))
end
if tickrange_x_lpf[end] >= 1E+06
	 x_label_fft="frequency(MHz)"
elseif tickrange_x_lpf[end] >= 1E+03
	x_label_fft="frequency(KHz)"
end
ticklabels_x_lpf=vec(ticklabels_x_lpf)



using Makie
@info "First, run will take some time, depending on system's performance, please be patient!"


scene,layout=layoutscene(resolution=(950,950))
scene_fft=layout[1,1]=LAxis(scene)
scene_lpf=layout[2,1]=LAxis(scene)


lines!(scene_fft,fft_freq,psd_fft,color=:blue,linewidth=0.5)
lines!(scene_lpf,lpf_freq,psd_lpf,color=:black,linewidth=0.5)

lineplot_fft=scene_fft.scene[end]
#lineplot_fft[1] for x and 2 for y
# then
# AbstractPlotting.update_limits!(scene)
# and,   AbstractPlotting.update!(scene)



scene_fft.title="FFT of signal"
scene_fft.titlesize=12.0f0
scene_fft.titlegap=4.0f0
scene_fft.titlealign=:left
scene_fft.xticks=(tickrange_x_fft,ticklabels_x_fft)
scene_fft.xticksize=5.0f0
scene_fft.xticklabelsize=8.0f0
scene_fft.xlabelpadding=4.0f0
scene_fft.yticks=(tickrange_y_fft,ticklabels_y_fft)
scene_fft.yticksize=5.0f0
scene_fft.yticklabelsize=8.0f0
scene_fft.xlabel="Frequency(MHz)"
scene_fft.xlabelsize=10.0f0
scene_fft.ylabel="Relative Power(dB)"
scene_fft.ylabelsize=10.0f0
scene_fft.ylabelpadding=4.0f0

scene_lpf.title="FFT of signal after LPF"
scene_lpf.titlesize=12.0f0
scene_lpf.titlegap=4.0f0
scene_lpf.titlealign=:left
scene_lpf.xticks=(tickrange_x_lpf,ticklabels_x_lpf)
scene_lpf.xticksize=5.0f0
scene_lpf.xticklabelsize=8.0f0
scene_lpf.xlabelpadding=4.0f0
scene_lpf.yticks=(tickrange_y_lpf,ticklabels_y_lpf)
scene_lpf.yticksize=5.0f0
scene_lpf.yticklabelsize=8.0f0
scene_lpf.xlabel="Frequency(MHz)"
scene_lpf.xlabelsize=10.0f0
scene_lpf.ylabel="Relative Power(dB)"
scene_lpf.ylabelsize=10.0f0
scene_lpf_ylabelpadding=4.0f0
sleep(0.05)
display(scene)



















lineplot=scene_fft[end]
axis = scene_fft[Axis]
axis.grid.linecolor = ((:black, 0.5), (:black, 0.5))
axis.names.textcolor = ((:black, 1.0), (:black, 1.0))
axis.names.axisnames = (x_label,y_label)
axis.names.textsize = (2.5,2.5)
axis.ticks.textsize=(1.5,1.5)
AbstractPlotting.ticks!(scene_fft,tickranges=(tickrange_x,tickrange_y_fft))
AbstractPlotting.ticks!(scene_fft,ticklabels=(ticklabels_x,ticklabels_y_fft))
display(scene_fft)



lin=lines!(scene_lpf,fft_freq,psd_lpf,color=:blue,linewidth=0.5)
lineplot=scene_lpf[end]
axis = scene_lpf[Axis]
axis.grid.linecolor = ((:black, 0.5), (:black, 0.5))
axis.names.textcolor = ((:black, 1.0), (:black, 1.0))
axis.names.axisnames = (x_label,y_label)
axis.names.textsize = (2.5,2.5)
axis.ticks.textsize=(1.5,1.5)
AbstractPlotting.ticks!(scene_lpf,tickranges=(tickrange_x,tickrange_y_lpf))
AbstractPlotting.ticks!(scene_lpf,ticklabels=(ticklabels_x,ticklabels_y_lpf))
display(scene_lpf)


for i=1:1:40
	read!(f_handle,buffer)
	
	reduced_fft=buffer[1:1:round(Int,num_of_samps-(num_of_samps%len_ham))]
	for i in 1:1:round(Int,(length(reduced_fft)/len_ham))
		start_=(len_ham*(i-1))+1
		end_=start_+len_ham-1
		reduced_fft[start_:1:end_]=h_fft .* reduced_fft[start_:1:end_]
	end
	
	psd=fftshift(abs2.(reduced_fft[1:decim:end])./num_of_samps^2)

	for i=1:1:length(psd)-M+1
		avrg_psd=0.0
		for j=i:1:i+M-1
			@inbounds avrg_psd+=psd[j]
		end
		avrg_psd/=M
		avrg_psd+=1.0
		@inbounds psd[i]=10*log10(avrg_psd)
	end
	
	psd[end-M+2:end].=psd[end-M+1]
	psd=vec(psd)

	tickrange_y=range(minimum(psd),maximum(psd),length=num_of_ticks)
	
	for i=1:num_of_ticks
		@inbounds ticklabels_y[i]=string(round(tickrange_y[i],digits=2))
	end

	ticklabels_y=vec(ticklabels_y)
	AbstractPlotting.ticks!(scene,tickranges=(tickrange_x,tickrange_y))
	AbstractPlotting.ticks!(scene,ticklabels=(ticklabels_x,ticklabels_y))
	lineplot[2]=psd
	AbstractPlotting.update_limits!(scene)
	AbstractPlotting.update!(scene)
	sleep(0.05)

end

close(f_handle)


