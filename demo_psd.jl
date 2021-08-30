module draw_psd
using FFTW,Makie
function draw_psd(file_name,samp_rate,center_freq,frame_rate)
	time_begin=0.0f0
	time_end=0.0f0
	time_to_sleep=0.0f0
	decimation_ratio=64
	M=100
	num_of_ticks=50
	num_of_samples=Int(samp_rate * (1/frame_rate))
	raw_iq_buffer=Array{Float32}(undef,1,2*num_of_samples)
	points_to_ignore=Int(num_of_samples/100)
	I_=Array{Float32}(undef,1,num_of_samples)
	Q_=Array{Float32}(undef,1,num_of_samples)
	IQ_complex=Array{Complex{Float32}}(undef,1,num_of_samples)
	tickrange_x=Array{Float32}(undef,1,num_of_ticks)
	tickrange_y=Array{Float32}(undef,1,num_of_ticks)
	ticklabels_x=Array{String}(undef,1,num_of_ticks)
	ticklabels_x.=""
	ticklabels_y=Array{String}(undef,1,num_of_ticks)
	ticklabels_y.=""
	x_label=""
	file_handle=open(file_name,"r")
	
	raw_iq_buffer.=0.0f0
	I_.=0.0f0
	Q_.=0.0f0
	IQ_complex.=0.0f0
	
	global freq_rec = range(0,num_of_samples-1,step=1) .* (samp_rate / num_of_samples)
	global freq_rec = freq_rec .- (samp_rate/2) .+ center_freq
	global freq_rec=freq_rec[1:decimation_ratio:end]
	freq_rec=vec(freq_rec)
	
	global avrg_psd_len=length(freq_rec)
	global avrg_psd=Array{Float32}(undef,1,avrg_psd_len)
	avrg_psd.=0.0f0
	avrg_psd=vec(avrg_psd)
	
	read!(file_handle,raw_iq_buffer)
	
	I_=raw_iq_buffer[1:2:end]
	Q_=raw_iq_buffer[2:2:end]
	
	for i=1:num_of_samples
		IQ_complex[i]=I_[i] + Q_[i]im
	end
	
	fft_IQ_complex = fft(IQ_complex)
	PSD_IQ_complex = abs2.(fft_IQ_complex) ./ num_of_samples
	
	for i=1:length(PSD_IQ_complex)
		if isnan(PSD_IQ_complex[i])
			PSD_IQ_complex[i]=0
		end
	end
	
	PSD_IQ_complex[1:1:points_to_ignore] .= 0.0f0
	PSD_IQ_complex[end-points_to_ignore:1:end] .= 0.0f0
	PSD_IQ_complex=PSD_IQ_complex[1:decimation_ratio:end]
	
	for i=1:1:avrg_psd_len-M+1
		for j=i:1:i+M-1
			avrg_psd[i]+=PSD_IQ_complex[j]
		end
	end
	
	tickrange_x=range(minimum(freq_rec),maximum(freq_rec),length=num_of_ticks)
	for i=1:num_of_ticks
		global ticklabels_x[i]=string(round((tickrange_x[i]/1E+06), digits=1))
	end
	if tickrange_x[1] >= 1E+06
		global x_label="frequency(MHz)"
	end
	ticklabels_x=vec(ticklabels_x)
	
	tickrange_y=range(minimum(avrg_psd),maximum(avrg_psd),length=num_of_ticks)
	global y_label="Relative Power(dB)"
	for i=1:num_of_ticks
		ticklabels_y[i]=string(round(tickrange_y[i],digits=2))
	end
	ticklabels_y=vec(ticklabels_y)
	
	global scene=Scene(resolution=(1200,600))
	lines!(scene,freq_rec,avrg_psd,color=:blue,linewidth=0.5)
	global lineplot=scene[end]
	global axis = scene[Axis]
	axis.grid.linecolor = ((:black, 0.5), (:black, 0.5))
	axis.names.textcolor = ((:black, 1.0), (:black, 1.0))
	axis.names.axisnames = (x_label,y_label)
	axis.names.textsize = (2.5,2.5)
	axis.ticks.textsize=(1.5,1.5)
	AbstractPlotting.ticks!(scene,tickranges=(tickrange_x,tickrange_y))
	AbstractPlotting.ticks!(scene,ticklabels=(ticklabels_x,ticklabels_y))
	avrg_psd.=0.0f0
	display(scene)
		
	for i=1:1:20
		global time_begin=time()
		read!(file_handle,raw_iq_buffer)
		
		global I_=raw_iq_buffer[1:2:end]
		global Q_=raw_iq_buffer[2:2:end]
		
		for i=1:num_of_samples
			IQ_complex[i]=I_[i] + Q_[i]im
		end
		
		global fft_IQ_complex = fft(IQ_complex)
		global PSD_IQ_complex = abs2.(fft_IQ_complex) ./ num_of_samples
		
		for i=1:length(PSD_IQ_complex)
			if isnan(PSD_IQ_complex[i])
				PSD_IQ_complex[i]=0
			end
		end
		
		PSD_IQ_complex[1:1:points_to_ignore] .= 0.0f0
		PSD_IQ_complex[end-points_to_ignore:1:end] .= 0.0f0
		
		PSD_IQ_complex=PSD_IQ_complex[1:decimation_ratio:end]
		
		for i=1:1:avrg_psd_len-M+1
			for j=i:1:i+M-1
				avrg_psd[i]+=PSD_IQ_complex[j]
			end
		end
		global tickrange_y=range(minimum(avrg_psd),maximum(avrg_psd),length=num_of_ticks)
		for i=1:num_of_ticks
			ticklabels_y[i]=string(round(tickrange_y[i],digits=2))
		end
		ticklabels_y=vec(ticklabels_y)
	
		AbstractPlotting.ticks!(scene,tickranges=(tickrange_x,tickrange_y))
		AbstractPlotting.ticks!(scene,ticklabels=(ticklabels_x,ticklabels_y))
		lineplot[2]=avrg_psd
		AbstractPlotting.update_limits!(scene)
		AbstractPlotting.update!(scene)
		avrg_psd.=0.0f0
		global time_end=time()
		global time_to_sleep=((1/frame_rate)-(time_end-time_begin))
		if time_to_sleep <= 0
			@info "System is experiencing performance throttle"
			time_to_sleep*=-1
		end
		sleep(time_to_sleep)
	end
	
	
	
	close(file_handle)
end
end
