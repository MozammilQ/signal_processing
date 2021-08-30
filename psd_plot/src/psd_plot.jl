module psd_plot
	using FFTW,Makie

	@doc """
		this is a multiline doc
		see
		multiline

	"""->

	function draw_psd(file_name::String,samp_rate::Int,center_freq::Int)
		

		#using FFTW
		#file_name="simple_fm_91.9MHz_FFT_120"
		#samp_rate=6000000
		#center_freq=91900000
		decim=400
		M=5
		num_of_ticks=50
		num_of_samps=samp_rate
		
		buffer=Array{Complex{Float64},1}(undef,num_of_samps)
		f_handle=open(file_name,"r")
		read!(f_handle,buffer)

		tickrange_x=Array{Float64}(undef,1,num_of_ticks)
		tickrange_y=Array{Float64}(undef,1,num_of_ticks)
		
		ticklabels_x=Array{String}(undef,1,num_of_ticks)
		ticklabels_x.=""
		ticklabels_y=Array{String}(undef,1,num_of_ticks)
		ticklabels_y.=""
		
		x_label=""
		
		psd=fftshift(abs2.(buffer[1:decim:end])./num_of_samps^2)
		fft_freq=fftshift((fftfreq(num_of_samps,samp_rate)[1:decim:end] .+ center_freq))

		for i=1:1:length(psd)-M+1
			avrg_psd=0.0
			for j=i:1:i+M-1
				avrg_psd+=psd[j]
			end
			avrg_psd/=M
			avrg_psd+=1.0
			psd[i]=10*log10(avrg_psd)			
		end
		psd[end-M+2:end].=minimum(psd)
		psd=vec(psd)
	
		tickrange_x=range(minimum(fft_freq),maximum(fft_freq),length=num_of_ticks)
		for i=1:num_of_ticks
			 ticklabels_x[i]=string(round((tickrange_x[i]/1E+06), digits=1))
		end
		if tickrange_x[end] >= 1E+06
			 x_label="frequency(MHz)"
		end
		ticklabels_x=vec(ticklabels_x)

		tickrange_y=range(minimum(psd),maximum(psd),length=num_of_ticks)
		y_label="Relative Power(dB)"
		for i=1:num_of_ticks
			ticklabels_y[i]=string(round(tickrange_y[i],digits=2))
		end
		ticklabels_y=vec(ticklabels_y)
		
		#using Makie
		@info "First, run will take some time, depending on system's performance, please be patient!"
		scene=Scene(resolution=(1200,600))
		@info "Just there, be patient"
		lines!(scene,fft_freq,psd,color=:blue,linewidth=0.5)
		lineplot=scene[end]
		axis = scene[Axis]
		axis.grid.linecolor = ((:black, 0.5), (:black, 0.5))
		axis.names.textcolor = ((:black, 1.0), (:black, 1.0))
		axis.names.axisnames = (x_label,y_label)
		axis.names.textsize = (2.5,2.5)
		axis.ticks.textsize=(1.5,1.5)
		AbstractPlotting.ticks!(scene,tickranges=(tickrange_x,tickrange_y))
		AbstractPlotting.ticks!(scene,ticklabels=(ticklabels_x,ticklabels_y))
		display(scene)

		for i=1:1:40

			read!(f_handle,buffer)
		
			psd=fftshift(abs2.(buffer[1:decim:end])./num_of_samps^2)
	
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


	end
end


