using Makie
using FFTW

struct	psd_data
	num_of_graph
	graphs=LAxis[]
	data_domain
end


function psd_init(num_of_graphs::Int,data_::Array{Array{Complex{Float64},1},1},center_freqs::Array{Int64,1};data_domain=:time,title::Array{String,1})
	
	if data_domain==:time
		for i in 1:num_of_graphs
		fft!(data_[i])
	elseif data_domain==:frequency
		nothing
	else
		@info "Set the data_domain to either :time or :frequency"
		return
	end
	
	samp_rate=Array{Int64,1}(undef,num_of_graphs)
	for i in 1:num_of_graphs
		samp_rate[i]=length(data_[i])
	end
	
	num_of_ticks=25
	num_of_points=10_000

	d_ratio=Array{Int64,1}(undef,num_of_graphs)
	for i in 1:num_of_graphs
		d_ratio[i]=round(Int,samp_rate[i]/num_of_points)
	end

	M=5


	num_of_samps=Array{Int64,1}(undef,num_of_graphs)
	for i in 1:num_of_graphs
		num_of_samps[i]=samp_rate[i]
	end







	scene,layout=layoutscene(resolution=(950,950))
	scene_fft=layout[1,1]=LAxis(scene)
	scene_lpf=layout[2,1]=LAxis(scene)

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

	psd_fft=fftshift(abs2.(data_fft[1:decim:end])./num_of_samps_fft^2)
	psd_lpf=fftshift(abs2.(data_lpf[1:decim:end])./num_of_samps_lpf^2)

	for i=1:1:length(psd_fft)-M+1
		avrg_psd=0.0
		for j=i:1:i+M-1
			avrg_psd+=psd_fft[j]
		end
		avrg_psd/=M
		psd_fft[i]=avrg_psd
	end
	psd_fft[end-M+2:end].=psd_fft[end-M+1]

	minimum_psd_fft=minimum(psd_fft)
	for i in 1:1:length(psd_fft)
		psd_fft[i]=10*log10(psd_fft[i]/minimum_psd_fft)
	end
	psd_fft=vec(psd_fft)
	
	for i=1:1:length(psd_lpf)-M+1
		avrg_psd=0.0
		for j=i:1:i+M-1
			avrg_psd+=psd_lpf[j]
		end
		avrg_psd/=M
		psd_lpf[i]=avrg_psd
	end
	psd_lpf[end-M+2:end].=psd_lpf[end-M+1]

	minimum_psd_lpf=minimum(psd_lpf)
	for i in 1:1:length(psd_lpf)
		psd_lpf[i]=10*log10(psd_lpf[i]/minimum_psd_lpf)
	end
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
	
	
	fft_freq=fftshift((fftfreq(num_of_samps_fft,samp_rate_fft)[1:decim:end] .+ center_freq_fft))
	lpf_freq=fftshift((fftfreq(num_of_samps_lpf,samp_rate_lpf)[1:decim:end] .+ center_freq_lpf))
	
	tickrange_x_fft=range(minimum(fft_freq),maximum(fft_freq),length=num_of_ticks)
	if tickrange_x_fft[end] >= 1E+06
		for i=1:num_of_ticks
			ticklabels_x_fft[i]=string(round((tickrange_x_fft[i]/1E+06), digits=1))
		end
		x_label_fft="frequency(MHz)"
	elseif tickrange_x_fft[end] >= 1E+03
		for i=1:num_of_ticks
			ticklabels_x_fft[i]=string(round((tickrange_x_fft[i]/1E+03), digits=1))
		end
		x_label_fft="frequency(KHz)"
	end
	ticklabels_x_fft=vec(ticklabels_x_fft)
	
	
	tickrange_x_lpf=range(minimum(lpf_freq),maximum(lpf_freq),length=num_of_ticks)
	if tickrange_x_lpf[end] >= 1E+06
		for i=1:num_of_ticks
			ticklabels_x_lpf[i]=string(round((tickrange_x_lpf[i]/1E+06), digits=1))
		end
		x_label_lpf="frequency(MHz)"
	elseif tickrange_x_lpf[end] >= 1E+03
		for i=1:num_of_ticks
			ticklabels_x_lpf[i]=string(round((tickrange_x_lpf[i]/1E+03), digits=1))
		end
		x_label_lpf="frequency(KHz)"
	end
	ticklabels_x_lpf=vec(ticklabels_x_lpf)

	scene,layout=layoutscene(resolution=(950,950))
	scene_fft=layout[1,1]=LAxis(scene)
	scene_lpf=layout[2,1]=LAxis(scene)
	lines!(scene_fft,fft_freq,psd_fft,color=:blue,linewidth=0.5)
	lines!(scene_lpf,lpf_freq,psd_lpf,color=:black,linewidth=0.5)
	
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
	scene_fft.xlabel=x_label_fft
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
	scene_lpf.xlabel=x_label_lpf
	scene_lpf.xlabelsize=10.0f0
	scene_lpf.ylabel="Relative Power(dB)"
	scene_lpf.ylabelsize=10.0f0
	scene_lpf.ylabelpadding=4.0f0
	display(scene)

	return scene_fft,scene_lpf

end

function psd_fft_update(scene_fft,scene_lpf,data_fft,data_lpf;data_domain=:time)

	if data_domain
	samp_rate_fft=length(data_fft)
	samp_rate_lpf=length(data_lpf)
	if !isopen(scene_fft)
		@info "Error! Graph is not open\nRun psd_fft_init() to set the graph"
		return
	end
	
	M=5
	num_of_ticks=25
	num_of_samps_fft=samp_rate_fft
	num_of_samps_lpf=samp_rate_lpf
	num_of_points=10_000

	decim=round(Int,samp_rate_fft/num_of_points)

	psd_fft=fftshift(abs2.(data_fft[1:decim:end])./num_of_samps_fft^2)
	psd_lpf=fftshift(abs2.(data_lpf[1:decim:end])./num_of_samps_lpf^2)


	ticklabels_y_fft=Array{String}(undef,num_of_ticks)
	ticklabels_y_fft.=""
	ticklabels_y_lpf=Array{String}(undef,num_of_ticks)
	ticklabels_y_lpf.=""

	for i=1:1:length(psd_fft)-M+1
		avrg_psd=0.0
		for j=i:1:i+M-1
			avrg_psd+=psd_fft[j]
		end
		avrg_psd/=M
		psd_fft[i]=avrg_psd
	end
	psd_fft[end-M+2:end].=psd_fft[end-M+1]

	minimum_psd_fft=minimum(psd_fft)
	for i in 1:1:length(psd_fft)
		psd_fft[i]=10*log10(psd_fft[i]/minimum_psd_fft)
	end
	psd_fft=vec(psd_fft)
	
	for i=1:1:length(psd_lpf)-M+1
		avrg_psd=0.0
		for j=i:1:i+M-1
			avrg_psd+=psd_lpf[j]
		end
		avrg_psd/=M
		psd_lpf[i]=avrg_psd
	end
	psd_lpf[end-M+2:end].=psd_lpf[end-M+1]

	minimum_psd_lpf=minimum(psd_lpf)
	for i in 1:1:length(psd_lpf)
		psd_lpf[i]=10*log10(psd_lpf[i]/minimum_psd_lpf)
	end
	psd_lpf=vec(psd_lpf)
			
	tickrange_y_fft=range(minimum(psd_fft),maximum(psd_fft),length=num_of_ticks)
	for i=1:num_of_ticks
		ticklabels_y_fft[i]=string(round(tickrange_y_fft[i],digits=2))
	end
	ticklabels_y_fft=vec(ticklabels_y_fft)
	
	tickrange_y_lpf=range(minimum(psd_lpf),maximum(psd_lpf),length=num_of_ticks)
	for i=1:num_of_ticks
		ticklabels_y_lpf[i]=string(round(tickrange_y_lpf[i],digits=2))
	end
	ticklabels_y_lpf=vec(ticklabels_y_lpf)
	
	scene_fft.yticks=(tickrange_y_fft,ticklabels_y_fft)
	scene_lpf.yticks=(tickrange_y_lpf,ticklabels_y_lpf)
	scene_fft.scene[end][2]=psd_fft
	scene_lpf.scene[end][2]=psd_lpf
	autolimits!(scene_fft)
	autolimits!(scene_lpf)
	sleep(0.05)
end


