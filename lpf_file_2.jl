using Makie
using FFTW

function decimate(data_::Array{Float64,1},decim_ratio::Int64)
        return data_[1:decim_ratio:end]
end

function interpolate_zero(data_::Array{Float64,1},inter_ratio::Int64)
	if inter_ratio==1
		return data_
	end
        dummy_=zeros(inter_ratio-1)
        buffer=Float64[]
        for i in 1:1:length(data_)
                append!(buffer,data_[i])
                append!(buffer,dummy_)
        end
        return buffer
end

function interpolate_same(data_::Array{Float64,1},inter_ratio::Int)
	if inter_ratio==1
		return data_
	end
	len_data=length(data_)
	dummy_=Array{Float64,1}(undef,len_data*inter_ratio)
	counter_::Int=1
	for i in 1:1:len_data
		dummy_[counter_]=data_[i]
		counter_+=1
		for j in 1:1:inter_ratio-1
			dummy_[counter_]=data_[i]
			counter_+=1
		end
	end
	return dummy_
end

function interpolate_gradual(data_::Array{Float64,1},inter_ratio::Int)
	if inter_ratio==1
		return data_
	end
	counter_::Int64=1
	len_data=length(data_)
	last_::Float64=0.0
	dummy_=Array{Float64,1}(undef,len_data*inter_ratio)
	for i in 1:1:len_data
		dummy_[counter_]=data_[i]
		counter_+=1
		dt=(data_[i]-last_)/inter_ratio
		for j in 1:1:inter_ratio-1
			dummy_[counter_]=data_[i] + (dt*j)
			counter_+=1
		end
	last_=data_[i]
	end
	return dummy_
end
	

function interpolate_same(data_::Array{Complex{Float64},1},inter_ratio::Int)
	@info "Interpolating Complex data"
	if inter_ratio==1
		return data_
	end
	len_data=length(data_)
	dummy_=Array{Complex{Float64},1}(undef,len_data*inter_ratio)
	counter::Int=1
	for i in 1:1:len_data
		dummy_[counter]=data_[i]
		counter+=1
		for j in 1:1:inter_ratio-1
			dummy_[counter]=data_[i]
			counter+=1
		end
	end
	return dummy_
end


function interpolate_gradual(data_::Array{Complex{Float64},1},inter_ratio::Int)
	@info "Interpolating Complex data"
	if inter_ratio==1
		return data_
	end
	dummy_=Array{Complex{Float64},1}(undef,length(data_)*inter_ratio)
	last_cm=Complex(0.0,0.0)
	dummy_.=last_cm
	counter::Int=1
	for i in 1:1:length(data_)
		dummy_[counter]=data_[i]
		counter+=1
		real_data_=real(data_[i])
		imag_data_=imag(data_[i])
		real_diff=real_data_-real(last_cm)
		imag_diff=imag_data_-imag(last_cm)
		len_loop=inter_ratio
		dt_real=real_diff/len_loop
		dt_imag=imag_diff/len_loop
		for j in 1:1:len_loop-1
			dummy_[counter]=Complex( (real_data_+(dt_real*j)), (imag_data_+(dt_imag*j)) )
			counter+=1
		end
		last_cm=data_[i]
	end
	return dummy_
end


function lpf_fft(data_fft::Array{Complex{Float64},1},f_p::Int64)
	f_=data_fft[1:1:f_p]
	neg_f=data_fft[end-f_p+1:1:end]
	append!(f_,neg_f)
	return f_
end

function polar_descr(data_t::Array{Complex{Float64}})
	len_data_t=length(data_t)
	buffer_=Array{Float64,1}(undef,len_data_t)
	buffer_[1]=0.0
	for i in 2:1:len_data_t
		dummy_=conj(data_t[i-1]) * data_t[i]
		buffer_[i]=atan(imag(dummy_)/real(dummy_))
	end
	return buffer_
end

function lpf_t(data_fft,samp_rate,f_p,trans_bw)

	f_samp=samp_rate
	num_of_samps=samp_rate
	f_s=f_p+trans_bw
	w_p=2pi*(f_p/f_samp)
	w_s=2pi*(f_s/f_samp)
	
	# Applying Hamming window
	len_ham=round(Int,8pi/(w_s-w_p))

	w_c=(w_p+w_s)/2
	h_lpf=Array{Float64,1}(undef,len_ham)
	
	# Constructing Hamming Window
	for n in 1:1:len_ham
		h_lpf[n]=(0.54-0.46*cos(2*(n-1)*pi/len_ham))*((w_c/pi)*sinc((w_c*((n-1)-(len_ham/2)))/pi)) 
	end
	
	h_fft=fft(h_lpf)

	ratio_=round(Int,samp_rate/len_ham)
	
	lpf_fft=Array{Complex{Float64},1}(undef,samp_rate)

	for i in 1:1:len_ham
		start_=(ratio_*(i-1))+1
		end_=start_+ratio_-1
		lpf_fft[start_:1:end_].=h_fft[i] 
	end

	return lpf_fft .* data_fft
end

function draw_psd_init(data_fft,data_lpf,samp_rate_fft,samp_rate_lpf,center_freq_fft,center_freq_lpf)

	num_of_ticks=25
	num_of_points=10_000

	decim=round(Int,samp_rate_fft/num_of_points)
	M=5
	num_of_samps_fft=samp_rate_fft
	num_of_samps_lpf=samp_rate_lpf

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

function draw_psd(scene_fft,scene_lpf,data_fft,data_lpf,samp_rate_fft,samp_rate_lpf)
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
	sleep(0.05)
end


