path="/home/red/Desktop/Julia_Workspaces/sig_proc"
include("$path/resampler.jl")
using FFTW




function de_emph(data_::Array{Float64,1};tau_::Float64=50E-6)
	len_data_=length(data_)
	x=exp(-1/(len_data_*tau_))
	b=[1-x]
	a=[1,-x]
	y=Array{Float64,1}(undef,len_data_)
	y.=0.0
	for i in 2:1:len_data_
		y[i] = ( b[1] * data_[i] ) - ( a[2]*y[i-1] )
	end
	return y
end













function lpf(data_lpf::Array{Complex{Float64},1},f_p::Int64,trans_bw::Int64;input_domain::Symbol=:time,output_domain::Symbol=:time,output_samp_rate::Symbol=:nyquist,decim_r::Int64=1)
	if input_domain==:time
		fft!(data_lpf)
	elseif input_domain==:frequency
		nothing
	else
		@info "Set input_domain to either :time or :frequency\n\tReturning to calling function"
		return
	end
	if output_domain==:time
		nothing
	elseif output_domain==:frequency
		nothing
	else
		@info "Set output_domain to :time or :frequency\n\tReturning to calling function"
		return
	end
	f_samp::Int64=length(data_lpf)
	if output_samp_rate==:nyquist
		global decim_r=round(Int,f_samp/(f_p*2))
	elseif output_samp_rate==:input_samp_rate
		nothing
	else
		@info "Set output_samp_rate to either :nyquist or :input_samp_rate\n\tReturning to calling function"
		return
	end
	f_s=f_p+trans_bw
	w_p=2pi*(f_p/f_samp)
	w_s=2pi*(f_s/f_samp)
	len_ham=round(Int,8pi/(w_s-w_p))
	if f_samp%len_ham!=0
		@info "Accuracy of filter may be compromised,\n\t Try changing transition bandwidth\n\tReturning to calling function"
		return
	end
	w_c=(w_p+w_s)/2
	h_t=Array{Float64,1}(undef,len_ham)
	for n in 1:1:len_ham
		h_t[n]=(0.54-0.46*cos(2*(n-1)*pi/len_ham))*((w_c/pi)*sinc((w_c*((n-1)-(len_ham/2)))/pi)) 
	end
	
	h_t=fft(h_t)

	h_t=interpolate(h_t,round(Int,f_samp/len_ham);method=:linear)
	#h_t=upsample(h_t,round(Int,f_samp/len_ham))
	data_lpf=h_t .* data_lpf




	if output_samp_rate==:nyquist && output_domain==:time
		return decimate(ifft(data_lpf),decim_r)
	elseif output_samp_rate==:input_samp_rate && output_domain==:time
		return ifft(data_lpf)
	elseif output_samp_rate==:nyquist && output_domain==:frequency
		return fft(decimate(ifft(data_lpf),decim_r))
	end
	return data_lpf
end

function polar_descr(data_t::Array{Complex{Float64},1})
	buffer_=Array{Float64,1}(undef,length(data_t))
	buffer_[1]=0.0
	for i in 2:1:length(data_t)
		buffer_[i]=angle(conj(data_t[i-1]) * data_t[i])
	end
	return buffer_
end

function polar_descr(data_t::Array{Complex{Float32},1})
	buffer_=Array{Float32,1}(undef,length(data_t))
	buffer_[1]=0.0f0
	for i in 2:1:length(data_t)
		buffer_[i]=angle(conj(data_t[i-1]) * data_t[i])
	end
	return buffer_
end


