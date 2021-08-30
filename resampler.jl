function decimate(data_::Array{Float64,1},decim_ratio::Int64)
        return data_[1:decim_ratio:end]
end

function decimate(data_::Array{Float32,1},decim_ratio::Int64)
        return data_[1:decim_ratio:end]
end

function decimate(data_::Array{Complex{Float64},1},decim_ratio::Int64)
	return data_[1:decim_ratio:end]
end

function decimate(data_::Array{Complex{Float32},1},decim_ratio::Int64)
	return data_[1:decim_ratio:end]
end

function upsample(data_::Array{Float64,1},inter_ratio::Int)
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

function upsample(data_::Array{Float32,1},inter_ratio::Int)
	if inter_ratio==1
		return data_
	end
	len_data=length(data_)
	dummy_=Array{Float32,1}(undef,len_data*inter_ratio)
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

function upsample(data_::Array{Complex{Float64},1},inter_ratio::Int)
	if inter_ratio==1
		return data_
	end
	len_data=length(data_)
	dummy_=Array{Complex{Float64},1}(undef,len_data*inter_ratio)
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

function upsample(data_::Array{Complex{Float32},1},inter_ratio::Int)
	if inter_ratio==1
		return data_
	end
	len_data=length(data_)
	dummy_=Array{Complex{Float32},1}(undef,len_data*inter_ratio)
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

function interpolate(data_::Array{Float64,1},inter_ratio::Int64;method::Symbol=:linear)
	if inter_ratio==1
		return data_
	end

	if method==:linear
		begin
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

	elseif method==:zeros
        	buffer=Array{Float64,1}(undef,length(data_)*inter_ratio)
		counter_=1
	        for i in 1:1:length(data_)
			buffer[counter_]=data_[i]
			counter_+=1
			for j in 1:1:inter_ratio-1
				buffer[counter_]=0.0
				counter_+=1               	
		        end
		end
	        return buffer
	
	else
		@info "Set method to one of the following :linear,:zeros\n\tReturning to the calling function"
		return
	end
end

function interpolate(data_::Array{Float32,1},inter_ratio::Int64;method::Symbol=:linear)
	if inter_ratio==1
		return data_
	end

	if method==:linear
		begin
		counter_::Int64=1
		len_data=length(data_)
		last_::Float32=0.0
		dummy_=Array{Float32,1}(undef,len_data*inter_ratio)
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

	elseif method==:zeros
        	buffer=Array{Float32,1}(undef,length(data_)*inter_ratio)
		counter_=1
	        for i in 1:1:length(data_)
			buffer[counter_]=data_[i]
			counter_+=1
			for j in 1:1:inter_ratio-1
				buffer[counter_]=0.0
				counter_+=1               	
		        end
		end
	        return buffer
	
	else
		@info "Set method to one of the following :linear,:zeros\n\tReturning to the calling function"
		return
	end
end

function interpolate(data_::Array{Complex{Float32},1},inter_ratio::Int64;method::Symbol=:linear)
	if inter_ratio==1
		return data_
	end

	if method==:linear
		begin
		counter_::Int64=1
		len_data=length(data_)
		last_::Complex{Float32}=0.0
		dummy_=Array{Complex{Float32},1}(undef,len_data*inter_ratio)
		for i in 1:1:len_data
			dummy_[counter_]=data_[i]
			counter_+=1
			dt_r=(real(data_[i])-real(last_))/inter_ratio
			dt_img=(imag(data_[i])-imag(last_))/inter_ratio
			for j in 1:1:inter_ratio-1
				dummy_[counter_]=Complex( ( real(data_[i]) + (dt_r*j)  ) , ( imag(data_[i]) + (dt_img*j) ) )
				counter_+=1
			end
		last_=data_[i]
		end
		return dummy_
		end

	elseif method==:zeros
        	buffer=Array{Complex{Float32},1}(undef,length(data_)*inter_ratio)
		counter_=1
	        for i in 1:1:length(data_)
			buffer[counter_]=data_[i]
			counter_+=1
			for j in 1:1:inter_ratio-1
				buffer[counter_]=Complex(0.0f0,0.0f0)
				counter_+=1               	
		        end
		end
	        return buffer
	
	else
		@info "Set method to one of the following :linear,:zeros\n\tReturning to the calling function"
		return
	end
end

function interpolate(data_::Array{Complex{Float64},1},inter_ratio::Int64;method::Symbol=:linear)
	if inter_ratio==1
		return data_
	end

	if method==:linear
		counter_::Int64=1
		len_data=length(data_)
		last_::Complex{Float64}=0.0
		dummy_=Array{Complex{Float64},1}(undef,len_data*inter_ratio)
		for i in 1:1:len_data
			dummy_[counter_]=data_[i]
			counter_+=1
			dt_r=(real(data_[i])-real(last_))/inter_ratio
			dt_img=(imag(data_[i])-imag(last_))/inter_ratio
			for j in 1:1:inter_ratio-1
				dummy_[counter_]=Complex( ( real(data_[i]) + (dt_r*j)  ) , ( imag(data_[i]) + (dt_img*j) ) )
				counter_+=1
			end
		last_=data_[i]
		end
		return dummy_

	elseif method==:zeros
        	buffer=Array{Complex{Float64},1}(undef,length(data_)*inter_ratio)
		counter_=1
	        for i in 1:1:length(data_)
			buffer[counter_]=data_[i]
			counter_+=1
			for j in 1:1:inter_ratio-1
				buffer[counter_]=Complex(0.0f0,0.0f0)
				counter_+=1               	
		        end
		end
	        return buffer
	
	else
		@info "Set method to one of the following :linear,:zeros\n\tReturning to the calling function"
		return
	end
end

function resample(data_::Array{Float64,1},out_samp_rate::Int64;inter_method::Symbol=:linear)
	len_data=length(data_)
	if out_samp_rate==len_data
		return data_
	end
	rational_::Rational{Int64}=out_samp_rate//len_data

	inter_r=numerator(rational_)
	decim_r=denominator(rational_)
	return decimate(interpolate(data_,inter_r;method=inter_method),decim_r)
end

function resample(data_::Array{Float32,1},out_samp_rate::Int64;inter_method::Symbol=:linear)
	len_data=length(data_)
	if out_samp_rate==len_data
		return data_
	end
	rational_::Rational{Int64}=out_samp_rate//len_data

	inter_r=numerator(rational_)
	decim_r=denominator(rational_)
	return decimate(interpolate(data_,inter_r;method=inter_method),decim_r)
end

function resample(data_::Array{Complex{Float64},1},out_samp_rate::Int64;inter_method::Symbol=:linear)
	len_data=length(data_)
	if out_samp_rate==len_data
		return data_
	end
	rational_::Rational{Int64}=out_samp_rate//len_data

	inter_r=numerator(rational_)
	decim_r=denominator(rational_)
	return decimate(interpolate(data_,inter_r;method=inter_method),decim_r)
end

function resample(data_::Array{Complex{Float32},1},out_samp_rate::Int64;inter_method::Symbol=:linear)
	len_data=length(data_)
	if out_samp_rate==len_data
		return data_
	end
	rational_::Rational{Int64}=out_samp_rate//len_data

	inter_r=numerator(rational_)
	decim_r=denominator(rational_)
	return decimate(interpolate(data_,inter_r;method=inter_method),decim_r)
end


