using FFTW,Makie

freq_list=[ 5, 10, 20, 40];
amp_list=[1.0, 1.0, 1.0, 1.0];

if length(freq_list) != length(amp_list)
        @info "Amplitude of any frequency is missing"
end
scene=Scene();
len_freq_list=length(freq_list)
samp_rate=maximum(freq_list)*2;
dt=1/samp_rate;
t=range(0,2*pi,step=dt);
y_t = zeros(length(t));

for i=1:len_freq_list
        global y_t = y_t + (amp_list[i] .* sin.( (2*pi) * freq_list[i] * t ));
end

n=length(t);
f_hat_t = fft(y_t);
PSD_t= ( f_hat_t .* conj.(f_hat_t) ) ./ n;
freq_rec_t= range(0,n-1,step=1) .* ( 1 / (dt * n) ) ;

PSD_t[1:1:Int(floor(length(t)/2))]=PSD_t[Int(floor(length(t)/2)):-1:1]

PSD_t[Int(ceil(length(t)/2)):1:end]=PSD_t[end:-1:Int(ceil(length(t)/2))]


#freq_rec_t_half=freq_rec_t[1:Int(floor(n/2))];
#PSD_t_half=PSD_t[1:Int(floor(n/2))];

lines!(scene,freq_rec_t,PSD_t)

