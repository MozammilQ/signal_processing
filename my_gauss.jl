

g_stop=1000
g_step=0.1
shift_right=0

sigma_=3
meu_=500

x=shift_right:g_step:g_stop

y=Array{Float64,1}(undef,length(x))
y.=0.0

for i in 1:1:length(y)
	y[i]= (  1  / ( sigma_*sqrt(2pi) ) ) * exp(  (-0.5)  *  (  ( (x[i]-meu_)  /  sigma_)^2  )  )
end

lines(x,y)


