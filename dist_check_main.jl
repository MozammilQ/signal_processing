
using Distributed
worker_ip="192.168.1.2"
worker_user="root"
num_cores=4
node_working_dir="/home/red/Desktop/Julia_Workspaces/sig_proc"
cd(node_working_dir)
workers_specs=("$worker_user@$worker_ip",num_cores)
machines=Vector{typeof(workers_specs)}(undef,1)
machines[1]=workers_specs
workers_pids=addprocs(machines,max_parallel=num_cores,dir=node_working_dir)


@everywhere include("dist_check_include.jl")

remotecall_wait(iq_file_open,node_file,file_name)

remotecall_wait(update_IQ_data!,node_file)

remotecall_wait(update_fft!,node_fft)

remotecall_wait(get_psd,node_psd,shared_fft)

