"""
Logging memory usage using different methods

# Example

* Include mem_logger
```include("mem_logger.jl")```

* Create a global JuliaMemLogger struct
```mem_logger = JuliaMemLogger()```

* Log the memory at different places in the code
```
function step_jules(...)
	...
	log_mem_usage(mem_logger)
	...
end
```
* Print results
```print_logs(mem_logger)```
"""
mutable struct JuliaMemLogger
	mem_usage
	gc_live_bytes
	jit_total_bytes
	maxrss

	function JuliaMemLogger()	
		return new([], [], [], [])
	end
end

function log_mem_usage(log::JuliaMemLogger)
	ret = read(`tasklist /fi "PID eq $(getpid())" /fo csv  /nh`, String)
	mem = filter.(isdigit, (split(ret, ",")[end]) )
	mem = parse(Int64, mem)
	push!(log.mem_usage, mem)
	push!(log.gc_live_bytes, Base.gc_live_bytes()/2^20)
	push!(log.jit_total_bytes, Base.jit_total_bytes()/2^20)
	push!(log.maxrss, Sys.maxrss()/2^20)	
end

function print_logs(log::JuliaMemLogger)
	println(log.mem_usage)
	println(log.gc_live_bytes)
	println(log.jit_total_bytes)
	println(log.maxrss)
end