using PythonCall
using Logging

const _py_logging = Ref{Py}()
function py_logging()
    if !isassigned(_py_logging)
        _py_logging[] = pyimport("logging")
    end
    return _py_logging[]
end

const _min_level = Ref{LogLevel}(Info)

struct PythonLogBridge <: AbstractLogger end
Logging.shouldlog(::PythonLogBridge, level, _module, group, id) = true
Logging.min_enabled_level(::PythonLogBridge) = _min_level[]
Logging.catch_exceptions(::PythonLogBridge) = true

function Logging.handle_message(
    ::PythonLogBridge, level, message, _module, group, id, filepath, line;
    kwargs...
)
    log = py_logging()
    name = _module !== nothing ? string(_module) : "julia"
    pylevel = if level == Debug
        10
    elseif level == Info
        20
    elseif level == Warn
        30
    elseif level == Error
        40
    else
        50
    end
    extra = pydict(Dict(string(k) => string(v) for (k, v) in kwargs))
    log.getLogger(name).log(pylevel, string(message), extra=extra)
end

function use_python_logging!()
    log = py_logging()
    pylevel = pyconvert(Int, log.root.level)
    _min_level[] = if pylevel <= 10
        Debug
    elseif pylevel <= 20
        Info
    elseif pylevel <= 30
        Warn
    else
        Error
    end
    global_logger(PythonLogBridge())
end