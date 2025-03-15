using CondaPkg
CondaPkg.withenv() do
    run(`python examples_python/Radau_Tether_04.py`)
end
