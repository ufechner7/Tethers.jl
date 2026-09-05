using CondaPkg
CondaPkg.withenv() do
    run(`python src/Tether_02.py`)
end
