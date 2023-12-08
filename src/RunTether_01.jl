using CondaPkg
CondaPkg.withenv() do
    run(`python src/Tether_01.py`)
end
