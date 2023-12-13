using CondaPkg
CondaPkg.withenv() do
    run(`python src/Radau_Tether_04.py`)
end
