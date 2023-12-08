using CondaPkg
CondaPkg.withenv() do
    run(`python src/IDA_Tether_03.py`)
end
