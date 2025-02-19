println("Running init...")

using ControlPlots
using REPL.TerminalMenus

options = ["include(\"src/Tether_01.jl\")",
           "include(\"src/Tether_02.jl\")",
           "include(\"src/Tether_03.jl\")",
           "include(\"src/Tether_03b.jl\")",
           "include(\"src/Tether_03c.jl\")",
           "include(\"src/Tether_04.jl\")",
           "include(\"src/Tether_05.jl\")",
           "include(\"src/Tether_06.jl\")",
           "include(\"src/Tether_06b.jl\")",
           "include(\"src/Tether_06c.jl\")",
           "include(\"src/Tether_07.jl\")",
           "include(\"src/Tether_08.jl\")",
           "include(\"src/Tether_09.jl\")",
           "quit()"]

function menu()
    active = true
    while active
        menu = RadioMenu(options, pagesize=8)
        choice = request("\nChoose function to execute or `q` to quit: ", menu)

        if choice != -1 && choice != length(options)
            eval(Meta.parse(options[choice]))
        else
            println("Left menu. Press <ctrl><d> to quit Julia!")
            active = false
        end
    end
end