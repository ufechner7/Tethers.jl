using Colors, FixedPointNumbers, Images, FFMPEG

folder = "video"
framerate = 20
gifname = joinpath(folder, "Tether.gif")
mp4name = joinpath(folder, "Tether.mp4")

A = Array{Array{RGB{Normed{UInt8,8}},2},1}()
files = readdir(folder)
pngfiles = filter(file->(occursin.("png",file)), files )

FFMPEG.ffmpeg_exe(`-framerate $(framerate) -f image2 -i $(folder)/img-%4d.png -vf palettegen -y $(gifname)`)
FFMPEG.ffmpeg_exe(`-framerate $(framerate) -f image2 -i $(folder)/img-%4d.png -c:v libx264 -pix_fmt yuv420p -y $(mp4name)`)

println("\nGif file $(gifname) created!")
println("Mp4 file $(mp4name) created!")