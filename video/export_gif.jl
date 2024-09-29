using Colors, FixedPointNumbers, Images, FFMPEG

folder = "video"
framerate = 20
gifname = joinpath(folder, "Tether.gif")
mp4name = joinpath(folder, "Tether.mp4")

A = Array{Array{RGB{Normed{UInt8,8}},2},1}()
files = readdir(folder)
pngfiles = filter(file->(occursin.("png",file)), files )
rm(gifname, force=true)
rm(mp4name, force=true)

FFMPEG.ffmpeg_exe(`-framerate $(framerate) -f image2 -i $(folder)/img-%4d.png -c:v libx264 -pix_fmt yuv420p -y $(mp4name)`)
FFMPEG.ffmpeg_exe(`-i $(mp4name) -filter_complex "[0]split[a][b]; [a]palettegen[palette]; [b][palette]paletteuse" -y $(gifname)`)

println("\nGif file $(gifname) created!")
println("Mp4 file $(mp4name) created!")