using FileIO, ImageIO, Colors, FixedPointNumbers, Images

folder = "video"

A = Array{Array{RGB{Normed{UInt8,8}},2},1}()
files = readdir(folder)
pngfiles = filter(file->(occursin.("png",file)), files )

for j in 0:(length(pngfiles)-1)
    global A
    local img
    img_path = folder*"/"*"img-"*lpad(j,4,"0")*".png"
    img=FileIO.load(img_path)
    if j==0
        A = img
    else
        A = cat(A, img, dims=3)
    end
    if j%10 == 0
        println("Processing image $img_path of $(length(pngfiles)-1)")
    end
end

FileIO.save(folder * "/Tether.gif", A)
println("\nGif file $(folder * "/Tether.gif") created!")  