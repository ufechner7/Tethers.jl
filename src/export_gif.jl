using Pkg; Pkg.add("FileIO"); Pkg.add("ImageIO"); Pkg.add("Colors"); Pkg.add("FixedPointNumbers"); Pkg.add("Images")
using FileIO, ImageIO, Colors, FixedPointNumbers, Images

folder = "video"

A = Array{Array{RGB{Normed{UInt8,8}},2},1}()

for j in 0:(600)
    global A
    local img
    img_path = folder*"/"*"img-"*lpad(j,4,"0")*".png"
    img=FileIO.load(img_path)
    if j==0
        A = img
    else
        A = cat(A, img, dims=3)
    end
end

FileIO.save(folder * "/Tether.gif", A)
Pkg.rm("FileIO"); Pkg.rm("ImageIO"); Pkg.rm("Colors"); Pkg.rm("FixedPointNumbers"); Pkg.rm("Images")  
println("\nGif file $(folder * "/Tether.gif") created!")  