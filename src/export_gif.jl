# Before running this script you must install the packages that are used
# with `using Pkg; Pkg.add("FileIO")` etc.
using FileIO, ImageIO, Colors, FixedPointNumbers, Images

folder = "video"

A = Array{Array{RGB{Normed{UInt8,8}},2},1}()

for j in 0:146
    global A
    local img
    img_path = folder*"/"*"img-"*lpad(j,4,"0")*".png"
    img=load(img_path)
    if j==0
        A = img
    else
        A = cat(A, img, dims=3)
    end
end

save(folder * "/Tether.gif", A)