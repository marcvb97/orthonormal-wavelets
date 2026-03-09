# Orthonormal-wavelets
Image compression using orthonormal wavelets is compared to image compression using the interpolating de la Vallée Poussin wavelets, Daubechies 2 and biorthogonal 3.5 wavelets.

The datasets with images that are used:
- KODAK - 24 - [768x512] - [512x768] - https://r0k.us/graphics/kodak/
- NY17 - 17 - [1200×1600] – [5430×3520] - https://www.gcc.tu-darmstadt.de/home/proj/dpid/index.en.jsp
- NY96 - 96 - [500×334] – [6394×3456] - https://www.gcc.tu-darmstadt.de/home/proj/dpid/index.en.jsp
- 13US - 13 - [241×400] – [400×310] - https://www.cl.cam.ac.uk/~aco41/Files/Sig15UserStudyImages.html
- URBAN100 - 100 - [1024×564] – [1024×1024] - https://paperswithcode.com/dataset/urban100
- PEXELS300 - 300 - [300×300] - https://github.com/ImgScaling/LCIscaling

Some of the Matlab files:
- compression_of_images_of_a_folder.m: main script to perform the compression
- analyze_results.m: main script to analyze the different .mat files
