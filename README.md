# Surface-Reconstruction-Toolbox
![alt text](/doc/scb_image.PNG)

## Introduction
This tool was mainly coded during my master thesis. It offers several algorithms to triangulate 3d scattered points for 3d scanning applications. After finding it again on a old hard disk I thought that github could give it a fresh restart.
During the thesis, in collaboration with Luca Di angelo and Paolo Di Stefano 3 papers were published

- Luca Di Angelo, Paolo Di Stefano, Luigi Giaccari: A new mesh-growing algorithm for fast surface reconstruction. Computer-Aided Design 43(6): 639-650 (2011).
https://pdfs.semanticscholar.org/4332/35ea0ca6b6db86131ceb206bc160aab65614.pdf
- Luca Di Angelo and Luigi Giaccari, 2011, An Efficient Algorithm for the Nearest Neighborhood Search for Point Cloud, IJCSI International Journal of Computer Science Issues, Vol. 8, Issue 5.
https://www.ijcsi.org/papers/IJCSI-8-5-1-1-11.pdf
- Luca Di Angelo, Paolo Di Stefano, Luigi Giaccari: A Fast Mesh-Growing Algorithm for Manifold Surface Reconstruction. Computer-Aided Design and Applications Pages 197-220 | Published online: 09 Aug 2013
http://www.improve2011.it/Full_Paper/33.pdf

## Downloads (binaries)
coming soon

## Description

The toolbox features at the moment 4 main algorithms:
- **SCBMesher** Fast adv front without Delaunay triangulation support. Can easly generate .5 M trias/secs on a laptop- Based on:
  A new mesh-growing algorithm for fast surface reconstruction. Computer-Aided Design 43(6): 639-650 (2011) 
- **RobustCrust** Based on the Powercrust always return a watertight surface therefore only supports closed shell-like surfaces
- **BallPivoting** Simulates a ball rolling on the surface, absed on the paper from Fausto Bernardini.
- **QuickTerrain** to triangula surface in the z=f(x,y) form, usually called "terrain".
 
Tools are designed to be used both interactively and from command prompt:

- scbmesher.exe -in bunny.dat –out bunny.stl
- robustcrust.exe –in bunny.dat –out bunny.stl
- ballpivoting.exe –in inputfile –out outputfile -r 1.0
- quickterrain.exe –in inputfile –out output file

All tools uses the -in and -out switch for I/O. The Ball Pivoting also asks for the ball radius.

## Input/Output
### Input file formats
- **.dat (binary)**
The file starts with an “int” typeindicating the number of points, next data are a list of doubles indicating the coordinates.
Points are stored this way: X1,Y1,Z1,X2,Y2,Z2…….
These are the Matlab commands to generate such file:
```Matlab
fid=fopen(name,'wb');%open input file
fwrite(fid, np, 'int');%first data=the numbers of points;
fwrite(fid,p,'double');%points;
fclose(fid);%close the file
```
And this are the C++ ones:
```C++
pFile =fopen("Input.dat", "wb");//open input file
nwritten=fwrite(&N, sizeof(int), 1, pFile);//first line (the number of points)
nwritten=fwrite(&p, sizeof(double), N*3, pFile);//Write the points coordinate
fclose(pFile);//close the file
```
- **.cgo ( ASCII)**
The first line contains the number of points, the followings the list of 3D coordinates. Coordinates are separted by space, the decimal separator can be a comma or a point.
Example:
```Ascii
6
0,060649 0,011160 0,059552
0,060351 0,012167 0,058159
0,060749 0,012164 0,059571
0,060267 0,011162 0,058137
0,060278 0,013168 0,058137
0,060653 0,013164 0,059557
```

### Output File 
- **Stl(binary)**
Standard binary stl file.


## Acknowledgements
Many thanks to:
- We used the robust predicates from  Jonathan Richard Shewchuk
- Amenta's directory of computational geometry software. We found great ideas in the Power Crust Paper. 
- Fausto Bernardini. Our Ball Pivoting implementation is based on his paper.
- Thanks to the Meshlab team for their software
- tetgen http://wias-berlin.de/software/tetgen/ saved a lot of time

## Contacts
For bugs,suggestions,questions refer to; giaccariluigi@msn.com

# Development
Build

Astyle formatting
./Astyle.exe  --options=astyle_options.txt ./*.cpp  ./*.h ./*.c

