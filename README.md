# TFC-Algo-implementation

8108 VITE's regular top probe can't correctly detect the thickness of thin film polyimide layer. So there exists a thin film probe to measure thickness of this layer more accurately.

Inputs:
n_poly = 1.7 // Refractive index of thin film
t=20000 // Thickness of polyimide [nanometers] %physical not optical
M=512 // spectrometer resolution %for VITE, its 512

Output:
correction_thickness // to be added to total thickness of VITE 8108

TFCAlgo.cpp is the main algorithm implementation file.
