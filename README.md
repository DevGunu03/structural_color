# Structural Color
## Aim
The code is created specifically to determine the number of layers of Polystyrene beads based photonic crystal array in nano-dimensions. Although the Transfer Matrix Model includes generation of color visuals for a wide variety of materials, the latter part of the code which deploys CIECAM-UCS based TCD for layer identification is linked to PS-beads.

## Using the codes - one step at a time
### 01. TRANSFER MATRIX MODEL
The first step in the beginning of the identification is generation of a reference model. I used TMM package created by George F. Burkhard and Eric T. Hoke, whose original source can be found at [McGehee Group](https://web.stanford.edu/group/mcgehee/transfermatrix/). I modified the code to be usable in our case, by adding the Effective Medium Approximation. See, the real TMM model is a basic tool which simulates a scattering event between a light source and a physical system containing an arbitrary number of cuboidal layers - no other shape is allowed. This is even more limiting since only the thickness of the layer is important and not even the actual x,y-dimensions. To simulate results of photonic crystal arrays, we need a model which can introduce spacings in the layer, such that not all space in the thickness we mention is taken up by the material, but some by air, since it is a void. This is taken into consideration by the Maxwell-Garnett Equation which looks like this:\
<p align="center">
<img width="500" height="214" alt="01ema" src="https://github.com/user-attachments/assets/8f9054cb-04b9-43a7-81a5-36c4e65112e5">  
</p>

This equation can simulate a matrix and inclusion, given a packing fraction, which we can calculate knowing the PS-beads configuration is HCP, as can be seen from the SEM images.\
<p align="center">
<img width="414" height="287" alt="02hcp" src="https://github.com/user-attachments/assets/35f8b96c-4f14-4728-adcc-69a7a845e6c6" />  
</p>

Next, I'll take you through the steps in which the code shall be used to create the reference dataset:
The original code has been divided into two codes - TransferMatrix_multiple and TransferMatrix_packing. The former should be used whenever **variable thickness** has to be simulated - thus, let's say a Air/PS/SiO2/Si system, with thickness for PS varying from 0nm to 1500nm with step size of 300nm. The user can also update thickness of multiple layer too. The latter should be used whenever **variable packing-fraction** has to be simulated - thus, for the same system as mentioned above, the *packing-fraction* of PS varies over a range, while the thickness is contant.
<p align="center">
<img width="477" height="194" alt="03TMM" src="https://github.com/user-attachments/assets/27170686-cd1d-4c75-ba9e-33edc3f6c7f1" />  
</p>

In the above image, the user has the ability to modify the initial thickness, the step-size and the final thickness of the PS beads to be simulated. The wavelength range is already mentioned, but can be modified - the ones currently written is the advisable range for converting a reflectance spectra into sRGB format. 
The layers is where the system has to be established. Starting with 'Air', whose thickness can be arbitrary, since it is not used for calculation but whose presence is essential for the code to run. The next is the layer which comes in contact with the light, followed by the rest in order. The names of the layers have to consulted from the file "Index_of_Refraction_library.xls", from which the complex refractive indices are actually taken for computation. Followed next is the thickness of each layer, in the same order in which the names have been written in the layers list. **incl** is a variable which dictates the packing fraction. The main function is hardcoded to only use EMA or packing fraction information only when **'PS-beads'** is mentioned. This could be changed in the main file TransferMatrix_multiple (line 93) or TransferMatrix_packing (line 98) as follows:
<p align="center">
<img width="703" height="143" alt="04emaInTMM" src="https://github.com/user-attachments/assets/06aa7a16-da6a-4a9c-9b0d-6cfff9284f79" />
</p>

Add more materials along with 'PS-beads' as per the interest.  

<p align="center">
<img width="634" height="261" alt="05TMM" src="https://github.com/user-attachments/assets/ef26a59b-a21c-4be1-b889-b639e3f9b712" />
</p>
This is where the physics and linear algebra takes over and scattering problem is solved using TMM formalism. The files are stored with filename as suggested by the variable with the same name. The user can keep *incl=1* if they don't want to use EMA. The main function called in the image above is TransferMatrix_multiple. This can be changed to TransferMatrix_packing and multiple *incl* values can be uploaded for each layer in order from near to light source to farther. Similar changes has to be done to the main TransferMatrix_packing code, to accomodate as much number of packing fraction for each layer.

If the user wants to have multiple packing-fraction for the same layer, a *for-loop* can be created which changes the *pf* every time a file is stored.

### 02. The colorbar
The next two codes titled *2_makes_colorbar_packing.ipynb* and *3_makes_colorbar_thickness.ipynb* is for creating colorbar, with the x-ticks showing packing fraction ranging from 0% to 100% for the former file and range of thickness for the latter file. Both the files in their first and second steps opens up a GUI for the user to select the folder containing a file name with a type as mentioned in step 3 below:\
<p align="center">
<img width="822" height="249" alt="image" src="https://github.com/user-attachments/assets/742ebbf1-d780-49a9-8ad7-164d19457f46" />
</p>

The filename should be of the type **PS_300nm4** or **PS_300nm** for former and latter files. The first number *300nm* suggests the thickness while the second number for *2_makes_colorbar_packing.ipynb* suggests packing, in percentage. IN step 5 for both the files, the user can modify the ticks for each files to be shown and how the final image generated should look. Both the files eventually store the images in *.png* format in the folder where the code resides.
<p align="center">
  <img width="1200" height="138" alt="image09" src="https://github.com/user-attachments/assets/5d76e175-45a5-4536-95de-2f04a20171c1" />
  <img width="1152" height="134" alt="image05" src="https://github.com/user-attachments/assets/1ce3c42e-d4e7-4336-9583-f021371489df" />
</p>
