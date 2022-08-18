# ITKTubeTK-CTTomosynthesisRegistration

This repository contains modules and notebooks (written in python) to register 3D pulmonary (lung) vascular networks, generated from CT scans, with the corresponding 2D vessels in tomosynthesis projection images, generated using stationary digital chest tomosynthesis (sDCT).\
A registration verification notebook is also included in this repository.  The verification notebook generates a combined image with transformed 3D CT vessels and the 3D tomosynthesis reconstructed volumes so the user can visually assess the success of the registration.\
**Clinical Significance**: registering pulmonary vessels from a preoperative CT scan with intraoperative Tomosynthesis images will enhance image guided lung biopsies by increasing the practicality of the procedure and reducing the patient’s exposure to high CT radiation [^1].\
[**Research Presentation Slides**](ReadMeFiles/ChestTomoPresentation.pdf)

## Folder Organization
### Data
**Folders**
- *CT_02*: Contains the single CT file
- *TomoAnnotation_02*: Contains 8 example vessel annotation files
- *TomoProjection_02*: Contains 29 tomosynthesis projection images and the geo.txt file
- *TomoReconstruction_02*: Contains 80 tomosynthesis reconstruction slices and one 3D tomosynthesis reconstruction volume
- *Vessels_02*: Contains segmented vessel file
### Results
**Folders**
- *Masks*: Folder where mask images are written
- *MaskedTomo*: Folder where masked tomosynthesis projections are written
- *SolutionTomo*: Folder were solution-transformed vessel overlays are written

**Files**
- combined_02: transformed combined CT and tomosynthesis reconstruction volume image file
- combined_02_untransformed: untransformed CT and tomosynthesis reconstruction volume combined image file
- CT-Lungs-Vessels_02_transformed: transformed vessel file
- CT-Lungs-Vessels_02_untransformed: unstransformed vessel file
- regSolution_02: text of 6 parameter solution to registration 
## Registration
| ![Registration](ReadMeFiles/ReadMe1.png)|
|:--:| 
| *Registration: Find 3D translation and rotation parameters that align 3D vessel perspective projection with 2D tomosynthesis projection image* |
___
**Input**
- *tomoProj_dir*: path to directory of tomosynthesis projection images
- *tomoRecon_dir*: path to directory of tomosynthesis reconstruction images
- *overlay_dir*: path to directory of hand-drawn vessel annotation files
- *emitterGeo_file*: path to file that contains emitter positions
- *vessel_file*: path to .tre file that contains tubes spatial objects of segmented pulmonary vasculature
- *ct_file*: path to CT volume
- *destDir*: path to directory where masked images should be written
- *parentDir*: path to directory where results should be written 
- *solution_output_filename*: string containing .txt name of solution file
- *x_init*: array that contains initial estimated transform to situate vessels prior to registration
- *tomoFileNumber*: number greater than the number of tomosynthesis projections
- *tomoReconFileNumber*: number greater than the number of tomosynthesis reconstruction slices

**Output**:
- *1x6 array*: [xtranslation, ytranslation, ztranslation, zrotation, yrotation, xrotation]

**How To Get .tre Vessel File**
- 	[Kitware Medical: Lung Vessel Segmentation](https://github.com/KitwareMedical/ITKTubeTK-CTLungs) (experiments folder)

**How To Make Annotated Vessel Files**
- [Kitware Medical: ImageViewer](https://github.com/KitwareMedical/ImageViewer) 
- Open image in ImageViewer
- **\\** button until PAINT2D appears in bottom right corner
- **\[** and **\]** buttons to adjust radius size
- Draw over prominant vessels in tomosynthesis reconstruction slices (~8 images)
- **Shift** + **”** buttons to save
- Save file as *“overlayMask_0XX.dcm”* where *XX* is mask number

**View Overlay**
- In command window: ```Imageviewer <path-to-image> -j <path-to-overlay>```

**Expected Ouput**
| ![OutputPlot](ReadMeFiles/ReadMe3.png)|
|:--:| 
| *Plot: Subeset of the un-registered data points (Yellow) and registered data points (Red) plotted* |
| ![OutputOverlay](ReadMeFiles/ReadMe4.png)|
|:--:| 
| *Overlay: The pipieline saves 29 projections of the CT points transformed with the registration solution for each emitter position, visualize with corresponding tomosynthesis image using the "View Overlay command listed above* |
## Verification

**Input**
- *pipelineDir*: path to directory where results should be written 
- *tube_file*: path to .tre file that contains tubes spatial objects of segmented pulmonary vasculature
- *tube_out_file*: path to directory and file name of transformed tube file
- *recon_file_3D*: path to 3D tomosynthesis reconstruction
- *ct_file*: path to CT volume
- *solution_output_filename*: string containing file name to registration solution
- *combined_image_filename*: string containing name of combined image file
- *x_init*: array that contains initial estimated transform to situate vessels prior to registration (from registration, except translation in direction may need to be modified and rotation all 0 because re-orientation is not needed)
- *img_3D*: path to tomosynthesis reconstruction with corrected spacing

**Output**
- Combined images
- Compare output of Registration Verification and Untransformed Registration Verification to see if vessel alignment was successful

---
| ![Verification](ReadMeFiles/ReadMe2.png)|
|:--:| 
| *Verification: Compare untransformed, pre-registration CT vessels with transformed, post-registration CT vessels. See how CT vessels align with visable vessels in tomosynthesis reconstruction* |

[^1]: https://www.researchgate.net/publication/269186336_Stationary_chest_tomosynthesis_using_a_carbon_nanotube_x-ray_source_array_A_feasibility_study
