# image-processing-for-psm

## psm.py is the python code the do an aberration analysis with PSM images

### Requirements: 
  -open cv,
  -numpy,
  -matplotlib,
  -os
  
 -The only input of the function is a folder with the PSM images with the following format:

 ![exemple_dossier](https://user-images.githubusercontent.com/73662195/177781194-564e7500-2d2d-4070-b502-2fc7e5c83a63.JPG)

-To use download psm.py and change the argument on line 395 to your folder name. Alternatively, create a new python file in
the same folder as psm.py a use the following code:

from psm import PSM

PSM.psm('your_folder_name')

## Other files description

### All images are duplicates of pictures present in the PDFs

### psm.pdf


