%% change current path to where the file is located
addpath('../landmark_detect/');
addpath('../matlab_ext_libs/MedicalImageProcessingToolbox/processing/IO');
addpath('../matlab_ext_libs/MedicalImageProcessingToolbox/class_image/');
% 
%% a window was already added to the cropped DRRs (T10-S1)

dir_DRR = '../data/auto_DRR/';
filename = 'DRRsur.mhd'; % the input is an abdominal DRR covers the subregion from Th10 to S1
dir_measure = '../data/auto_measures/';
ref_size = [100; 100]; % if you want all the measurements scaled to as DRR is 100 * 100 dimension

landmarks = Rib_detection(dir_DRR, filename,1,"N",ref_size);% "N" means not scaled/ switch to "Y" if you want to calculate the scaled one
write_to_file(strcat(dir_measure, 'measurement_ref','.txt'),landmarks); % write measurements to a txt file
   

                           




