% Batch run the RunRegistration script to register all images for each patient

%% PATIENT 1: AIDR3D STRONG FC03 DATA
% We DO NOT register CFA data
% convert mat file to dicom
% Write DICOM Segmentation

thresh_background = -1024; % WE SHOULD NEVER HAVE A PERFUSION VALUE EVEN CLOSE TO THIS
voxel_scale_factor = 100;

% Load mask
im_mask_rest = loadMAT('\\polaris.radsci.uci.edu\Data5_6\HUMAN_PERFUSION_DATA\05_11_18_data_py_auto\PATIENT_5\AIDR_3D_STR_FC03\REST_FULL\ASSIGN_DCM\myocardium_perfusion_bed_FULL.mat');
im_mask_rest = double(logical(im_mask_rest));

% Load rest perfusion
rest_map = loadMAT('\\polaris.radsci.uci.edu\Data5_6\HUMAN_PERFUSION_DATA\05_11_18_data_py_auto\PATIENT_5\AIDR_3D_STR_FC03\REST_FULL\perf_map.mat');
outputPath = '\\polaris.radsci.uci.edu\Data4_New\Seth\coronary_flow_capacity\patient_data\PATIENT 5\DICOM\PERF_DCM_01';
dcmPath = '\\polaris.radsci.uci.edu\Data5_6\HUMAN_PERFUSION_DATA\05_11_18_data_py_auto\PATIENT_5\AIDR_3D_STR_FC03\REST_FULL\DICOM' ;
dcmTitle = 'rest_perf_map_dcm';

% Threshold the perfusion map - Set background voxels to a VERY LOW number
rest_map_thresh = rest_map;
rest_map_thresh(im_mask_rest == 0) = thresh_background;
rest_map_thresh(im_mask_rest == 1) = rest_map_thresh(im_mask_rest == 1) * voxel_scale_factor;
try
    WriteDICOM( rest_map_thresh, dcmPath, outputPath, dcmTitle );
catch
    warning('Issue while writing rest DICOMs')
end


% Load mask
im_mask_stress = loadMAT('\\polaris.radsci.uci.edu\Data5_6\HUMAN_PERFUSION_DATA\05_11_18_data_py_auto\PATIENT_5\AIDR_3D_STR_FC03\STRESS_FULL\ASSIGN_DCM\myocardium_perfusion_bed_FULL.mat');
im_mask_stress = double(logical(im_mask_stress));

% Load stress perfusion
stress_map = loadMAT('\\polaris.radsci.uci.edu\Data5_6\HUMAN_PERFUSION_DATA\05_11_18_data_py_auto\PATIENT_5\AIDR_3D_STR_FC03\STRESS_FULL\perf_map.mat');
outputPath = '\\polaris.radsci.uci.edu\Data4_New\Seth\coronary_flow_capacity\patient_data\PATIENT 5\DICOM\PERF_DCM_02';
dcmPath = '\\polaris.radsci.uci.edu\Data5_6\HUMAN_PERFUSION_DATA\05_11_18_data_py_auto\PATIENT_5\AIDR_3D_STR_FC03\STRESS_FULL\DICOM' ;
dcmTitle = 'stress_perf_map_dcm';

stress_map_thresh = stress_map;
stress_map_thresh(im_mask_stress == 0) = thresh_background;
stress_map_thresh(im_mask_stress == 1) = stress_map_thresh(im_mask_stress == 1) * voxel_scale_factor;
try
    WriteDICOM( stress_map_thresh, dcmPath, outputPath, dcmTitle );
catch
    warning('Issue while writing stress DICOMs')
end



%% Register the Stress and Rest images
BB = [1, 1, 1, 510, 510, 230]; 
study = '\\POLARIS\Data 1_2_4\CT_Perfusion\WLA_VA_DATA\Subject005\';
acq_dirs = {...
     [study 'TEMP'],...
    };

% % Define master directory to register all others against.
master_dir = [study 'TEMP']; 
reference_fname = 'nii1'; % full scan V2

for idx = 1 : numel(acq_dirs)
    output_text = [ 'Running registration on: ' acq_dirs{idx} ];
    disp( output_text );
    
    acq_dir = acq_dirs{idx};
    RunRegistrationDICOM_mod(acq_dir, master_dir, reference_fname, BB)
    
end

%% PATIENT 1:  Calculate CFR
path_study = '\\polaris.radsci.uci.edu\Data4_New\Seth\coronary_flow_capacity\patient_data\PATIENT 5\REGISTERED\';
path_im_rest = [path_study  'mat\01.mat'];
path_im_stress = [path_study 'mat\02.mat'];


im_rest = loadMAT(path_im_rest);
im_stress = loadMAT(path_im_stress);

im_stress = (im_stress - 15000);
im_rest = (im_rest - 15000);

im_rest(im_rest <= thresh_background*0.80) = thresh_background;
im_stress(im_stress <= thresh_background*0.80) = thresh_background;

im_rest(im_rest > thresh_background) = im_rest(im_rest > thresh_background) / voxel_scale_factor;
im_stress(im_stress > thresh_background) = im_stress(im_stress > thresh_background) / voxel_scale_factor;

%% Crop images
%min_perf_cutoff = -10; % We noted that most perfusion values should be greater than -7

% Define bounding box for stress and rest images - BB is the same for
% rest/stress since both images are registered

BB1 = [30, 30, 30, 455, 455, 190];
cropped_im_rest = imcrop3D(im_rest, BB1, thresh_background);
cropped_im_stress = imcrop3D(im_stress, BB1,thresh_background);


%% Create myocardial mask
im_mask_registered = ones(size(cropped_im_rest));
im_mask_registered((cropped_im_stress <= min_perf_cutoff) | (cropped_im_rest <= min_perf_cutoff)) = 0;

%% Calculate CFR

% Clean up erroneous voxel values - it likely occurs due to small
% misregistration errors, e.g. we are dividing by zero in some cases,
% giving us Inf or "Not A Number" values.
im_cfr = cropped_im_stress ./ cropped_im_rest;
im_cfr(im_mask_registered ~= 1) = 0;

im_cfr(isnan(im_cfr)) = 0;
im_cfr((im_cfr == Inf) | (im_cfr == -Inf)) = 0;

%% Calculate CFC Index
% https://www.sciencedirect.com/science/article/pii/S1936879815012911?via%3Dihub
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5489577/
% https://www.sciencedirect.com/science/article/pii/S1936878X12000952?via%3Dihub

cfc_idx_1 = ((im_cfr <= 1)) &...
            ((im_stress <= 0.91)) &...
            (im_mask_registered == 1);

cfc_idx_2 = ((im_cfr > 1) & (im_cfr <= 1.74)) &...
            ((im_stress > 0) & (im_stress <= 0.91)) &...
            (im_mask_registered == 1);
        
cfc_idx_3 = ((im_cfr > 1.74) & (im_cfr <= 2.03)) &...
            ((im_stress > 0.91) & (im_stress <= 1.12)) &...
            (im_mask_registered == 1);
        
cfc_idx_4 = ((im_cfr > 2.03) & (im_cfr <= 2.70)) &...
            ((im_stress > 1.12) & (im_stress <= 1.76)) &...
            (im_mask_registered == 1);
        
cfc_idx_5 = ((im_cfr > 2.70) & (im_cfr <= 3.37)) &...
            ((im_stress > 1.76) & (im_stress <= 2.39)) &...
            (im_mask_registered == 1);
        
cfc_idx_6 = ((im_cfr > 3.37)) &...
            ((im_stress > 2.39)) &...
            (im_mask_registered == 1);

% Create CFC indices - worse CFC is lower in number
im_cfc = zeros(size(im_mask_registered));
im_cfc(cfc_idx_1 == 1) = 1;
im_cfc(cfc_idx_2 == 1) = 2;
im_cfc(cfc_idx_3 == 1) = 3;
im_cfc(cfc_idx_4 == 1) = 4;
im_cfc(cfc_idx_5 == 1) = 5;
im_cfc(cfc_idx_6 == 1) = 6;



%% Image Data Visualization

% NOTE - LARGE CLUSTER OF 0 VOXELS IS DUE TO POOR RV SEGMENTATION IN THE
% UNREGISTERED REST MASK!
f_handle = figure(); % f_handle can be used to change things like axis limits
title('REST - Normalized Count');
hold on;
histogram(cropped_im_rest(im_mask_registered == 1),'Normalization', 'probability','BinWidth', 0.3);
hold on;
histogram(rest_map(im_mask_rest == 1),'Normalization', 'probability', 'BinWidth', 0.3);
hold on;
legend('Registered', 'Unregistered');
xlabel('Perfusion (g/mL/min)');
ylabel('Fraction of Voxels (# in bin / total # of voxels)');
xlim([-10, 10]);
ylim([0, 0.20]);

%% Image Data Visualization
f_handle = figure(); % f_handle can be used to change things like axis limits
title('STRESS - Normalized Count');
hold on;
histogram(cropped_im_stress(im_mask_registered == 1),'Normalization', 'probability','BinWidth', 0.3);
hold on;
histogram(stress_map(im_mask_stress == 1),'Normalization', 'probability', 'BinWidth', 0.3);
hold on;
legend('Registered', 'Unregistered');
xlabel('Perfusion (g/mL/min)');
ylabel('Fraction of Voxels (# in bin / total # of voxels)');
xlim([-10, 10]);
ylim([0, 0.20]);

%% Image Data Visualization
f_handle = figure(); % f_handle can be used to change things like axis limits
title('STRESS vs. REST - Normalized Count');
hold on;
histogram(cropped_im_stress(im_mask_registered == 1),'Normalization', 'probability','BinWidth', 0.3);
hold on;
histogram(cropped_im_rest(im_mask_registered == 1),'Normalization', 'probability', 'BinWidth', 0.3);
hold on;
legend('STRESS', 'REST');
xlabel('Perfusion (g/mL/min)');
ylabel('Fraction of Voxels (# in bin / total # of voxels)');
xlim([-10, 10]);
ylim([0, 0.20]);

%% Image Data Visualization
f_handle = figure(); % f_handle can be used to change things like axis limits
title('CFR - Normalized Count');
hold on;
histogram(im_cfr(im_mask_registered == 1),'Normalization', 'probability','BinWidth', 0.3);
xlim([-15, 15]);

%% Image Data Visualization
f_handle = figure(); % f_handle can be used to change things like axis limits
title('CFC - Normalized Count');
hold on;
histogram(im_cfc(im_cfc > 0),'Normalization', 'probability','BinWidth', 1);
xlim([1, 6]);

%% Image Data Visualization
f_handle = figure();
title('Stress vs CFR - Perfusion');
hold on;
X = cropped_im_stress(im_mask_registered == 1);
Y = im_cfr(im_mask_registered == 1);

data = [X, Y];
hh3 = hist3(data, 'NBins', [10, 10]);
colorbar
image(flipud(hh3));
xlabel('STRESS');
ylabel('CFR');


%% Scatter plot
figure();
title('CFR vs Stress Perfusion');
hold on;
scatter(X, Y);
xlabel('STRESS (mL/min/g)');
ylabel('CFR');




%% CALCUALTE MEAN OF STRESS / REST PERFUSION BEFORE AND AFTER REGISTRATION - MAKE SURE THERE WAS NO SIGNIFICANT CHANGE IN PERFUSION MEASUREMENTS
% This needs to be refined.  im_mask is not a perfect segmentation of the
% left ventricle myocardium.  This is causing our calculation to be
% erroneous - we are include additional voxels in our means of the
% registered data.  Consider redoing the segmentations using the registered
% data.

%mean
perf_rest_original = mean(rest_map(im_mask_rest ~= 0));
perf_stress_original =  mean(stress_map(im_mask_stress ~= 0));

perf_rest_registered = mean(cropped_im_rest(im_mask_registered ~= 0));
perf_stress_registered = mean(cropped_im_stress(im_mask_registered ~= 0));

perf_cfr = mean(im_cfr(im_mask_registered ~= 0));
perf_cfc = mean(im_cfc(im_mask_registered ~= 0));

%stdev
perf_rest_original_stdev = std(rest_map(im_mask_rest ~= 0));
perf_stress_original_stdev = std(stress_map(im_mask_stress ~= 0));

perf_rest_registered_stdev = std(cropped_im_rest(im_mask_registered ~= 0));
perf_stress_registered_stdev = std(cropped_im_stress(im_mask_registered ~= 0));

perf_cfr_stdev = std(im_cfr(im_mask_registered ~= 0));
perf_cfc_stdev = std(im_cfc(im_mask_registered ~= 0));

%MIN

perf_rest_original_min = min(rest_map(im_mask_rest ~= 0));
perf_stress_original_min =  min(stress_map(im_mask_stress ~= 0));

perf_rest_registered_min = min(cropped_im_rest(im_mask_registered ~= 0));
perf_stress_registered_min = min(cropped_im_stress(im_mask_registered ~= 0));

perf_cfr_min = min(im_cfr(im_mask_registered ~= 0));
perf_cfc_min = min(im_cfc(im_mask_registered ~= 0));

%MAX

perf_rest_original_max = max(rest_map(im_mask_rest ~= 0));
perf_stress_original_max = max(stress_map(im_mask_stress ~= 0));

perf_rest_registered_max = max(cropped_im_rest(im_mask_registered ~= 0));
perf_stress_registered_max = max(cropped_im_stress(im_mask_registered ~= 0));

perf_cfr_max = max(im_cfr(im_mask_registered ~= 0));
perf_cfc_max = max(im_cfc(im_mask_registered ~= 0));

% INCLUDE CFC AS WELL

% NOTE:  Add STDEV, MIN, and MAX values here too; it'd be nice to put this
% in a table

%% SCREENSHOT SLICE 150 OF REST_MAP TO SHOW WHERE THE ZERO VOXELS ARE COMING FROM ON THE HISTOGRAM
rest_map_masked = rest_map;
rest_map_masked(im_mask_rest == 0) = -1024;
imtool3D2(rest_map_masked);