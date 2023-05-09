%% Simple Approach

% import patient rest

rest_map1 = loadMAT('\\polaris.radsci.uci.edu\Data4_New\Seth\coronary_flow_capacity\patient_data\SUBJECT002\REST_MAP_FILTERED\20211106204819_MEDIAN_FILTER_5_MM.mat');
rest_map2 = loadMAT('\\polaris.radsci.uci.edu\Data4_New\Seth\coronary_flow_capacity\patient_data\SUBJECT004\REST_MAP_FILTERED\20211108084611_MEDIAN_FILTER_5_MM.mat');
rest_map3 = loadMAT('\\polaris.radsci.uci.edu\Data4_New\Seth\coronary_flow_capacity\patient_data\SUBJECT006\REST_MAP_FILTERED\20211108111616_MEDIAN_FILTER_5_MM.mat');
rest_map4 = loadMAT('\\polaris.radsci.uci.edu\Data4_New\Seth\coronary_flow_capacity\patient_data\SUBJECT007\REST_MAP_FILTERED\20211108120525_MEDIAN_FILTER_5_MM.mat');
rest_map5 = loadMAT('\\polaris.radsci.uci.edu\Data5_6\HUMAN_PERFUSION_DATA\05_11_18_data_py_auto\PATIENT_9\AIDR_3D_STR_FC03\REST_FULL\perf_map.mat');


% calculate mean rest

rest_mean1 = mean(rest_map1(rest_map1 ~= 0));
rest_mean2 = mean(rest_map2(rest_map2 ~= 0));
rest_mean3 = mean(rest_map3(rest_map3 ~= 0));
rest_mean4 = mean(rest_map4(rest_map4 ~= 0));
rest_mean5 = mean(rest_map5(rest_map5 ~= 0));


% import patient stress

stress_map1 = loadMAT('\\polaris.radsci.uci.edu\Data4_New\Seth\coronary_flow_capacity\patient_data\SUBJECT002\STRESS_MAP_FILTERED\20211107195152_MEDIAN_FILTER_5_MM.mat');
im_stress2 = loadMAT('\\polaris.radsci.uci.edu\Data4_New\Seth\coronary_flow_capacity\patient_data\SUBJECT004\STRESS_MAP_FILTERED\20220302203429_MEDIAN_FILTER_5_MM.mat');
im_stress3 = loadMAT('\\polaris.radsci.uci.edu\Data4_New\Seth\coronary_flow_capacity\patient_data\SUBJECT006\STRESS_MAP_FILTERED\20211108115142_MEDIAN_FILTER_5_MM.mat');
im_stress = loadMAT('\\polaris.radsci.uci.edu\Data4_New\Seth\coronary_flow_capacity\patient_data\SUBJECT007\STRESS_MAP_FILTERED\20211108121501_MEDIAN_FILTER_5_MM.mat');
stress_map5 = loadMAT('\\polaris.radsci.uci.edu\Data5_6\HUMAN_PERFUSION_DATA\05_11_18_data_py_auto\PATIENT_9\AIDR_3D_STR_FC03\STRESS_FULL\perf_map.mat');

% divide stress by rest mean

cfr1 = stress_map1 ./ rest_mean1;
im_cfr2 = im_stress2 ./ rest_mean2;
im_cfr3 = im_stress3 ./ rest_mean3;
im_cfr = im_stress ./ rest_mean4;
cfr5 = stress_map5 ./ rest_mean5;

%% Averages

cfr1_mean = mean(cfr1(cfr1 ~= 0));
cfr2_mean = mean(cfr2(cfr2 ~= 0));
cfr3_mean = mean(cfr3(cfr3 ~= 0));
cfr4_mean = mean(cfr4(cfr4 ~= 0));
cfr5_mean = mean(cfr5(cfr5 ~= 0));

stress_mean1 = mean(stress_map1(stress_map1 ~= 0));
stress_mean2 = mean(stress_map2(stress_map2 ~= 0));
stress_mean3 = mean(stress_map3(stress_map3 ~= 0));
stress_mean4 = mean(stress_map4(stress_map4 ~= 0));
stress_mean5 = mean(stress_map5(stress_map5 ~= 0));

%% Create Mask

im_mask = double(logical(stress_map1));

%% Calculate CFC

cfc_idx_1 = ((im_cfr <= 1)) &...
            ((im_stress <= 0.91)) &...
            (im_mask == 1);

cfc_idx_2 = ((im_cfr > 1) & (im_cfr <= 1.74)) &...
            ((im_stress > 0) & (im_stress <= 0.91)) &...
            (im_mask == 1);
        
cfc_idx_3 = ((im_cfr > 1.74) & (im_cfr <= 2.03)) &...
            ((im_stress > 0.91) & (im_stress <= 1.12)) &...
            (im_mask == 1);
        
cfc_idx_4 = ((im_cfr > 2.03) & (im_cfr <= 2.70)) &...
            ((im_stress > 1.12) & (im_stress <= 1.76)) &...
            (im_mask == 1);
        
cfc_idx_5 = ((im_cfr > 2.70) & (im_cfr <= 3.37)) &...
            ((im_stress > 1.76) & (im_stress <= 2.39)) &...
            (im_mask == 1);
        
cfc_idx_6 = ((im_cfr > 3.37)) &...
            ((im_stress > 2.39)) &...
            (im_mask == 1);

% Create CFC indices - worse CFC is lower in number
im_cfc = zeros(size(im_mask));
im_cfc(cfc_idx_1 == 1) = 1;
im_cfc(cfc_idx_2 == 1) = 2;
im_cfc(cfc_idx_3 == 1) = 3;
im_cfc(cfc_idx_4 == 1) = 4;
im_cfc(cfc_idx_5 == 1) = 5;
im_cfc(cfc_idx_6 == 1) = 6;
%% Image Data Visualization

save('im_cfc');
im_cfc = load('im_cfc.mat');
im_vol = loadMAT('\\POLARIS\Data4_New\Seth\coronary_flow_capacity\patient_data\SUBJECT007\im_cfc.mat', 7);
%im_vol = im_vol.im_vol;
im_vol = int16(im_vol);
im_vol = im_vol.*10; % Rescale by 10x to improve dynamic range of color in Vitrea
dcmpathIN = '\\POLARIS\Data4_New\Seth\coronary_flow_capacity\patient_data\DUMMY_DCM'; % Steal path from unused acq
dcmpathOUT = '\\POLARIS\Data4_New\Seth\coronary_flow_capacity\patient_data\SUBJECT007\CFC_Vitrea\'; % Save in main folder for study --> use for all ACQs
dcmTitle = 'Subject007_CFC';
WriteDICOM(im_vol, dcmpathIN, dcmpathOUT, dcmTitle, im_vol);

%1

f_handle = figure();
X1 = stress_map1(stress_map1 ~= 0);
Y1 = cfr1(cfr1 ~= 0);
% Scatter plot 
figure();
title('CFR vs Stress Perfusion');
hold on;
scatter(X1, Y1);
xlabel('STRESS (mL/min/g)');
ylabel('CFR');

%2

X2 = stress_map2(stress_map2 ~= 0);
Y2 = cfr2(cfr2 ~= 0);
% Scatter plot 
figure();
title('CFR vs Stress Perfusion');
hold on;
scatter(X2, Y2);
xlabel('STRESS (mL/min/g)');
ylabel('CFR');

%3 

X3 = stress_map3(stress_map3 ~= 0);
Y3 = cfr3(cfr3 ~= 0);
% Scatter plot 
figure();
title('CFR vs Stress Perfusion');
hold on;
scatter(X3, Y3);
xlabel('STRESS (mL/min/g)');
ylabel('CFR');


