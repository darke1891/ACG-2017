obsFrequencies = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 1234, 1216, 1171, 1186, 1200, 1230, 1175, 1181, 1201, 15572, 46539, 3264, 1152, 1180, 1151, 1219, 1156, 1208, 1260, 1180; 3671, 3623, 3733, 3462, 3536, 3656, 3544, 3600, 4423, 27270, 52391, 11284, 3627, 3655, 3672, 3598, 3552, 3664, 3533, 3621; 6102, 6119, 5985, 5986, 6080, 5992, 5957, 6121, 9619, 32744, 47782, 19382, 7029, 6085, 5994, 5953, 5889, 6051, 5963, 6107; 8377, 8093, 8317, 8544, 8683, 8393, 8405, 9197, 14459, 29123, 35649, 22066, 11525, 8560, 8360, 8394, 8242, 8281, 8424, 8386; 10671, 10925, 10783, 10795, 10977, 10778, 11279, 12286, 15134, 18898, 20004, 17227, 13497, 11814, 11074, 10945, 10899, 10801, 10752, 10818 ];
expFrequencies = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 1200, 1200, 1200, 1200, 1200, 1200, 1200, 1200, 1252.6, 15630.4, 46622.6, 3312.74, 1201.36, 1200, 1200, 1200, 1200, 1200, 1200, 1200; 3600, 3600, 3600, 3600, 3600, 3600, 3600, 3603.7, 4481.79, 27220.8, 52140, 11238.9, 3701.74, 3600.2, 3600, 3600, 3600, 3600, 3600, 3600; 6000, 6000, 6000, 6000, 6000, 6000.01, 6001.72, 6124.38, 9541.24, 32661.1, 47948.8, 19270, 6943.74, 6020.99, 6000.21, 6000, 6000, 5999.99, 6000, 6000; 8400, 8400, 8400, 8400.02, 8400.28, 8404.38, 8468.47, 9241.18, 14566.8, 29134.4, 35696.4, 22010.9, 11202.5, 8694.4, 8420.89, 8401.3, 8400.09, 8400, 8400, 8400; 10817.5, 10819.4, 10825, 10838.2, 10871.8, 10968.8, 11294.9, 12360.4, 15005.8, 18675.1, 19898, 17123.1, 13622.1, 11749.3, 11102.8, 10912.7, 10853, 10831, 10821.9, 10818.2 ];
colormap(jet);
clf; subplot(2,1,1);
imagesc(obsFrequencies);
title('Observed frequencies');
axis equal;
subplot(2,1,2);
imagesc(expFrequencies);
axis equal;
title('Expected frequencies');