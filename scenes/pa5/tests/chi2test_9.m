obsFrequencies = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 1163, 1147, 1156, 1168, 1229, 1139, 1107, 1228, 1220, 1186, 1221, 1239, 1170, 1283, 71812, 8424, 1215, 1198, 1230, 1236; 3685, 3627, 3558, 3647, 3480, 3561, 3590, 3496, 3603, 3514, 3661, 3596, 3544, 5372, 62453, 20656, 3606, 3620, 3649, 3642; 6014, 5956, 5999, 6004, 6106, 6033, 6162, 5988, 5877, 6070, 6011, 6037, 6008, 11096, 45358, 24500, 6624, 6030, 5988, 5974; 8419, 8401, 8381, 8584, 8406, 8240, 8430, 8411, 8334, 8418, 8223, 8342, 8882, 13785, 26254, 19906, 9921, 8571, 8408, 8388; 10754, 10760, 10715, 10670, 10805, 10951, 10902, 10689, 10859, 10784, 10837, 10659, 11278, 12634, 14377, 13582, 11921, 11020, 10836, 10763 ];
expFrequencies = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 1200, 1200, 1200, 1200, 1200, 1200, 1200, 1200, 1200, 1200, 1200, 1200, 1200, 1261.78, 71375, 8346.58, 1200.02, 1200, 1200, 1200; 3600, 3600, 3600, 3600, 3600, 3600, 3600, 3600, 3600, 3600, 3600, 3600, 3600.42, 5276.57, 62747, 20602.9, 3630.53, 3600, 3600, 3600; 6000, 6000, 6000, 6000, 6000, 6000, 6000, 6000, 6000, 6000, 6000, 6000.07, 6049.99, 11110.6, 45118.2, 24413.6, 6542.21, 6001.78, 6000, 6000; 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400.01, 8400.18, 8410.94, 8800.94, 13921.5, 26330.6, 19908.6, 9948.84, 8461.85, 8401.21, 8400.02; 10803.3, 10802.6, 10802.5, 10802.4, 10802.4, 10802.4, 10802.4, 10802.5, 10802.9, 10804.7, 10815.3, 10883.1, 11282.4, 12663.6, 14260.5, 13538.5, 11764, 10989.5, 10831.9, 10807.4 ];
colormap(jet);
clf; subplot(2,1,1);
imagesc(obsFrequencies);
title('Observed frequencies');
axis equal;
subplot(2,1,2);
imagesc(expFrequencies);
axis equal;
title('Expected frequencies');