obsFrequencies = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 817, 789, 834, 863, 951, 1194, 2146, 4630, 8630, 12993, 14562, 12079, 7510, 3688, 1839, 1134, 943, 829, 808, 791; 2447, 2520, 2676, 2787, 3165, 3921, 5726, 8489, 11651, 14796, 15597, 14009, 11053, 7455, 5058, 3741, 3089, 2667, 2570, 2569; 4689, 4814, 4885, 5417, 6087, 7256, 9091, 11402, 14160, 15547, 16250, 15470, 13446, 10855, 8617, 6821, 5852, 5123, 4849, 4783; 7735, 7701, 7980, 8545, 9637, 10652, 12217, 14029, 15420, 16500, 17001, 16414, 15152, 13658, 11887, 10386, 9373, 8590, 8002, 7607; 11783, 12040, 12079, 12584, 13319, 14013, 14780, 15625, 16223, 16510, 16685, 16774, 16012, 15322, 14551, 13654, 13089, 12615, 12188, 11704 ];
expFrequencies = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 807.995, 810.145, 818.359, 844.539, 932.119, 1237.97, 2186.41, 4559.19, 8716.88, 13055.7, 14587.3, 12094.7, 7524.92, 3770.97, 1841.96, 1121.12, 898.628, 834.656, 815.284, 809.211; 2525.42, 2546, 2613.08, 2774.65, 3147.36, 3970.89, 5612.97, 8337.14, 11767.5, 14645.1, 15565.6, 14046, 10877.3, 7538.76, 5092.32, 3697.37, 3020.87, 2719.37, 2589.59, 2537.27; 4678.31, 4746.24, 4946.11, 5347.23, 6073.46, 7291.11, 9127.68, 11497.6, 13949.2, 15761.1, 16309.2, 15397.2, 13350.2, 10856.6, 8596.11, 6922.45, 5847.07, 5219.15, 4879.25, 4717.89; 7644.87, 7762.1, 8083.5, 8652.89, 9527.27, 10745.3, 12272.8, 13948, 15478.2, 16524.1, 16829.7, 16318.6, 15118.4, 13518.2, 11857.4, 10399.8, 9270.76, 8479.88, 7979.45, 7713.74; 11829.3, 11933.5, 12204.9, 12643.5, 13241.3, 13970.3, 14773.2, 15560.3, 16221.2, 16649.6, 16771.9, 16566.7, 16069.8, 15365.7, 14564.3, 13773, 13073.7, 12515.1, 12119.2, 11890.9 ];
colormap(jet);
clf; subplot(2,1,1);
imagesc(obsFrequencies);
title('Observed frequencies');
axis equal;
subplot(2,1,2);
imagesc(expFrequencies);
axis equal;
title('Expected frequencies');