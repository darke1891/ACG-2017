obsFrequencies = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 58, 300, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 40, 10096, 34277, 2538, 3, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 7, 3630, 153049, 380404, 58091, 583, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 2, 165, 7789, 101376, 193866, 51394, 2281, 39, 1, 0, 0, 0, 0, 0, 0 ];
expFrequencies = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 2.59619e-30, 1.01572e-16, 3.80678e-07, 0.0312799, 0.23803, 0.00211255, 6.16923e-10, 5.25889e-21, 1.02887e-35, 0, 0, 0, 0, 0, 0; 0, 0, 0, 2.71083e-39, 5.0934e-28, 1.72872e-17, 1.07765e-08, 0.0227182, 59.1441, 288.031, 8.56479, 0.000308585, 1.77284e-11, 5.71627e-21, 7.06583e-32, 5.32625e-43, 0, 0, 0, 0; 2.39833e-40, 5.94872e-36, 5.90775e-30, 5.13781e-23, 7.85394e-16, 4.36517e-09, 0.00231919, 40.4231, 10057.6, 34168.5, 2518.64, 2.18243, 3.44326e-05, 2.53233e-11, 2.71935e-18, 1.79855e-25, 3.65719e-32, 1.11165e-37, 2.32072e-41, 5.32625e-43; 3.76851e-21, 9.45082e-19, 2.17612e-15, 1.89195e-11, 2.50969e-07, 0.00225195, 6.9536, 3589.95, 153126, 380198, 57783.4, 532.198, 0.512491, 0.000106599, 9.57423e-09, 7.55825e-13, 1.23997e-16, 1.04839e-19, 1.0969e-21, 3.25182e-22; 0.166096, 0.166096, 0.166115, 0.167023, 0.222157, 3.24992, 161.476, 7834.06, 101492, 194393, 51444.6, 2250.65, 39.7376, 0.959095, 0.179682, 0.166329, 0.166102, 0.166096, 0.166096, 0.166096 ];
colormap(jet);
clf; subplot(2,1,1);
imagesc(obsFrequencies);
title('Observed frequencies');
axis equal;
subplot(2,1,2);
imagesc(expFrequencies);
axis equal;
title('Expected frequencies');