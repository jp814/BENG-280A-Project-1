% filtered_vs_convolution_backprojection.m
% Use matlab radon/iradon with filter methods to demonstrate
% difference between filtered inverse backprojection and convoution
% reconstruction
%
% Alex Postlmayr, BENG280A, 2020



% ------------------ README ------------------------
% filtered inverse
% 1. project image 
% 2. take fourier
% 3. multiply by image
% 4. take inverse fourier
% 5. backproject (iradon w/ 'none' as filter)

% convolution 
% 1. project image 
% 2. create filter
% 3. fourier filter
% 4. convolute image projection w/ fourier filter
% 5. backproject (iradon w/ 'none' as filter)




%
clear all;
close all;
% Use coronary image

I = imread('example_coronary_patient_000.tif');
I = imresize(I, [256,256]); 
I = double(I);

%I = phantom('Modified Shepp-Logan',200);

figure('Name','original image','Position',[1 420 400 400]);
imagesc(I); title('original image'); axis('square'); colormap('gray');

% ADJUSTABLE PARAMETERS FOR IRADON - the set of projections used to produce the images

first_projection_angle=0;
last_projection_angle=179;
delta_theta=1;

% set the number of theta views
theta=first_projection_angle:delta_theta:last_projection_angle;
%
% step 1 create projection of image
rad_I=radon(I,theta);
figure('Name','g(l,theta)'); imagesc(rad_I); title('projections g(l,theta) : Radon Transform'); axis('square'); xlabel('projection angle theta'); ylabel('linear displacement - l'); colormap('gray');

figure('Name','g(l,pi/2)');
plot(1:numel(squeeze(rad_I(:,round(90/delta_theta)))),squeeze(rad_I(:,round(90/delta_theta))))
title('g(l,pi/2)');
xlabel('L')
ylabel('Signal')

[N_l,N_theta]=size(rad_I);

tic

%% Add noise (adopted from Elliot's Code)
%change the level of the noise on raw data...
data_peak_to_noise_ratio=20;

data_peak=max(max(rad_I));
raw_noise=rand(N_l,N_theta);
figure('Name','raw noise','Position',[20 20 400 400]); imagesc(raw_noise); title('Raw Noise to add to projections'); axis('square');

scaled_noise=raw_noise*(data_peak/data_peak_to_noise_ratio); 
rad_I=rad_I+scaled_noise;

     % step 2: create filter
     %step 3: fourier transform filter 
     
      %y = sinc(-4*pi:.1:4*pi);
      
      
      hanning_filter = hann(N_theta/2);
      hanning_filter = [hanning_filter(1:end-1); flip(hanning_filter)];
      
      %adopted from fltered_back_test (need to cite)
      freqs = linspace(-1,1,N_theta);
      ram_lak = abs( freqs ); 
      figure, plot(hanning_filter)
      title('Filter in Frequency Domain') 
      xlabel('Frequency') 
      ylabel('Magnitude') 
      
      
      filter = ifft(fftshift(ram_lak));
      filter = real(fftshift(filter)); 
      figure, plot(filter);
      title('Filter in Spatial Domain') 
%       xlabel('n') 
      ylabel('g(n)') 
     
     % step 4: convolute image projection with fourier transform filter
     
     rad_I = (rad_I); % I added this
     
     conv_2d = zeros(367,361);
     for theta_cnt = 1:N_theta
         conv_back(:,theta_cnt) = conv(squeeze(rad_I(:,theta_cnt)), filter);
         %conv_back_cropped = conv_back(length(filter)/2: length(conv_back)-length(filter)/2);
        %conv_2d(:, theta_cnt) = conv_back_cropped;
     end
     
     figure('Name','Convolution Backprojection','Position',[1 420 400 400]);
    imagesc(conv_back); title('Convolution Backprojection'); axis('square'); colormap('gray');
    
    figure('Name','g(l,pi/2)*c(l)');
    plot(1:numel(squeeze(conv_back(:,round(90/delta_theta)))),squeeze(conv_back(:,round(90/delta_theta))))
    title('g(l,pi/2)*c(l)');
    xlabel('L')
    ylabel('Signal')
    
     %figure('Name','g(l,theta)'); imagesc(conv_2d); title('projections g(l,theta) : Radon Transform + Convolution'); axis('square'); xlabel('projection angle theta'); ylabel('linear displacement - l');
     
     %figure('Name','g(l,theta)'); imagesc(iradon(rad_I, theta, 'none')); title('projections g(l,theta) : Backprojection'); axis('square'); xlabel('projection angle theta'); ylabel('linear displacement - l');
   % figure('Name','g(l,theta)'); imagesc(iradon(conv_2d, theta, 'none')); title('projections g(l,theta) : Convolution Backprojection'); axis('square'); xlabel('projection angle theta'); ylabel('linear displacement - l');
     %   figure, imagesc(iradon(conv_back_2d,theta, 'none')); title('convolution reconstruction');
    
%%
% step 5. create back projection 


inv_rad_I=iradon(conv_back,theta,'linear','none');
%inv_rad_I = iradon(rad_I,theta,'hamming');
toc
crop_iradI = inv_rad_I(round(end/2)-125:round(end/2)+125,round(end/2)-125:round(end/2)+125);

figure('Name','Reconstructed Image','Position',[1 420 400 400]);
imagesc(abs(crop_iradI)); title('Reconstructed image'); axis('square'); colormap('gray');

%figure('Name','Original - Reconstructed'); imagesc(I - inv_rad_I(2:257, 2:257)); title('Original - Reconstructed'); axis('square'); colormap('gray');

%% calculate SNR and other program parameters

background = crop_iradI(225:250,1:25);
figure('Name','Background','Position',[1 420 400 400]);
imagesc(abs(background)); title('Background'); axis('square'); colormap('gray');

max_signal_average = crop_iradI(112:137,75:100);

SNRdiff = (mean(mean(max_signal_average))-mean(mean(background)))/std(reshape(background,[],1))
%SNRdiff = (mean(mean(max_signal_average)))/std(reshape(crop_iradI,[],1))
% 
% 
% %% elliott's code again for adding noise & filters 
% % try different filter kernels for the inverse radon transform, for example
% % inv_rad_I=iradon(rad_I,theta,'hamming');
% 
% 
% figure('Name','noise free iradon{g(l,theta)}','Position',[1 420 400 400]); imagesc(inv_rad_I); title('Inverse Radon of noise free Projections'); axis('square'); colormap('gray')
% %
% % add noise to the "raw data"
% %
% data_peak=max(max(rad_I));
% raw_noise=rand(N_l,N_theta);
% figure('Name','raw noise','Position',[20 20 400 400]); imagesc(raw_noise); title('Raw Noise to add to projections'); axis('square');
% %
% scaled_noise=raw_noise*(data_peak/data_peak_to_noise_ratio); 
% nse_rad_I=rad_I+scaled_noise;
% %
% % try smoothing the raw data to increase the SNR
% % 
% % invent your own filters to experiment with smoothing the data in
% % different directions...
% %
% smoothing=1
% if smoothing
% h=ones(4,4)/16;
% filt_nse_rad_I=conv2(nse_rad_I,h,'same');
% % note: the 'same' parameter cuts off extra data post convolution
% nse_rad_I=filt_nse_rad_I;
% figure('Name','g(l,theta) * filter','Position',[401 420 400 400]); 
% imagesc(filt_nse_rad_I); title('3x3 Filtered Projections'); axis('square');
% end
%     
% figure('Name','g(l,theta) + scaled_noise','Position',[801 420 400 400]); 
% imagesc(nse_rad_I); title('Noisy Projections'); axis('square');
% %
% inv_nse_rad_I=iradon(nse_rad_I,theta,'ram-lak');
% %
% %
% figure('Name','iradon{ g(l,theta) + scaled_noise}','Position',[401 420 400 400]); 
% imagesc(inv_nse_rad_I); title('Inverse Radon of Noisy Projections'); axis('square'); colormap(gray);
% 
% 
% %%