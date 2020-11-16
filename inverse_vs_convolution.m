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
last_projection_angle=180;
delta_theta=0.1;


%% from elliott's code
% %change the level of the noise on raw data...
% %
% data_peak_to_noise_ratio=2.0;

% change text to reflect parameters of the particular run
% figure('Name','parameters of run','Position',[1201 604 240 201]);
% text(0.1,0.9,['delta theta = ',num2str(delta_theta)],'Fontsize', 14)
% text(0.1,0.7,'angles 0 through 180','Fontsize', 14)
% text(0.1,0.5,['data SNR = ',num2str(data_peak_to_noise_ratio)],'Fontsize', 14)
% text(0.1,0.3,'recon kernel = hamming','Fontsize', 14)
%%


filtinv = 0; 
convolution = 1;

% set the number of theta views
theta=first_projection_angle:delta_theta:last_projection_angle;
%
% step 1 create projection of image
rad_I=radon(I,theta);
figure('Name','g(l,theta)'); imagesc(rad_I); title('projections g(l,theta) : Radon Transform'); axis('square'); xlabel('projection angle theta'); ylabel('linear displacement - l'); colormap('gray');
[N_l,N_theta]=size(rad_I);
%%
if filtinv == 1 
    % step 2: fourier transform of projected image
    rad_I_FFT = fft(rad_I);   
    Pr = real(rad_I_FFT); 
    figure, plot(Pr);
    
    % step 3: TODO- multiply by filter 
    
    
    %step 4: inverse fourier transform 
    rad_I_transformed = ifft(rad_I_FFT);
    %plot difference between original and inverse
    figure, imagesc(rad_I_transformed-rad_I);
end 

%%
if convolution == 1 
     % step 2: create filter
     %step 3: fourier transform filter 
     
      %y = sinc(-4*pi:.1:4*pi);
      
      
      hanning_filter = hann(64);
      hanning_filter = [hanning_filter(1:end-1); flip(hanning_filter)];
      
      %adopted from fltered_back_test
      freqs = linspace(-1,1,N_theta);
      ram_lak = abs( freqs ); 
      figure, plot(ram_lak)
      title('Filter in Frequency Domain') 
      xlabel('Frequency') 
      ylabel('Magnitude') 
      
      
      filter = ifft(fftshift(ram_lak));
      filter = real(fftshift(filter)); 
      figure, plot(filter);
      title('Filter in Spatial Domain') 
      xlabel('n') 
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
   
     %figure('Name','g(l,theta)'); imagesc(conv_2d); title('projections g(l,theta) : Radon Transform + Convolution'); axis('square'); xlabel('projection angle theta'); ylabel('linear displacement - l');
     
     %figure('Name','g(l,theta)'); imagesc(iradon(rad_I, theta, 'none')); title('projections g(l,theta) : Backprojection'); axis('square'); xlabel('projection angle theta'); ylabel('linear displacement - l');
   % figure('Name','g(l,theta)'); imagesc(iradon(conv_2d, theta, 'none')); title('projections g(l,theta) : Convolution Backprojection'); axis('square'); xlabel('projection angle theta'); ylabel('linear displacement - l');
     %   figure, imagesc(iradon(conv_back_2d,theta, 'none')); title('convolution reconstruction');
    
end 
%%
% step 5. create back projection 


inv_rad_I=iradon(conv_back,theta,'None');

crop_iradI = inv_rad_I(round(end/2)-125:round(end/2)+125,round(end/2)-125:round(end/2)+125);

figure('Name','Reconstructed Image','Position',[1 420 400 400]);
imagesc(crop_iradI); title('Reconstructed image'); axis('square'); colormap('gray');

%figure('Name','Original - Reconstructed'); imagesc(I - inv_rad_I(2:257, 2:257)); title('Original - Reconstructed'); axis('square'); colormap('gray');




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