%TO DO: add an amazing simulation

% ft_sim.m
% Jacob Prince, BENG280A, 2020
% Alex Postlmayr, BENG 280A, 2020

%adopted from Elliot's Code
clear all;
close all;
% Use coronary image

I = imread('example_coronary_patient_000.tif');
I = imresize(I, [256,256]); 
I = double(I);

figure('Name','original image','Position',[1 420 400 400]);
imagesc(I); title('original image'); axis('square'); colormap('gray');

% ADJUSTABLE PARAMETERS FOR IRADON - the set of projections used to produce the images

% adopted from Elliot's Code
first_projection_angle=0;
last_projection_angle=179;
delta_theta=1;  

% set the number of theta views
theta=first_projection_angle:delta_theta:last_projection_angle; 

%
% step 1 create projection of image
rad_I=(radon(I,theta)); %rows of rad_I correspond to L colms correspond to theta
figure('Name','g(l,theta)'); imagesc(rad_I); title('projections g(l,theta) - Radon Transform'); axis('square'); xlabel('projection angle theta'); ylabel('linear displacement - l');
[N_l,N_theta]=size(rad_I);

%% Add noise (adopted from Elliot's Code)
%change the level of the noise on raw data...
data_peak_to_noise_ratio=20;

data_peak=max(max(rad_I));
raw_noise=rand(N_l,N_theta);
figure('Name','raw noise','Position',[20 20 400 400]); imagesc(raw_noise); title('Raw Noise to add to projections'); axis('square');

scaled_noise=raw_noise*(data_peak/data_peak_to_noise_ratio); 
rad_I=rad_I+scaled_noise;

%% inverse FT sim 

% step 1: fourier transform of image projections
 
 [Num_L, Num_theta] = size(rad_I);
 
 tic
 
  for i = 1:N_theta
     rad_I_FFT_shift(:,i) = fftshift(fft(fftshift(rad_I(:,i))));
  end
 
 % Neded L_array from -183 to 183 to index all 367 "detectors"
 L_array = (-(Num_L-1)/2:(Num_L-1)/2)';
 
 figure('Name','FFT Original')
 plot(1:numel(rad_I_FFT_shift(:,round(30/delta_theta))),1+log(abs(rad_I_FFT_shift(:,round(30/delta_theta)))))
 title('log(FFT Original) at g(l,pi/6)')
 
 %step 2: insert this value at the correct line in the 2D FT Plan
    % create a zero padded 2D array 
    % create an x,y value for each L,theta value 
    %place these values into the 2D array 
    theta_reshape = repmat(theta,Num_L,1);
    
    L_array_reshape = repmat(L_array,1, Num_theta);
    
    %calculate x and y values 
    x = L_array_reshape .* cosd(theta_reshape);
    %x = reshape(x,[],1);
    
    %x_nan = sum(isnan(x));
    
    y = L_array_reshape .* sind(theta_reshape);
    
    
    %y = reshape(y,[],1); 
    
    %y_nan = sum(isnan(y));
    
    %reshape rad_I_FFT to match x and y
    v = rad_I_FFT_shift;
    %v = reshape(rad_I_FFT_shift,[],1); 
    %v_nan = sum(isnan(v));
    
  [xq, yq] = meshgrid(L_array(1):delta_theta:L_array(end),L_array(1):delta_theta:L_array(end));
  
%     ft_plane = griddata(x,y,v,theta_reshape,L_array_reshape);
    ft_plane = griddata(x,y,v,xq,yq);    

    ft_plane(isnan(ft_plane)) = 0;
    
    %plot 2D_ft with interpolation
    %mesh(xq,yq,abs(ft_plane))
    %hold on
   % plot3(x,y,abs(v),'o')
    
    figure('Name','ft_plane','Position',[1 420 400 400]);
    imagesc(1+log(abs(ft_plane))); axis('square'); colormap('gray');
    title('Log of Interpolated 2D Fourier Plane')
    xlabel('U')
    ylabel('V')
    
    reconstructed_image = fftshift(ifft2(fftshift(ft_plane)));
    
    toc
    
    crop_image = reconstructed_image(round(end/2)-125:round(end/2)+125,round(end/2)-125:round(end/2)+125);
    
    pretty_reconstruct = flip(imrotate(abs(crop_image),180),2); 
    
    figure('Name','Reconstructed Image','Position',[1 420 400 400]);
    imagesc(pretty_reconstruct); title('Reconstructed image'); axis('square'); colormap('gray');
    
%% calculate SNR and other program parameters

background = pretty_reconstruct(225:250,1:25);
figure('Name','Background','Position',[1 420 400 400]);
imagesc(abs(background)); title('Background'); axis('square'); colormap('gray');

max_signal_average = pretty_reconstruct(112:137,75:100);

SNRdiff = (mean(mean(max_signal_average))-mean(mean(background)))/std(reshape(background,[],1))
%SNRdiff = (mean(mean(max_signal_average)))/std(reshape(background,[],1))
 