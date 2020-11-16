%TO DO: add an amazing simulation

% ft_sim.m
% Jacob Prince, BENG280A, 2020
%

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

%% inverse FT sim 

% step 1: fourier transform of image projections
 rad_I_FFT = ((fft((fftshift(rad_I))))); 
 
 %rad_I_FFT_shifted = fftshift(rad_I);
 
 [Num_L, Num_theta] = size(rad_I_FFT);
 
  for i = 1:N_theta
     rad_I_FFT_shift(:,i) = fftshift(rad_I_FFT(:,i));
 end
 
 %L_array = (1:Num_L - (Num_L+1)/2)'; %check this
 L_array = (1:Num_L)';
 L_array = L_array - (Num_L+1)/2;
 
 figure('Name','FFT Original')
 %plot(1:numel(rad_I_FFT_shift(:,round(90/delta_theta))),abs(rad_I_FFT_shift(:,round(90/delta_theta))))
 title('FFT Original at g(l,pi/2)')
 
 %step 2: insert this value at the correct line in the 2D FT Plan
    % create a zero padded 2D array 
    % create an x,y value for each L,theta value 
    %place these values into the 2D array 
    theta_reshape = repmat(theta,Num_L,1);
    theta_reshape = theta_reshape;
    
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
    
    [xq, yq] = meshgrid(L_array(1):0.1:L_array(end),L_array(1):0.1:L_array(end));
    
    ft_plane = griddata(x,y,v,xq,yq);
    
    ft_plane(isnan(ft_plane)) = 0;
    
    %plot 2D_ft with interpolation
    %mesh(xq,yq,abs(ft_plane))
    %hold on
   % plot3(x,y,abs(v),'o')
    
    figure
    imagesc(abs(ft_plane))
    
    reconstructed_image = fftshift(ifft2((fftshift(ft_plane))));
    
    figure('Name','Reconstructed Image','Position',[1 420 400 400]);
    imagesc(abs(reconstructed_image)); title('Reconstructed image'); axis('square'); colormap('gray');
 
 
 
 