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
last_projection_angle=180;
delta_theta=0.5;  

% set the number of theta views
theta=first_projection_angle:delta_theta:last_projection_angle; 

%
% step 1 create projection of image
rad_I=radon(I,theta); %rows of rad_I correspond to L colms correspond to theta
figure('Name','g(l,theta)'); imagesc(rad_I); title('projections g(l,theta) - Radon Transform'); axis('square'); xlabel('projection angle theta'); ylabel('linear displacement - l');
[N_l,N_theta]=size(rad_I);

%% inverse FT sim 

% step 1: fourier transform of image projections
 rad_I_FFT = fftshift(fft(rad_I)); 
 %rad_I_FFT_shifted = fftshift(rad_I);
 
 [Num_L, Num_theta] = size(rad_I_FFT);
 
 %L_array = (1:Num_L - (Num_L+1)/2)'; %check this
 L_array = (1:Num_L)';
 
 figure('Name','FFT Original')
 plot(1:numel(rad_I_FFT(:,round(90/delta_theta))),abs(rad_I_FFT(:,round(90/delta_theta))))
 title('FFT Original at g(l,pi/2)')
 
 %step 2: insert this value at the correct line in the 2D FT Plan
    % create a zero padded 2D array 
    % create an x,y value for each L,theta value 
    %place these values into the 2D array 
    theta_reshape = repmat(theta,Num_L,1);
    theta_reshape = theta_reshape .* (pi/180);
    
    L_array_reshape = repmat(L_array,1, Num_theta);
    
    %calculate x and y values 
    x = L_array_reshape .* cos(theta_reshape);
    x = reshape(x,[],1);
    
    y = L_array_reshape .* sin(theta_reshape);
    y = reshape(y,[],1); 
    
    %reshape rad_I_FFT to match x and y
    v = reshape(rad_I_FFT,[],1); 
    
    [xq, yq] = meshgrid(L_array,L_array);
    
    ft_plane = griddata(x,y,v,xq,yq);
    
    %plot 2D_ft with interpolation
    %mesh(xq,yq,ft_plane)
    %hold on
    %plot3(x,y,v,'o')
    
    reconstructed_image = ifft2(ft_plane);
    
    figure('Name','Reconstructed Image','Position',[1 420 400 400]);
    imagesc(real(reconstructed_image)); title('Reconstructed image'); axis('square'); colormap('gray');
 
 
 
 