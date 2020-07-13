%Mateusz Polak

clc 
clear all
close all

%% IN DATA

% refractive index of a surrounding medium
n0 = 1;

% wavelength
lambda = 0.633;

% pixel size after scaling by transverse magnification of an imaging system
dx = 3.45/10;

% parameters of axial scanning
z_min = 20;%60
z_max = 160;
dz = 4;

% size of zeropadded image (P>=max(Ny,Nx))
% P = 2048*1.5;
P = 2500;

% size of filter kernel
filter_size = 9;

% type of object
objtype = 'phase';
% objtype ='amplitude';

%% main body

%%% loading of an object wave
folder_name = 'C:\Users\hologram\';
data1 = double(imread([folder_name,'holo1.bmp']));
data2 = double(imread([folder_name,'holo2.bmp']));
data3 = double(imread([folder_name,'holo3.bmp']));
data4 = double(imread([folder_name,'holo4.bmp']));
data5 = double(imread([folder_name,'holo5.bmp']));
ph = atan2(2*(data2-data4),(2*data3-data1-data5));
amp = (1/4)*sqrt(4*(data2-data4).^2 + (2.*data3-data1-data5).^2);
u=amp.*exp(1i*ph);


% ROI selection
ROI = ROISelection(abs(u));

dz_vec = z_min:dz:z_max;
steps_count = length(dz_vec);
focus_curve = zeros(steps_count,1);

%%% zeroppading and enveloping with Tukey window 
[Ny,Nx] = size(u);
wy = tukeywin(Ny,0.001);wx = tukeywin(Nx,0.001);Tukey = wy*wx';
uin = ones(P,P);
uin( P/2-Ny/2+1 : P/2+Ny/2, P/2-Nx/2+1 : P/2+Nx/2 ) = u.*Tukey;

x=(-Nx/2+1:Nx/2)*dx;
y=(-Ny/2+1:Ny/2)*dx;

%%% autofocusing loop
for ii = 1:steps_count

    % propagation
%     uout = PropagatePlaneWaveDecC2DNxNy(uin,dz_vec(ii),n0,lambda,dx);
     uout = AS_propagate(uin,dz_vec(ii),n0,lambda,dx);

    uout = uout( P/2-Ny/2+1:P/2+Ny/2 , P/2-Nx/2+1:P/2+Nx/2 );
    amp_out = abs(uout(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1));

    % filtering
%   amp_out = medfilt2(abs(uout),[filter_size,filter_size]);%median filter
    amp_out  =  ConvFiltr2D(amp_out,filter_size,'hann');% convolution filter
    
    %figure(1);imagesc(amp_out);colormap gray;title(num2str(dz_vec(ii)));pause(0.1);
    figure(1)
    %figure('Name','Focusing')
    imagesc(amp_out);colormap gray;title(num2str(dz_vec(ii)));pause(0.1);
    
    % evaluation of focusing condition by calculating of the image sharpness
    switch objtype
        case 'amplitude'
%             focus_curve(ii) = DERIVE(amp_out);
%             focus_curve(ii) = RANGE(amp_out);
            focus_curve(ii) = LogDerive(amp_out);
            % othe focus evaluation method: RANGE, LogDerive
        case 'phase'
%             focus_curve(ii) = VAR(amp_out);
            focus_curve(ii) = NORM_L1(amp_out);
            % othe focus evaluation method: NORM_L1
    end
    
end

% determination of minimum/maximum value 
[maxVal, Imax] = max(focus_curve);
[minVal, Imin] = min(focus_curve);
focus_curve = ( focus_curve - ones(steps_count,1)*minVal )./( ones(steps_count,1)*(maxVal - minVal) );  

switch objtype
    case 'amplitude'
        [junk peak_pos] = max(focus_curve);
    case 'phase'
        [junk peak_pos] = min(focus_curve);
end

z_rec = dz_vec(peak_pos);
figure;plot(dz_vec,focus_curve);axis tight;xlabel('dz \mum');ylabel('norm. sharpness');
saveas(gcf,'img_save\focus_curve.bmp');

%%propagation to the image plane
%uout = PropagatePlaneWaveDecC2DNxNy(uin,z_rec,n0,lambda,dx);
 uout = AS_propagate(uin,z_rec,n0,lambda,dx);
uout = uout( P/2-Ny/2+1:P/2+Ny/2 , P/2-Nx/2+1:P/2+Nx/2 );

% displaying
amp_in = abs(u(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1));
amp_out = abs(uout(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1));

figure('Name','Hologram and image plane')
subplot(1,2,1);imagesc(amp_in);colormap gray;title('hologram plane');axis off;axis image
subplot(1,2,2);imagesc(amp_out);colormap gray;title('image plane'); axis off;axis image
saveas(gcf,'img_save\hol_img_plane.bmp');

amp = abs(uout);
ph = angle(uout);

tic

%%%%%%%%%%%%%%Choose phase unwrapping method%%%%%%%%%%%%%%%%%%%%%%%%%%%%

unph = Miguel_unwrap(ph);
unwrap_name = "Miguel";

%unph = TIE_unwrap2DPdct(ph,dx,lambda);
%unwrap_name = "TIE";

time_of_exec = toc

figure('Name','Amplitude')
imagesc(x,y,amp);axis image;title(unwrap_name);
saveas(gcf,'img_save\amp_' + unwrap_name + '.bmp');

figure('Name','Phase')
imagesc(x,y,ph);axis image;title(unwrap_name);
saveas(gcf,'img_save\ph_' + unwrap_name + '.bmp');

figure('Name','Unwrapped phase')
imagesc(x,y,unph);axis image;title(unwrap_name);
saveas(gcf,'img_save\unph_' + unwrap_name + '.bmp');

figure('Name','Unwrapped phase 3D')
title(unwrap_name);
surfl(unph);colormap(pink);shading interp

%%%%%%%%%%%%%%%%%%%%%%%%% Circle detection %%%%%%%%%%%%%%%%%%%%%%%
image = mat2gray(unph);

figure('Name','Grayscale')
title(unwrap_name);
imshow(image)
saveas(gcf,'img_save\grayscale_' + unwrap_name + '.bmp');

img_test = image;
img_test = imbinarize(img_test);
figure('Name','Binary image')
title(unwrap_name);
imshow(img_test)
saveas(gcf,'img_save\binary_' + unwrap_name + '.bmp');
 
 
figure('Name','ROI circle detection')
img223 = imagesc(image);colormap gray;title('ROI,choose at least 2x2');
[Xroi,Yroi,Wroi,PosROI] = imcrop(img223);
PosROI = round(PosROI);
if ( mod(PosROI(3),2) == 1 ), PosROI(3) = PosROI(3)+1; end; 
if ( mod(PosROI(4),2) == 1 ), PosROI(4) = PosROI(4)+1; end;

image = imbinarize(image);
image(1:PosROI(2),:)=0;
image(PosROI(2)+PosROI(4):end,:)=0;
image(:,1:PosROI(1))=0;
image(:,PosROI(1)+PosROI(3):end)=0;

image = mat2gray(image);
sigma = 2; % standard deviation 
im_gauss = imgaussfilt(image, sigma); % gaussian filtering

figure('Name','Binary cleaned image')
title(unwrap_name);
imshow(image)
saveas(gcf,'img_save\binary_cleaned_' + unwrap_name + '.bmp');

figure('Name','Image after Gaussian filtering')
title(unwrap_name);
imshow(im_gauss)
saveas(gcf,'img_save\gauss_' + unwrap_name + '.bmp');

% Laplacian of Gaussian
ker = ones(3); % 
ker = ker *(-1);
ker(2,2) = 8;

lap = conv2(im_gauss, ker, 'same');
figure('Name','Edge detection')
title(unwrap_name);
imshow(lap)
saveas(gcf,'img_save\edge_' + unwrap_name + '.bmp');

% range of predicted radiuses values
r_min = 90;
r_max = 130;

[X, Y] = size(image); % accumulator array

step = 5; % step
accumulator = zeros(X, Y, (r_max-r_min)/step+1); %proper accumulator with intensity

% iteration over every non-black ( =edge) pixel of an image
% for every such pixel, a circle is drawn in accumulator space
% with chosen R

for rad = r_min:step:r_max
    for x = 1:X
        for y= 1:Y
            for theta = 0:2:360
                if  lap(x, y) ~= 0
                    n = floor(x - rad * cos(theta));
                    m = floor(y - rad * sin(theta));
                    if( n > 0 && n < X && m > 0 && m < Y) %antiexciding of boundaries
                        accumulator(n, m, (rad-r_min)/step+1) =  accumulator(n, m, (rad-r_min)/step+1) + 1;
                    end
                end
            end
        end
    end
    figure(2)
    imshow(mat2gray(accumulator(:, :,(rad-r_min)/step+1)));title(num2str(rad));pause(0.1);
    if rad == r_min || rad == r_max || rad == 100
        saveas(gcf, 'img_save\acc_' + unwrap_name + '_rad_' + rad + '.bmp');
    end
end

accumulator_max = accumulator > 179; % thresholding

CNTRS = []; % centers
r_val = []; % radius value 

i = 1;
% voting for centers
for rad = r_min:step:r_max
    reg = regionprops(accumulator_max(:,:,(rad-r_min)/step+1) , 'Centroid');
    cntr = vertcat(reg.Centroid);
    [g, h] = size(cntr);
    r_val(i:i+g-1,1) = rad;
    CNTRS = [CNTRS; cntr];
    i = i+g;
end

all_centers = floor([CNTRS, r_val]);
[g, h] = size(all_centers);

coord = []; % circles radiuses and coordinates container
for i = 1:g
    cn = all_centers(i, :); % considered center
    
    if cn == [0,0,0] % skips rest of  current loop goes to next loop iteration
        continue 
    end
    
    cent_counter = []; % array for repetitions
    for j = i+1:g
         comparator = all_centers(j, :);
         cn_x = abs(cn(1) - comparator(1));
         cn_y = abs(cn(2) - comparator(2));
         
         if (sqrt(cn_x*cn_x+cn_y*cn_y) < cn(3)) % is current centre inside circle
            cent_counter = [cent_counter; comparator];
            all_centers(j, :) = [0,0,0];
         end
         if j == g
             counter = size(cent_counter);
             if counter(1) == 0 % only one center inside circle = add coordinates to container
                 coord = [coord; cn];
             else % if there is more than one center in circle calculates average and add to container
                 coord = [coord; floor(mean([cn;cent_counter]))];
             end
         end
    end
end


[q, p] = size(coord); 
result = zeros(X, Y);
sorted = sort(coord(:, 1));
differences = diff(sorted);
av_dist = floor(mean(differences(differences > 50))/2);
for i = 1:q % obtained circles plot 
    for theta= 1:0.5:360
            g = floor(coord(i,1) - coord(i,3) * cos(theta));
            h = floor(coord(i,2) - coord(i,3) * sin(theta));
            result(h, g) = 1;
    end
    xx = coord(i,2); 
    yy = coord(i,1);
    result(xx,yy) = 1; % marks centers on circles plot

    final_r = coord(i,3); %radius
    
    unph_max = unph(xx, yy);
    coord(i, 4) = unph_max; % unph maximal values
    
    av_min = (unph(xx-av_dist, yy)+ unph(xx, yy-av_dist) + unph(xx+av_dist, yy) + unph(xx, yy+av_dist))/4;
    coord(i, 5) = av_min;% average of minimal values
    
    final_h = (unph_max-av_min)/(2*pi);
    coord(i, 6) = final_h; % heigth
    
    ROC = final_h/2 + (final_r*dx)*(final_r*dx)/(2*final_h);
    coord(i, 7) = ROC; % radius of curvature
    
    NA = final_h/(final_r*dx);
    coord(i, 8) = NA; % numerical aperture
    
    diam = final_r*dx*2;
    coord(i, 9) = diam; % diameter
    
    final_f = ROC/(1.45701-1);
    coord(i, 10) = final_f; % focal length

end

writematrix(coord, 'img_save\coords' + unwrap_name + '.xlsx');

figure('Name','Obtained circles and centers')
title(unwrap_name);
imshow(result)
imwrite(result, 'img_save\results' + unwrap_name + '.bmp');
