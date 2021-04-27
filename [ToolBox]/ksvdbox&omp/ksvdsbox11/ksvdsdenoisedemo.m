function ksvdsdenoisedemo
%KSVDSDENOISEDEMO Sparse K-SVD denoising demonstration.
%  KSVDSDENISEDEMO reads an image, adds random white noise and denoises it
%  using Sparse K-SVD denoising. The input and output PSNR are compared,
%  and the trained dictionary is displayed.
%
%  To run the demo, type KSVDSDENISEDEMO from the Matlab prompt.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  October 2009


disp(' ');
disp('  **********  Sparse K-SVD Denoising Demo  **********');
disp(' ');
disp('  This demo reads a 3-D volume, adds random Gaussian noise, and denoises it');
disp('  using 3-D Sparse K-SVD denoising. The function displays a slice from');
disp('  the original, noisy, and denoised volumes, and computes the output PSNR.');
disp(' ');


%% prompt user for volume %%

pathstr = fileparts(which('ksvdsdenoisedemo'));
dirname = fullfile(pathstr, 'volumes', '*.mat');
imglist = dir(dirname);

disp('  Available volumes:');
disp(' ');
for k = 1:length(imglist)
  load(fullfile(pathstr, 'volumes', imglist(k).name), 'vname');
  printf('  %d. %s (%s)', k, vname, imglist(k).name);
end
disp(' ');

imnum = 0;
while (~isnumeric(imnum) || ~iswhole(imnum) || imnum<1 || imnum>length(imglist))
  imnum = input(sprintf('  Volume to denoise (%d-%d): ', 1, length(imglist)), 's');
  imnum = sscanf(imnum, '%d');
end

imgname = fullfile(pathstr, 'volumes', imglist(imnum).name);



%% generate noisy volume %%

sigma = 50;

disp(' ');
disp('Generating noisy volume...');

load(imgname);
s = double(s);
s = s/max(s(:))*255;   % normalize to [0,255]

n = randn(size(s)) * sigma;
sn = s + n;
clear n;


%%% 2-D denoising parameters %%%

params2d.x = sn(:,:,slicenum);              % slice to denoise
params2d.blocksize = 8;                     % size of block to operate on
params2d.dictsize = 100;                    % size of trained dictionary

params2d.basedict{1} = odctdict(8,10);      % the separable base dictionary -
params2d.basedict{2} = odctdict(8,10);      % basedict{i} is the base dictionary
                                            % for the i-th dimension

params2d.initA = speye(params2d.dictsize);  % initial A matrix (identity)

params2d.Tdict = 6;                         % sparsity of each trained atom
params2d.sigma = sigma;                     % noise power
params2d.maxval = 255;                      % maximum sample intensity value
params2d.trainnum = 30000;                  % number of training blocks
params2d.iternum = 15;                      % number of training iterations
params2d.stepsize = 1;                      % distance between adjecent blocks to denoise
params2d.lambda = 0;                        % lagrange multiplier
params2d.gain = 1.15;                       % noise gain



%%% 3-D denoising parameters %%%

params3d.x = sn;                            % volume to denoise
params3d.blocksize = 8;                     % size of block to operate on
params3d.dictsize = 1000;                   % size of trained dictionary

params3d.basedict{1} = odctdict(8,10);      % the separable base dictionary -
params3d.basedict{2} = odctdict(8,10);      % basedict{i} is the base dictionary
params3d.basedict{3} = odctdict(8,10);      % for the i-th dimension

params3d.initA = speye(params3d.dictsize);  % initial A matrix (identity)

params3d.Tdict = 16;                        % sparsity of each trained atom
params3d.sigma = sigma;                     % noise power
params3d.maxval = 255;                      % maximum sample intensity value
params3d.trainnum = 80000;                  % number of training blocks
params3d.iternum = 15;                      % number of training iterations
params3d.stepsize = 2;                      % distance between adjecent blocks to denoise
params3d.lambda = 0;                        % lagrange multiplier
params3d.gain = 1.04;                       % noise gain


% denoise!

printf('\n\nPerforming 3-D Sparse K-SVD denoising...\n');
sout3 = ksvdsdenoise(params3d);

printf('\nPerforming 2-D Sparse K-SVD denoising of slice %d...', slicenum);
sout2 = ksvdsdenoise(params2d,0);


%% show results %%

figure; imshow(s(:,:,slicenum)/params2d.maxval);
title(sprintf('Original volume (showing slice %d)', slicenum));

figure; imshow(sn(:,:,slicenum)/params2d.maxval); 
title(sprintf('Noisy volume, PSNR = %.2fdB (showing slice %d)', 20*log10(params2d.maxval * sqrt(numel(s)) / norm(s(:)-sn(:))), slicenum ));

figure; imshow(sout2/params2d.maxval);
title(sprintf('2-D Sparse K-SVD, PSNR: %.2fdB (showing slice %d)', 20*log10(params2d.maxval * sqrt(numel(sout2)) / norm(s(:,:,slicenum)-sout2,'fro')), slicenum ));

figure; imshow(sout3(:,:,slicenum)/params2d.maxval);
title(sprintf('3-D Sparse K-SVD, PSNR: %.2fdB (showing slice %d)', 20*log10(params2d.maxval * sqrt(numel(s)) / norm(s(:)-sout3(:))), slicenum ));
