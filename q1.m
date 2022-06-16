clear all;
currentPath = pwd;
data_folder = [currentPath,'\data\'];
localization_folder = [currentPath,'\localization\'];
IQfiles = dir([data_folder '*.mat']);

iq = load([IQfiles(1).folder filesep IQfiles(1).name],'x').x;

load('template.mat','PSF');

SizeOfBloc = size(iq);

Nbuffers = numel(IQfiles);
res = 10;
clear iq

coimage_tot = {};
parfor hhh = 1:min(999,Nbuffers)
    fprintf('Processing bloc %d/%d\n',hhh,Nbuffers);
    % Load IQ data (or other kind of images without compression)
    iq = load([IQfiles(hhh).folder filesep IQfiles(hhh).name],'x').x;
    coimage = my_local(abs(iq),PSF,res);
    %coimage_tot{hhh} = coimage;
    [r,c,f] = ind2sub(size(coimage),find(coimage == 1));
    sr = size(r);
    mattracking = nan(sr(1),3);
    mattracking(:,1) = r;
    mattracking(:,2) = c;
    mattracking(:,3) = f;
    parsave(localization_folder+"mattracking"+string(hhh)+".mat",mattracking);
end

%imshow(coimage_tot{1}(:,:,1));

function coimage = my_local(image, psf, res)
    
    coimage = imresize(image,res,"bilinear");
    ims = size(coimage);
    
    for ii = 1:ims(3)
        coimage_t = normxcorr2e(psf,coimage(:,:,ii));
        coimage(:,:,ii) = imregionalmax(coimage_t);
        
    end
end

function I = normxcorr2e(template, im)

  pad = floor(size(template)./2);
  center = size(im);

  I = normxcorr2(template, im);
  I = I(pad(1)+1:pad(1)+center(1), ...
        pad(2)+1: pad(2)+center(2));
  
  threshold = 0.8;
  I = I.*imbinarize(I,threshold);

end

function parsave(fname, mattracking)
  save(fname, 'mattracking');
end


