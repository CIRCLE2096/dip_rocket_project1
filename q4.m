clear all;
matout = load('matout.mat').MatOut;
matout_pala_local = load('matout_pala_local.mat').MatOut;
matout_pala = load('matout_pala.mat').MatOut;
refimage = load('ref.mat').refImgData;

% image = uint8(image);
% figure
% imshow(image,[]);
% figure
% image_pala_local = matout_pala_local.^IntPower;
% imshow(image_pala_local,[]);


%% adjust position due to different localization methods
% matsize = size(matout);
% matout_match_position = zeros(matsize);
% proj = 5;
% for i = 1:matsize(1)
%     for j = 1:matsize(2)
%         ii = i+proj;
%         jj = j+proj;
%         if ii >=1 && ii<=matsize(1) && jj >=1 && jj<=matsize(2) 
%             matout_match_position(i,j) = matout(ii,jj);
%         end
%     end
% end

%%

image = adjust(matout);
image_pala_local = adjust(matout_pala_local);
image_pala = adjust(matout_pala);



figure(1);imshow(image);
figure(2);imshow(image_pala_local);
figure(3);imshow(image_pala);



psnr1 = psnr(image,refimage);
psnr2 = psnr(image_pala_local,refimage);
psnr3 = psnr(image_pala,refimage);

disp([psnr1,psnr2,psnr3]);

ssim1 = ssim(image,refimage);
ssim2 = ssim(image_pala_local,refimage);
ssim3 = ssim(image_pala,refimage);

disp([ssim1,ssim2,ssim3]);
%%
% IntPower = 1/3;
% image = matout.^IntPower;
% image = normalize(image);
% 
% image_pala_local = matout_pala_local.^IntPower;
% image_pala_local = normalize(image_pala_local);
% 
% image_pala = matout_pala.^IntPower;
% image_pala = normalize(image_pala);
% 
% 
% 
% figure
% imshow(image);
% figure
% imshow(image_pala_local);
% figure
% imshow(image_pala);

%%
function image_nor = normalize(a)
    maxi = max(max(a));
    mini = min(min(a));
    image_nor = uint8(double((a-mini).*255./(maxi-mini)));
end

%%
function imagedata = adjust(matout)
    f = figure('visible','off');
    load referenceImg
    clear MatOut
    addpath(['..\PALA_addons']);

    clf,set(gcf,'Position',[652 393 941 585]);
    IntPower = 1/3;
    im=imagesc(llx,llz,matout.^IntPower);axis image

    title('ULM intensity display')
    colormap(gca,gray(128))
    clbar = colorbar;caxis(caxis*.8)  % add saturation in image
    clbar.Label.String = 'number of counts';
    clbar.TickLabels = round(clbar.Ticks.^(1/IntPower),1);xlabel('\lambda');ylabel('\lambda')
    ca = gca;ca.Position = [.05 .05 .8 .9];
    WriteTif(im.CData,ca.Colormap,['temp.tif'],'caxis',caxis,'Overwrite',1)

    imagedata = imread('temp.tif');
    imagedata = imagedata(:,:,1);
%     figure(2);imshow(imagedata);
%     set(f, 'visible', 'on'); 
end
