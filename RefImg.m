load referenceImg
addpath(['..\PALA_addons']);

figure(1);clf,set(gcf,'Position',[652 393 941 585]);
IntPower = 1/3;
im=imagesc(llx,llz,MatOut.^IntPower);axis image

title('ULM intensity display')
colormap(gca,gray(128))
clbar = colorbar;caxis(caxis*.8)  % add saturation in image
clbar.Label.String = 'number of counts';
clbar.TickLabels = round(clbar.Ticks.^(1/IntPower),1);xlabel('\lambda');ylabel('\lambda')
ca = gca;ca.Position = [.05 .05 .8 .9];
WriteTif(im.CData,ca.Colormap,['referenceImg.tif'],'caxis',caxis,'Overwrite',1)

refImgData = imread('referenceImg.tif');
refImgData = refImgData(:,:,1);
figure(2);imshow(refImgData);
save('ref.mat','refImgData');
