function [nc]=mkclrmap(lowc,midc,hic)

%
% make fancy colormap
%
% John Bradford gave me this and god knows what
% it does.
%
% example:
%
% lowc = [0 0 1];
% midc = [.99 .99 .99];
% hic = [.5 .12 0];
% cmap = mkclrmap(lowc,midc,hic);
% colormap(cmap)

nc = zeros(3,256);

nc(1,1:128) = [lowc(1):(midc(1)-lowc(1))/127:midc(1)];
nc(1,128:256) = [midc(1):(hic(1)-midc(1))/128:hic(1)];
nc(2,1:128) = [lowc(2):(midc(2)-lowc(2))/127:midc(2)];
nc(2,128:256) = [midc(2):(hic(2)-midc(2))/128:hic(2)];
nc(3,1:128) = [lowc(3):(midc(3)-lowc(3))/127:midc(3)];
nc(3,128:256) = [midc(3):(hic(3)-midc(3))/128:hic(3)];
nc = nc.';

end
