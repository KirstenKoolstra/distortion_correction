function [decomposed_field] = decompose_field(spherHarm_res,B0temp)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

mask = ones(size(B0temp));

fieldMap = B0temp; 
fieldMap(round(end/2),round(end/2))= 0;
% figure; imagesc(fieldMap);

mask(fieldMap==0)=NaN;
flattenField = fieldMap(mask == 1);

% spherHarm_res = spherHarm_res.*repmat(mask,[1 1 size(spherHarm_res,3)]);

% figure;
% imagesc(spherHarm_res(:,:,4)); colorbar;

for index=1:size(spherHarm_res,3)
    temp = spherHarm_res(:,:,index);
    spherHarm_res_mask(:,index) = temp(mask==1);
end

fitResult = pinv(spherHarm_res_mask)*flattenField;

decomposed_field = zeros(size(spherHarm_res,1),size(spherHarm_res,2));
for index = 1:length(fitResult)
    decomposed_field = decomposed_field + spherHarm_res(:,:,index)*fitResult(index);
end
decomposed_field(isnan(decomposed_field)) = 0;

% figure;
% subplot(1,2,1); imagesc(B0temp,[-800 800]);
% subplot(1,2,2); imagesc(decomposed_field, [-800 800]);

end

