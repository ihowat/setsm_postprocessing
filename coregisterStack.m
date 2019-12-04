function out=coregisterStack(x,y,z,mask)


N = size(z,3); % # of strips

%      maskFile=strrep(subsetFiles{file_n},'subset.mat','subset_mask.mat');
%     if exist(maskFile,'file')
%         mask=load(maskFile);
%         mask=repmat(mask.A,1,1,N);
%         z(~mask)=NaN;
%     end
%

if N < 2
    out=[];
    fprintf('less than two DEMs in the stack, skipping\n')
    return
end

I = nchoosek(1:N,2); % indices of unique pairs

i=I(:,1); % first strip in pair
j=I(:,2); % second strip in pair

Npairs=size(I,1); % # of pairs

dz = nan(Npairs,1); % initialize pair offset vector
dx = nan(Npairs,1); % initialize pair offset vector
dy = nan(Npairs,1); % initialize pair offset vector

dze = nan(Npairs,1); % initialize pair offset sigma vectors
dxe = nan(Npairs,1); % initialize pair offset sigma vectors
dye = nan(Npairs,1); % initialize pair offset sigma vectors

mean_dz_uncoreg =  nan(Npairs,1);
median_dz_uncoreg =  nan(Npairs,1);
sigma_dz_uncoreg =  nan(Npairs,1);

mean_dz_coreg =  nan(Npairs,1);
median_dz_coreg =  nan(Npairs,1);
sigma_dz_coreg =  nan(Npairs,1);

% Pair coregistration loop
for pair_n=1:Npairs
    
   % fprintf('i:%d, j:%d,pair %d of %d\n',i(pair_n),j(pair_n),pair_n,Npairs)
    
    % get overlap stats before coregistration
    p = z(:,:,i(pair_n)) - z(:,:,j(pair_n));
    p(~mask) = NaN;
    mean_dz_uncoreg(pair_n) = nanmean(p(:));
    
    if isnan(mean_dz_uncoreg(pair_n))
        continue
    end
    
    median_dz_uncoreg(pair_n) = nanmedian(p(:));
    sigma_dz_uncoreg(pair_n) = nanstd(p(:));
    
    % coregister pair
    [zj,p,perr] = coregisterdems(x,y,z(:,:,i(pair_n)),x,y,z(:,:,j(pair_n)),mask);
    
    dz(pair_n) = p(1);
    dx(pair_n) = p(2);
    dy(pair_n) = p(3);
    
    dze(pair_n) = perr(1);
    dxe(pair_n) = perr(2);
    dye(pair_n) = perr(3);
    
    
    if ~isnan(dz(pair_n))
        p = z(:,:,i(pair_n)) - zj;
        p(~mask) = NaN;
        mean_dz_coreg(pair_n) = nanmean(p(:));
        median_dz_coreg(pair_n) = nanmedian(p(:));
        sigma_dz_coreg(pair_n) = nanstd(p(:));
    else
        mean_dz_coreg(pair_n) =NaN;
        median_dz_coreg(pair_n) = NaN;
        sigma_dz_coreg(pair_n) = NaN;
    end
end





out.i=i;
out.j=j;
out.dz = dz;
out.dx = dx;
out.dy = dy;

out.dze = dze;
out.dxe = dxe;
out.dye = dye;

out.mean_dz_uncoreg=mean_dz_uncoreg;
out.median_dz_uncoreg=median_dz_uncoreg;
out.sigma_dz_uncoreg=sigma_dz_uncoreg;

out.mean_dz_coreg=mean_dz_coreg;
out.median_dz_coreg=median_dz_coreg;
out.sigma_dz_coreg=sigma_dz_coreg;


