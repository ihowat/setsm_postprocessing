function batch_batchMergeTileBuffer(outdir,tileList)

%% Blend tile edges
% get list of tiles recursively where the parent dir matches an item in the
% cell attay tileList
% Usage: batch_batchMergeTileBuffer('/path/to/tile/dirs',{'14_51','14_52'})

f=dir([outdir,'/*']);
tilef={};

for i=1:length(f)
    f1=f(i);
    if exist([outdir,'/',f1.name],'dir')
        if strmatch(f1.name,tileList)
            g=dir([outdir,'/',f1.name,'/*2m_reg_dem.mat']);
            gp=[outdir,'/',f1.name,'/',g.name];
            if ~isdir(gp)
                tilef = [tilef, gp];
            end
        end
    end
end

batchMergeTileBuffer(tilef);

end

