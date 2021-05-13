function batch_batchMergeTileBuffer(outdir,quadList)

%% Blend tile edges
% get list of tiles recursively where the parent dir matches an item in the
% cell attay tileList
% Usage: batch_batchMergeTileBuffer('/path/to/tile/dirs',{'14_51_1_1','14_51_1_2'})

f=dir([outdir,'/*']);
tilef={};


for i=1:length(f)
    if (exist([outdir,'/',f(i).name],'dir') & length(f(i).name)>=5)
        q=dir([outdir,'/',f(i).name,'/*_2m.mat']);
        for j=1:length(q)
            qn=strrep(q(j).name,'_2m.mat','');
            if strmatch(qn,quadList)
                qp=[outdir,'/',f(i).name,'/',q(j).name];
                fprintf('Found quad %s\n',qp)
                if ~isdir(qp)
                    tilef = [tilef, qp];
                end
            end
        end
    end
end

batchMergeTileBuffer(tilef);

end