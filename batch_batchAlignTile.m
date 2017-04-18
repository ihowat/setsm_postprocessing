function batch_batchAlignTile(outdir,tileList)
%% Align unregistered tiles

% find unregistered tiles recursively where the parent dir matches an item in the
% cell attay tileList and intersect/remove with registered list
% Usage: batch_batchAlignTile('/path/to/tile/dirs',{'14_51','14_52'})

f=dir([outdir,'/*']);
m={};
mreg={};

for i=1:length(f);
    f1=f(i);
    if exist([outdir,'/',f1.name],'dir')
        if strmatch(f1.name,tileList)
            g=dir([outdir,'/',f1.name,'/*2m_dem.mat']);
            gp=[outdir,'/',f1.name,'/',g.name];
            if ~isdir(gp)
                m = [m, gp];
            end

            g2=dir([outdir,'/',f1.name,'/*2m_reg_dem.mat']);
            gp2=[outdir,'/',f1.name,'/',g2.name];
            if ~isdir(gp2)
                mreg = [mreg, gp2];
            end
        end
    end
end

[~,IA]=intersect(m,strrep(mreg,'_reg',''));
m(IA)=[];

batchAlignTile(m,mreg);