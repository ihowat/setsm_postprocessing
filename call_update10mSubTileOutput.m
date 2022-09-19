function call_update10mSubTileOutput(subtilepath)

% subtilepath is a filepath to the folder containing subtiles

fprintf('Updating tile %s\n',subtilepath);
subtilelist=dir([subtilepath,'/*10m*']);
sts = {subtilelist.name};
sts = cellfun( @(x) [subtilepath,'/',x], sts, 'uniformoutput', false);
update10mSubTileOutput(sts);