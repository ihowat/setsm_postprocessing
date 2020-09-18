function projstr = getTileProjection(tileDefFile)

tileDefs=load(tileDefFile);

if ~isfield(tileDefs,'projstr')
    error("'projstr' field must be added to tile definition structure: %s", tileDefFile);
end

projstr = tileDefs.projstr;
