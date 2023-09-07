function run_mosaicSubTiles(...
    tileName,quadrant,tileDefFile,tileParamListFile,...
    subTileDir,resolution,outMatFile,...
    projection,version,exportTif...
)

if exportTif
    outRasterType = 'full-COG';
else
    outRasterType = 'browse-COG';
end

try
    if ~isempty(version)
        version = strrep(version, ',', '|');
    end
    if isempty(quadrant)
        [x0,x1,y0,y1] = getTileExtents(tileName,tileDefFile,'tileParamListFile',tileParamListFile);
        mosaicSubTiles(subTileDir,resolution,outMatFile,...
            'projection',projection,...
            'version',version,...
            'outRasterType',outRasterType,...
            'extent',[x0,x1,y0,y1]...
        )
    else
        [x0,x1,y0,y1] = getTileExtents(tileName,tileDefFile,'quadrant',quadrant,'tileParamListFile',tileParamListFile);
        mosaicSubTiles(subTileDir,resolution,outMatFile,...
            'projection',projection,...
            'version',version,...
            'quadrant',quadrant,...
            'outRasterType',outRasterType,...
            'extent',[x0,x1,y0,y1]...
        )
    end
catch e
    disp(getReport(e)); exit(1)
end

disp("End of run_mosaicSubTiles.m");
