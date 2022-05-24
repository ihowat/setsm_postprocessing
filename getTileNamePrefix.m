function prefix = getTileNamePrefix(tileFileFull)
[tileFileDir, tileFileName_noext, tileFileExt] = fileparts(tileFileFull);
[prefix, supertile, quadtile, suffix] = parseTileNameParts(tileFileName_noext, true);
if isempty(prefix)
    prefix = 'none';
end
