function supername = get10mTileName(tileFileName)
[tileFileDir, tileFileName_noext, tileFileExt] = fileparts(tileFileName);
[parsed_tilename, supername, quadname, suffix] = parseTileNameWithSuffix(tileFileName_noext);
