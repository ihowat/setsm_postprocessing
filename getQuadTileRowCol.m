function [quad_row, quad_col] = getQuadTileRowCol(tileFileName)
[tileFileDir, tileFileName_noext, tileFileExt] = fileparts(tileFileName);
[prefix, super_row, super_col, quad_row, quad_col, suffix] = parseTileNameAll(tileFileName_noext, true);
