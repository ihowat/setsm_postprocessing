function combineDatabase4_triple(dbase1,dbase2,dbase3,dbase_out)

fprintf("Loading dbase1 as meta1: %s\n", dbase1);
meta1 = load(dbase1);
fprintf("Loading dbase2 as meta2: %s\n", dbase2);
meta2 = load(dbase2);
fprintf("Loading dbase3 as meta3: %s\n", dbase3);
meta3 = load(dbase3);

fprintf("Appending meta2 and meta3 to meta1\n");
flds = fields(meta1);
i=1;
for i = 1:length(flds)
    fld = flds{i};
    eval(['meta1.',fld,'= [meta1.',fld,',meta2.',fld,',meta3.',fld,'];']);
end

fprintf("Clearing meta2\n");
clear meta2;
fprintf("Clearing meta3\n");
clear meta3;
fprintf("Saving combined meta1 database: %s\n", dbase_out);
save(dbase_out,'-struct','meta1','-v7.3');
fprintf("Clearing meta1\n");
clear meta1;

rlist1 = strrep(dbase1, '.mat', '_reproject_list.txt');
rlist2 = strrep(dbase2, '.mat', '_reproject_list.txt');
rlist3 = strrep(dbase3, '.mat', '_reproject_list.txt');
rlist_out = strrep(dbase_out, '.mat', '_reproject_list.txt');

fprintf("Combining reproject lists to: %s\n", rlist_out);
system(sprintf('cat %s %s %s > %s', rlist1, rlist2, rlist3, rlist_out));

fprintf("Done!\n")
