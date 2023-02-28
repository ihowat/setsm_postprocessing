function combineDatabase4_quad(dbase1,dbase2,dbase3,dbase4,dbase_out)

fprintf("Loading dbase1 as meta1: %s\n", dbase1);
meta1 = load(dbase1);
flds = fields(meta1);

rlist1 = strrep(dbase1, '.mat', '_reproject_list.txt');
rlist_out = strrep(dbase_out, '.mat', '_reproject_list.txt');
fprintf("Appending %s > %s\n", rlist1, rlist_out);
system(sprintf('cat %s > %s', rlist1, rlist_out));

dbase_list = {dbase2, dbase3, dbase4};
for i=1:length(dbase_list)
    dbase = dbase_list{i};
    fprintf("Loading dbase as meta2: %s\n", dbase);
    meta2 = load(dbase);
    fprintf("Appending meta2 to meta1\n");
    j=1;
    for j = 1:length(flds)
        fld = flds{j};
        eval(['meta1.',fld,'= [meta1.',fld,',meta2.',fld,'];']);
    end
    fprintf("Clearing meta2\n");
    clear meta2;
    rlist = strrep(dbase, '.mat', '_reproject_list.txt');
    fprintf("Appending %s >> %s\n", rlist, rlist_out);
    system(sprintf('cat %s >> %s', rlist, rlist_out));
end

fprintf("Saving combined meta1 database: %s\n", dbase_out);
save(dbase_out,'-struct','meta1','-v7.3');
fprintf("Clearing meta1\n");
clear meta1;

fprintf("Done!\n")
