
// Hi! Swift's documentation is awful, so here are a few important notes:
// 1. If your task is going to be executing a bash script, you can pass around any necessary
//    arguments as string-type variables and list those in the bashcoaster call like you would
//    when calling a normal bash script.
// 2. The advantage of using "file"-type objects instead of string-type objects to hold your
//    file path variables is that whenever these files need to be read or written, Swift will
//    do so by {copying those files to be read to} or {writing those files to be written into}
//    Swift's working directory (AKA "swiftwork" in PGC lingo).
//    So what's the big advantage of letting file I/O be managed by Swift?
//    You have the option to let I/O for certain files happen on the compute nodes!
//    Read more on Swift's "Collective Data Management" (CDM) documentation here:
//    http://swift-lang.org/guides/release-0.94/userguide/userguide.pdf
// 3. CDM provides 3 different options for file I/O:
//    1) DEFAULT:
//       - Just use file staging as provided by Swift. Identical to behavior if no CDM file is given.
//    2) DIRECT:
//       - Allows for direct I/O to the parallel FS without staging.
//    3) LOCAL:
//       - Allows for client-directed input copy to the compute node.
// 4. IMPORTANT!!
//    The "@" prefix is usually required for file-type variable references when I/O is expected as
//    a result of the method being called on them.
//    The Swift documentation would indicate this prefix as being part of calls to Swift's standard
//    functions (like 'length'), but use of the prefix in that context has been deprecated.


type file;

type logstruct {
    file out;
    file err;
};


app (file outlog) bashlocal (string cmd) {
    bashlocal "-c" cmd stdout=@outlog;
}

app (logstruct log) bashcoaster (string jobscript, string task_str) {
    bashcoaster jobscript task_str stdout=filename(log.out) stderr=filename(log.err);
}


// Custom swift script command line args
string tasklist_file = arg("tasklist_file");
string jobscript = arg("jobscript");
string task_description = arg("task_description");
string logdir = arg("logdir");

// Get list of tasks for swift jobs to manage
string load_tasklist_cmd = sprintf("cat %s", tasklist_file);
string tasklist_arr[] = readData(bashlocal(load_tasklist_cmd));


if (length(tasklist_arr) == 0) {
    tracef("%s\n", "Tasklist is empty");
} else {

    tracef("Submitting %i %s\n", length(tasklist_arr), task_description);
    foreach task_str, task_num in tasklist_arr {
        tracef("%s,%i\n", task_str, task_num);

        //string task_args[] = strsplit(task_str, "\\%");
        //string task_logpath = task_args[1];
        string tilename = task_str;
        string task_logpath = strcat(logdir, "/", tilename);

        logstruct log <simple_mapper; prefix=strcat(task_logpath,".")>;

        log = bashcoaster(jobscript, task_str);
    }
}
