
import argparse
import configparser
import os
import math
import platform
import subprocess
import sys
import textwrap
import time
import warnings
from collections import defaultdict
from datetime import datetime
import pathlib


# Script paths
SCRIPT_FILE = os.path.abspath(os.path.realpath(__file__))
SCRIPT_FNAME = os.path.basename(SCRIPT_FILE)
SCRIPT_NAME, SCRIPT_EXT = os.path.splitext(SCRIPT_FNAME)
SCRIPT_DIR = os.path.dirname(SCRIPT_FILE)
HOME_DIR = os.path.expanduser('~/')

# Paths relative to this script
JOBSCRIPT_INIT_ENV_SCRIPT = os.path.join(SCRIPT_DIR, 'jobscript_init_env.sh')
DEFAULT_JOBSCRIPT_FILE = os.path.join(SCRIPT_DIR, 'jobscript_<scheduler>.sh')
DEFAULT_CONFIG_FILE = os.path.join(SCRIPT_DIR, 'jobscript_config.ini')

# Global vars
DEFAULT_TASK_BUNDLE_DIR = os.path.join(os.path.expanduser('~'), 'scratch', 'task_bundles')
SCHEDULER_PBS = 'pbs'
SCHEDULER_SLURM = 'slurm'
_sched_name_checkcmd_dict = {
    SCHEDULER_PBS: 'pbsnodes',
    SCHEDULER_SLURM: 'sinfo',
}
SCHEDULER_AVAIL_LIST = None
_sched_opt_flags_list = None


class UnsupportedMethodError(Exception):
    def __init__(self, msg=""):
        super(Exception, self).__init__(msg)


def wrap_multiline_str(text, width=float('inf')):
    """Format a multiline string, preserving indicated line breaks.

    Wraps the `text` (a string) so every line is at most `width` characters long.
    Common leading whitespace from every line in `text` is removed.
    Literal '\n' are considered line breaks, and area treated as such in wrapping.

    Args:
        text (str): A multiline string to be wrapped.

    Returns:
        str: The wrapped string.

    Example:
        animal_a = "Cats"
        animal_b = "Dogs"
        text = wrap_multiline_rfstr(
            rf\"""
            Cats and dogs are the most popular pets in the world.
            \n  1) {animal_a} are more independent and are generally
            cheaper and less demanding pets.
            \n  2) {animal_b} are loyal and obedient but require more
            attention and exercise, including regular walks.
            \""", width=40
        )
        >>> print(text)
        Cats and dogs are the most popular pets
        in the world.
          1) Cats are more independent and are
        generally cheaper and less demanding
        pets.
          2) Dogs are loyal and obedient but
        require more attention and exercise,
        including regular walks.
    """
    s_in = textwrap.dedent(text.strip('\n'))

    p_in = [p for p in s_in.split(r'\n')]
    p_out = [textwrap.fill(p, width=width) for p in p_in]

    return '\n'.join(p_out)


def get_available_sched():
    global SCHEDULER_AVAIL_LIST

    locate_cmd = 'where' if platform.system() == 'Windows' else 'which'

    sched_avail_list = []

    for sched_name in sorted(list(_sched_name_checkcmd_dict.keys())):
        sched_checkcmd = _sched_name_checkcmd_dict[sched_name]
        test_cmd_args = [locate_cmd, sched_checkcmd]
        try:
            proc = subprocess.Popen(test_cmd_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            if proc.returncode == 0:
                sched_avail_list.append(sched_name)
        except OSError:
            pass

    if not sched_avail_list:
        sched_avail_list = [None]

    SCHEDULER_AVAIL_LIST = sched_avail_list.copy()

    return sched_avail_list


def get_config_value(config_dict, key, default_val=None):
    if key not in config_dict:
        val = default_val
    else:
        val = config_dict[key]
        if val == '':
            val = default_val
    return val


def argparse_add_job_scheduler_group(
        parser,
        config_file=None,
        config_group=None,
):
    global _sched_opt_flags_list
    add_to_sched_opt_flags_list = (_sched_opt_flags_list is None)
    if add_to_sched_opt_flags_list:
        _sched_opt_flags_list = []

    if config_file is None and DEFAULT_CONFIG_FILE is not None:
        if os.path.isfile(DEFAULT_CONFIG_FILE):
            config_file = DEFAULT_CONFIG_FILE
    
    available_sched_list = get_available_sched()
    jobscript = None
    jobname_prefix = None
    ncpus = None
    mem = None
    walltime = None
    queue = None
    node_list = None
    tasks_per_job = None
    tasks_per_job_mode = None
    task_bundle_dir = None

    provided_argstr_set = set()
    for token_idx, token in enumerate(sys.argv):
        if token.startswith('-'):
            provided_argstr_set.add(token.split('=')[0])
        if token.startswith('--job-config-file'):
            if token.startswith('--job-config-file='):
                config_file = token.split('=')[0]
            else:
                config_file = sys.argv[token_idx + 1]

    if '--job-config-on' in provided_argstr_set:
        use_config = True
    elif '--job-config-off' in provided_argstr_set:
        use_config = False
    elif '--job-script' in provided_argstr_set and '--job-config-file' not in provided_argstr_set:
        use_config = False
    else:
        use_config = True

    config_was_used = False
    if use_config and config_file is not None:
        config = configparser.ConfigParser()
        config.read(config_file)

        if config_group is None or config_group not in config:
            config_group = 'DEFAULT'
        if config_group in config:
            config_dict = config[config_group]
            config_was_used = True

            jobscript = get_config_value(config_dict, 'jobscript')
            jobname_prefix = get_config_value(config_dict, 'jobname-prefix')
            ncpus = get_config_value(config_dict, 'ncpus')
            mem = get_config_value(config_dict, 'mem')
            walltime = get_config_value(config_dict, 'walltime')
            queue = get_config_value(config_dict, 'queue')
            if 'nodelist' in config_dict:
                node_list = get_config_value(config_dict, 'nodelist')
                if node_list is not None:
                    node_list = [name.strip() for name in node_list.split(',')]
            tasks_per_job = get_config_value(config_dict, 'tasks-per-job')
            tasks_per_job_mode = get_config_value(config_dict, 'tasks-per-job-mode')
            task_bundle_dir = get_config_value(config_dict, 'task-bundle-dir')

    if jobscript is None:
        jobscript = DEFAULT_JOBSCRIPT_FILE
    if node_list is None:
        node_list = [None]
    if tasks_per_job is None:
        tasks_per_job = 1
    if tasks_per_job_mode is None:
        tasks_per_job_mode = 'serial'
    if task_bundle_dir is None:
        task_bundle_dir = DEFAULT_TASK_BUNDLE_DIR

    group = parser.add_argument_group(
        title='Job Scheduler',
        description=textwrap.fill(textwrap.dedent("""
            Job scheduler arguments.
        """))
    )

    if config_was_used:
        config_group_note = r'\n' + wrap_multiline_str(rf"""
            >>> Note that the config group '{config_group}' was used to populate
            job argument default settings <<<
        """)
    else:
        config_group_note = ''
    group.add_argument(
        '--job-config-file',
        type=str,
        default=(config_file if config_file is not None else DEFAULT_CONFIG_FILE),
        help=wrap_multiline_str(rf"""
            Configuration file used to populate script-specific default job options.
            {config_group_note}
            \nProvide the --job-config-off option to disable automatic use of the
            default config file.
            \nIf the --job-script option is provided, the default config file
            will not be used unless the --job-config-file option is explicitly provided
            or the --job-config-on option is provided.
        """)
    )
    group.add_argument(
        '--job-config-off',
        action='store_true',
        help=wrap_multiline_str("""
            Do not use the --job-config-file, leaving job settings to be
            set from the --job-script or provided through relevant job setting
            script arguments.
        """)
    )
    group.add_argument(
        '--job-config-on',
        action='store_true',
        help=wrap_multiline_str("""
            Force useage of the --job-config-file in situations where
            the default config file would not be automatically used.
        """)
    )

    sched_name_check_cmd_table = r'\n   '.join(
        ['{}: {}'.format(sched_name, check_cmd)
         for sched_name, check_cmd in _sched_name_checkcmd_dict.items()]
    )
    group.add_argument(
        '--job-scheduler',
        type=str,
        choices=available_sched_list,
        default=None,
        help=wrap_multiline_str(rf"""
            Submit processing tasks to the designated job scheduler.
            \nChoices are schedulers identified as available by checking
            for existance of a relevant command using the following lookup table.
            \n   {sched_name_check_cmd_table}
            \n
        """)
    )
    # for sched in available_sched_list:
    for sched in sorted(list(_sched_name_checkcmd_dict.keys())):
        if add_to_sched_opt_flags_list:
            _sched_opt_flags_list.append(sched)
        group.add_argument(
            f'--{sched}',
            action='store_true',
            help=f"Shortcut for --job-scheduler='{sched}'."
        )
    group.add_argument(
        '--job-script',
        type=str,
        default=jobscript,
        help="Jobscript to use in scheduler job submission."
    )
    group.add_argument(
        '--job-name-prefix',
        type=str,
        default=jobname_prefix,
        help="String prefix affixed to job names for sumitted jobs."
    )
    group.add_argument(
        '--job-log-dir',
        type=pathlib.Path,
        default=pathlib.Path.home(),
        help=wrap_multiline_str("""
            Directory to write slurm logs (combined stdout & stderr).
        """)
    )

    group.add_argument(
        '--job-ncpus',
        type=int,
        default=ncpus,
        help="Number of CPUs requested per job."
    )
    group.add_argument(
        '--job-mem',
        type=int,
        default=mem,
        help="Memory requested per job (GB)."
    )
    group.add_argument(
        '--job-walltime',
        type=int,
        default=walltime,
        help="Wallclock time per job (hours)."
    )

    group.add_argument(
        '--job-queue',
        type=str,
        default=queue,
        help="Job queue to submit jobs to."
    )
    group.add_argument(
        '--job-node-list',
        type=str,
        nargs='+',
        default=node_list,
        help="Specific cluster node names to submit jobs to."
    )
    group.add_argument(
        '--job-hold',
        type=str,
        choices=['on', 'off'],
        default='off',
        help=wrap_multiline_str("""
            Submit jobs in held (H) state.
            Jobs will need to be released manually after submission
            before they will run.
        """)
    )

    # TODO: Implement --max-jobs argument,
    # -t    which will become --max-processes in later multiprocessing addition.

    group.add_argument(
        '--tasks-per-job',
        type=int,
        default=tasks_per_job,
        help="Run multiple tasks within each scheduler job."
    )
    group.add_argument(
        '--tasks-per-job-mode',
        type=str,
        choices=['serial', 'parallel'],
        default=tasks_per_job_mode,
        help="Mode for running multiple tasks within a scheduler job."
    )
    group.add_argument(
        '--task-bundle-dir',
        type=str,
        default=task_bundle_dir,
        help=wrap_multiline_str("""
            Directory where task bundle textfiles will be automatically
            created if --tasks-per-job > 1.
        """)
    )


def adjust_args(script_args, arg_parser):
    sched_avail_list = SCHEDULER_AVAIL_LIST
    if sched_avail_list is None:
        sched_avail_list = get_available_sched()

    if _sched_opt_flags_list is not None:
        sched_flags_provided = [opt for opt in _sched_opt_flags_list if hasattr(script_args, opt) and getattr(script_args, opt) is True]
        num_sched_flags_provided = len(sched_flags_provided)
        if num_sched_flags_provided == 0:
            pass
        elif num_sched_flags_provided > 1:
            arg_parser.error("Only one job scheduler shortcut flag can be provided")
        else:
            sched_opt = sched_flags_provided[0]
            if script_args.job_scheduler is not None:
                arg_parser.error(f"Job scheduler shortcut flag --{sched_opt} and --job-scheduler option are mutually exclusive")
            if sched_opt not in sched_avail_list:
                warnings.warn(wrap_multiline_str(f"""
                    Job scheduler selected by shortcut flag --{sched_opt}
                    was not identified as an available scheduler on this system.
                    Identified schedulers are the following: {','.join(sched_avail_list)}
                """))
            script_args.job_scheduler = sched_opt

    jobscript_path = script_args.job_script
    if jobscript_path is not None:
        if script_args.job_scheduler is not None:
            jobscript_path = jobscript_path.replace('<scheduler>', script_args.job_scheduler)
        script_args.job_script = jobscript_path

    jobname_prefix = script_args.job_name_prefix
    if jobname_prefix is not None:
        jobname_prefix = jobname_prefix.replace('<resolution>', str(script_args.resolution))
        script_args.job_name_prefix = jobname_prefix


def create_dirs(script_args, arg_parser):
    if script_args.job_scheduler is not None and script_args.tasks_per_job > 1:
        if not os.path.isdir(script_args.task_bundle_dir):
            print(f"Creating directory where task bundle files will be created: {script_args.task_bundle_dir}")
            os.makedirs(script_args.task_bundle_dir)


def verify_args(script_args, arg_parser):
    if script_args.job_scheduler is not None:
        if script_args.job_script is None:
            arg_parser.error("--job-script file must be specified when using --scheduler option")
        elif not os.path.isfile(script_args.job_script):
            arg_parser.error(f"--job-script file does not exist: {script_args.job_script}")
        if script_args.tasks_per_job > 1:
            if not os.path.isdir(script_args.task_bundle_dir):
                arg_parser.error(f"--task-bundle-dir directory does not exist: {script_args.task_bundle_dir}")


def escape_problem_jobsubmit_chars(
        str_item,
        escape_single_quotes=True,
        escape_double_quotes=True,
        escape_comma=True,
        escape_space=True,
):
    if escape_single_quotes:
        str_item = str_item.replace("'", "\\'")
    if escape_double_quotes:
        str_item = str_item.replace('"', '\\"')
    if escape_comma:
        str_item = str_item.replace(',', '@COMMA@')
    if escape_space:
        str_item = str_item.replace(' ', '@SPACE@')
    return str_item


def get_jobsubmit_cmd(
        script_args,
        job_name=None,
        node_name=None,
        envvar_list_str=None,
        hold_jobs=False
):
    scheduler = script_args.job_scheduler
    if scheduler is None:
        return None

    resource_arg_values = [
        script_args.job_ncpus,
        script_args.job_mem,
        script_args.job_walltime,
        node_name
    ]
    resource_args_provided = any([val is not None for val in resource_arg_values])

    if scheduler == SCHEDULER_PBS:

        resource_list_str = None
        if resource_args_provided:
            resource_list = []
            if script_args.job_ncpus is not None or node_name is not None:
                resource_list.append('nodes={}{}'.format(
                    node_name if node_name is not None else 1,
                    ':ppn={}'.format(script_args.job_ncpus) if script_args.job_ncpus is not None else ''
                ))
            if script_args.job_mem is not None:
                resource_list.append('mem={}gb'.format(script_args.job_mem))
            if script_args.job_walltime is not None:
                resource_list.append('walltime={}:00:00'.format(script_args.job_walltime))
            resource_list_str = ','.join(resource_list)

        cmd_str = ' '.join([
            'qsub',
            '-N {}'.format(job_name) if job_name is not None else '',
            '-q {}'.format(script_args.job_queue) if script_args.job_queue is not None else '',
            '-l {}'.format(resource_list_str) if resource_list_str is not None else '',
            '-v "{}"'.format(envvar_list_str) if envvar_list_str is not None else '',
            '-h' if script_args.job_hold == 'on' else '',
            script_args.job_script
        ])

    elif scheduler == SCHEDULER_SLURM:

        cmd_str = ' '.join([
            'sbatch',
            '--job-name {}'.format(job_name) if job_name is not None else '',
            '--partition {}'.format(script_args.job_queue) if script_args.job_queue is not None else '',
            '--ntasks {}'.format(script_args.job_ncpus) if script_args.job_ncpus is not None else '',
            '--mem {}G'.format(script_args.job_mem) if script_args.job_mem is not None else '',
            '--time {}:00:00'.format(script_args.job_walltime) if script_args.job_walltime is not None else '',
            '--export "{}"'.format(envvar_list_str) if envvar_list_str is not None else '',
            f'-o {script_args.job_log_dir / "%x.o%j"}',
            f'-e {script_args.job_log_dir / "%x.o%j"}',
            script_args.job_script
        ])

    else:
        raise UnsupportedMethodError("Unsupported job scheduler name: '{}'".format(job_name))

    return cmd_str


def matlab_cmdstr_to_shell_cmdstr(cmdstr):
    return """matlab -nojvm -nodisplay -nosplash -r "try; {}; catch e; disp(getReport(e)); exit(1); end; exit(0);" """.format(cmdstr)


def yield_loop(iterable):
    """Generator to unendingly yield items from a non-exhaustive iterable.

    Items are yielded from the iterable through a for-loop. Once all items
    have been yielded, the for-loop repeats to yield additional items.

    Args:
        iterable: An non-exhaustive iterable collection containing items
          to be retrieved.

    Yields:
        An item in the `iterable` collection.

    Examples:
        >>> item_gen = yield_loop([0, 1])
        >>> for i in range(5):
        >>>     print(next(item_gen))
        0
        1
        0
        1
        0
    """
    while True:
        for item in iterable:
            yield item


def get_left_zero_pad_fmtstr(max_num, min_digits=3):
    return '{:0>' + str(max(min_digits, len(str(max_num)))) + '}'
            
            
class JobHandler(object):
    
    def __init__(
            self,
            script_args,
            num_tasks,
            num_tasks_is_estimate=False,
            init_env_requests=None,
    ):
        self.script_args = script_args
        self.init_env_requests = '' if init_env_requests is None else init_env_requests
        self.submit_to_scheduler = (self.script_args.job_scheduler is not None)

        self.bundle_tasks = None
        self.last_bundle_file = None
        self.last_bundle_timestamp = None
        self._current_bundle_fp = None
        self.reset_bundle_attr()

        self.num_tasks_is_estimate = num_tasks_is_estimate
        self.num_tasks = None
        self.num_jobs = None
        self.task_num = None
        self.job_num = None
        self.job_num_fmtstr = None
        self.reset_counts(num_tasks)

        self.gen_job_node = yield_loop(self.script_args.job_node_list)

    def reset_bundle_attr(self):
        self.bundle_tasks = (self.submit_to_scheduler and self.script_args.tasks_per_job > 1)
        self.last_bundle_file = None
        self.last_bundle_timestamp = None
        if self._current_bundle_fp is not None:
            self._current_bundle_fp.close()
            self._current_bundle_fp = None
        
    def reset_counts(
            self,
            num_tasks,
            task_num=0,
            job_num=0
    ):
        self.num_tasks = num_tasks
        if self.bundle_tasks:
            self.num_jobs = int(math.ceil(num_tasks / float(self.script_args.tasks_per_job)))
        else:
            self.num_jobs = num_tasks
        self.task_num = task_num
        self.job_num = job_num
        self.job_num_fmtstr = get_left_zero_pad_fmtstr(self.num_jobs)

    def add_task_cmd(self, task_cmd, job_id, force_submit=False):
        self.task_num += 1
        
        if self.submit_to_scheduler and self.bundle_tasks:
            if self._current_bundle_fp is None:
                self.job_num += 1
                self.create_bundle_file()
            self._current_bundle_fp.write(task_cmd+'\n')
            
            if (    self.task_num % self.script_args.tasks_per_job == 0
                or  (self.task_num == self.num_tasks and not self.num_tasks_is_estimate)
                or  force_submit):
                self._current_bundle_fp.close()
                self._current_bundle_fp = None
                submit_cmd = self.get_jobsubmit_cmd()
            else:
                submit_cmd = None
            
        else:
            self.job_num += 1
            if self.submit_to_scheduler:
                submit_cmd = self.get_jobsubmit_cmd(task_cmd, job_id)
            else:
                submit_cmd = task_cmd
            
        return submit_cmd

    def get_last_bundle_submit_cmd(self):
        submit_cmd = None
        if self._current_bundle_fp is not None:
            self._current_bundle_fp.close()
            self._current_bundle_fp = None
            submit_cmd = self.get_jobsubmit_cmd()
        return submit_cmd
        
    def get_jobsubmit_cmd(self, task_cmd=None, job_id=None):
        envvar_dict = {
            'JOBSCRIPT_INIT_ENV_SCRIPT_PATH': JOBSCRIPT_INIT_ENV_SCRIPT,
            'JOBSCRIPT_INIT_ENV_REQUESTS': self.init_env_requests,
        }
        
        if task_cmd is None:
            envvar_dict['JOBSCRIPT_TASK_BUNDLE_FILE'] = self.last_bundle_file
            envvar_dict['JOBSCRIPT_TASK_BUNDLE_MODE'] = self.script_args.tasks_per_job_mode
        else:
            envvar_dict['JOBSCRIPT_TASK_CMD'] = task_cmd

        escape_kwargs = dict()
        if self.script_args.job_scheduler == SCHEDULER_SLURM:
            escape_kwargs['escape_single_quotes'] = False
            # escape_kwargs['escape_double_quotes'] = False

        envvar_list_str = ','.join(
            ['{}={}'.format(key, escape_problem_jobsubmit_chars(val, **escape_kwargs))
             for key, val in envvar_dict.items()]
        )

        submit_cmd = get_jobsubmit_cmd(
            self.script_args,
            self.get_job_name(job_id),
            node_name=next(self.gen_job_node),
            envvar_list_str=envvar_list_str,
        )
        
        return submit_cmd
    
    def create_bundle_file(self):
        if self.script_args.job_name_prefix is None:
            bundle_fname_prefix = os.path.splitext(os.path.basename(self.script_args.job_script))[0]
        else:
            bundle_fname_prefix = self.script_args.job_name_prefix

        while True:
            if self.last_bundle_timestamp is None:
                timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
                check_for_existing_file = True
            else:
                timestamp = self.last_bundle_timestamp
                check_for_existing_file = False

            bundle_fname = "{}_{}_{}.txt".format(
                bundle_fname_prefix,
                timestamp,
                self.job_num_fmtstr.format(self.job_num)
            )
            bundle_file = os.path.join(self.script_args.task_bundle_dir, bundle_fname)

            if check_for_existing_file:
                assert self.job_num == 1
                if os.path.isfile(bundle_file):
                    time.sleep(1)
                    continue
                else:
                    self.last_bundle_timestamp = timestamp

            break

        self.last_bundle_file = bundle_file
        self._current_bundle_fp = open(bundle_file, 'w')
        
    def get_job_name(self, job_id=None):
        if self.script_args.job_name_prefix is None:
            job_name = None
        else:
            job_name = "{}{}".format(
                self.script_args.job_name_prefix,
                '_{}'.format(job_id) if job_id is not None else self.job_num_fmtstr.format(self.job_num)
            )
        return job_name
