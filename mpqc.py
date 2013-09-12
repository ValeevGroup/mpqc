#!/usr/bin/env python
#
# Python program to run MPI programs.
#

import math
import sys
import tempfile
import os
from os.path import expandvars
import subprocess
from subprocess import Popen
import re
from string import Template

import traceback

import optparse
formatter = optparse.IndentedHelpFormatter(max_help_position=40)
parser = optparse.OptionParser(formatter=formatter)

env = {}
# default environment variables
# env["I_MPI_DEVICE"] = "sock"
# env["I_MPI_DEVICE"] = "rdma"
# env["I_MPI_WAIT_MODE"] = "enable"
# env["I_MPI_NETMASK"] = "ib0"
# env["I_MPI_PIN"] = "1"
# env["I_MPI_DEBUG"] = "4"
# env["CUDA_CACHE_DISABLE"] = "1"
# env["CUDA_CACHE_PATH"] = "/tmp/$USER-ptx-cache"

env["MV2_ENABLE_AFFINITY"] = "0"
env["PATH"] = expandvars("$PATH")
env["LD_LIBRARY_PATH"] = expandvars("$LD_LIBRARY_PATH")

# set defaults here
parser.set_defaults(cmd=expandvars("$PWD/src/bin/mpqc/mpqc"))
#parser.set_defaults(wdir=expandvars("$SCRATCH"))

#parser.set_defaults(mpi="mpich")
#parser.set_defaults(rsh="ssh") # enable to use rsh to set up job

def main():

    MPI.classes = [MPICH, OpenMPI]

    # usage and description
    parser.set_usage("%prog [options] input...")
    parser.set_description("Submit an MPI/SGE program.")

    # parse options
    add_options(parser)

    (options, args) = init(sys.argv, parser)

    # input files have to be specified
    if not args:
        parser.error("input(s) must be specified.")
        sys.exit(1)
        
    # submit or job jobs
    for a in args:
        try:
            job = MPQC(a, options)
            (cmd, script) = options.scheduler.get_script(job, options)

            if options.noop:
                print script
            else:
                fh = opentemp(prefix="%s." % (job.name), suffix="."+options.scheduler.key)
                fh.write(script)
                fh.close()
                p = Popen(cmd + " " + fh.name, shell=True)
                sts = os.waitpid(p.pid, 0)
                #os.system(cmd + " " + script) 


        except Exception, e:
            traceback.print_exc()

    sys.exit(0)



def add_options(parser):
    # common options
    parser.add_option("-N", "--nodes", type="int", help="number of nodes", metavar="X")
    parser.add_option("-p", "--ppn", type="int", help="processors per node", metavar="X")
    parser.add_option("-n",   "--np", type="int", help="number of processes", metavar="X")
    parser.add_option("-t", "--threads", type="int", help="number of threads", metavar="X")
    parser.add_option("", "--gpus", type="int", help="number of gpus", metavar="X")
    parser.add_option("", "--stack",type="int", help="stacksize")
    parser.add_option("-d", "--wdir", help="working directory", metavar="X")
    parser.add_option("-o", "--output", help="output", metavar="X")
    parser.add_option("-e", "--error", help="error", metavar="X")
    parser.add_option("-E", "--env", action="append", help="environment", metavar="X")

    parser.add_option("-z", "--noop", action="store_true", help="no op")
    parser.add_option("-s", "--sh", help="shell", default="/bin/bash")
    parser.add_option("-r", "--rsh", help="remote shell")

    parser.add_option("-P", "--path", action="append", help="path")
    parser.add_option("-L", "--libpath", action="append", help="Library path")

    parser.add_option("", "--mpi", help="MPI command")

    parser.add_option("-S", "--scheduler",
                      choices=["sge", "pbs"],
                      help="Scheduler")

    parser.add_option("-w", "--wtime", help="wall clock time", metavar="X")
    parser.add_option("-q", "--queue", help="queue")
    parser.add_option("-A", "--account", help="account")
    parser.add_option("-J", "--name", help="job name")
    
    for s in [ SGE, PBS ]:
        s.add_options(parser)

    parser.add_option("", "--valgrind", action="store_true", help="run valgrind")
    parser.add_option("", "--operf", action="store_true", help="run operf")
    parser.add_option("", "--opgprof", action="store_true", help="run opgprof")
    parser.add_option("-g", "--debug", action="store_true", help="enable debug")
    parser.add_option("-G", "--debugger",
                      choices=["gdb", "idb", "tv", "dbx"],
                      default="gdb",
                      help="debugger:gdb, idb, totalview, dbx")
    # parser.add_option("", "--gcmd", help="debugger command")    

    parser.add_option("-X", "--cmd",
                      help="executable", metavar="X")


def get_global_env():
    global env
    return env
    
def init(argv, options):

    arglen = len(argv)
    if "--" in argv: arglen = sys.argv.index("--") + 1
    # debugger args = argv[arglen:]
    (options, args) = parser.parse_args(argv[1:arglen])

    scheduler = None
    mpi = None
    env = {}

    for (k,v) in get_global_env().items():
	env[k] = v

    path = (options.path or [])
    path.append(env.get("PATH", ""))
    path = filter(None, path)
    if path: env["PATH"] = ":".join(path)

    path = (options.libpath or [])
    path.append(env.get("LD_LIBRARY_PATH", ""))
    path = filter(None, path)
    if path: env["LD_LIBRARY_PATH"] = ":".join(path)

    if options.stack:
        env["OMP_STACK_SIZE"] = options.stack*1024
        env["KMP_STACKSIZE"] = options.stack*1024
        env["XLMSOPTS"] = "stack="+str(options.stack*1024)
        env["MP_STACK_SIZE"] =  options.stack*1024
        
    if options.threads:
	env["OMP_NUM_THREADS"] = options.threads

    # if options.bg_flat_profile: env["FLAT_PROFILE"] = "yes"
    # if options.xt_symm_heap: env["XT_SYMMETRIC_HEAP"] = options.xt_symm_heap

    # if options.hostlist: options.hostlist = options.hostlist.split(',')
    options.hostlist = None
    options.hostfile = None
    
    if options.hostlist and not options.nodes:
	options.nodes = len(options.hostlist)

    #options.pbs = options.pbs or []
    #options.sge = options.sge or []

    if not options.np:
	if options.nodes and options.ppn:
	    options.np = options.nodes*options.ppn
	elif options.nodes: options.np = options.nodes
	elif options.ppn: options.np = options.ppn

    options.serial = not (options.np or options.hostfile)
    if not (options.mpi or options.serial):
	options.mpi = MPI.cmd

    if options.env:
        for kv in options.env:
            (k,v) = kv.split("=", 2)
            env[k] = v


    if options.scheduler:
        scheduler = None
        for s in [ SGE, PBS ]:
            if s.key == options.scheduler:
                scheduler = s()
                break
        if not scheduler:
            raise RuntimeError("Unknown scheduler '%s'" % options.scheduler)
    else:
	scheduler = Shell()


    if options.mpi:
        mpi = None
        options.mpicmd = options.mpi
        for C in MPI.classes:
            if C.matches(options.mpicmd):
                mpi = C(options)
                break
        assert mpi, "Undefined MPI: %s" % options.mpi

    if options.debug:
        if options.debugger == "gdb": options.gcmd = "gdb"
        # options.interactive = True
        # if not options.gcmd:
        #     if options.debugger == "gdb": options.gcmd = "gdb"
        #     elif options.debugger == "idb": options.gcmd = "idb"
        #     elif options.debugger == "tv": options.gcmd = "totalview"
        #     elif options.debugger == "dbx": options.gcmd = "dbx"

    options.mpi = mpi
    options.scheduler = scheduler
    options.env = env

    return (options, args)


# Scheduler class
class Scheduler:
    cmd = None
    arg = None
    
    def getOptionGroup(parser):
        return None 
    getOptionGroup = staticmethod(getOptionGroup)

    def __init__(self):
        global options
        #if options.qcmd: self.cmd = options.qcmd
        #if options.qarg: self.arg = options.qarg

    def get_script(self, job, options, stdout=None, stderr=None):

        script = ""

	if not job.hostlist:
	    job.hostlist = options.hostlist or []
	if not job.hostfile:
	    job.hostfile = options.hostfile

	# #script += "set -x\n"
	# script += "echo $(hostname)\n"
	# script += "mkdir -p %s\n" % job.wdir

	hosts = None
	hostfile = None

	if job.hostfile:
	    hosts = "$(cat %s | awk '{print $1}' | grep -v '^#')" % job.hostfile
	if job.hostlist:
	    hosts = " ".join(job.hostlist)

        if hosts and options.rsh:
            script += "for host in %s; do\n" % (hosts)
            script += "%s $host \"%s\";\n" % (options.rsh, job.getProlog("\\\\\\"))
            script += "done\n"
	else:
	    script += "%s\n" % (job.getProlog("\\\\\\"))
        script += "\n"


	if job.hostlist:
	    job.hostfile = job.wdir + "/machinefile"
	    script += "echo \"%s\" > %s\n" % ("\n".join(job.hostlist), job.hostfile)
	    job.hostlist = ()

	if job.hostfile:
	    job.hostfile = os.path.abspath(job.hostfile)
	    hosts = "$(cat %s | awk '{print $1}' | grep -v '^#')" % job.hostfile
	    hostfile = os.path.join(job.wdir, "machines")

        script += "\n"


        if not hostfile:
	    script += "# no hostfile/hostlist\n"
            script += job.getProlog()
        else:
	    script += "hostfile=%s\n" % hostfile
            hostfile = "$hostfile"
	    filt = "cat"
	    if options.nodes: filt = "sort | uniq"
            if options.ppn: filt += "| sed '%s'" % ((options.ppn-1)*"p;") 
	    script += "cat %s | %s > %s\n" % (job.hostfile, filt, hostfile)
	    job.hostfile = hostfile

        script += "\n"

	# # generate input files
        # for inp in job.files():
        #     script += "echo '%s' | cat > %s\n" % (inp.string, inp.name)
	# #     if hosts and options.rsh:
	# # 	script += "for host in %s; do\n" % (hosts)
	# # 	script += "echo '%s' | %s $host \"cat > %s\";\n" % (inp.string, options.rsh, inp.name)
	# # 	script += "done\n"
	# #     else:
        # script += "\n"

        if len(filter(None, (options.valgrind, options.operf, options.opgprof))) > 1:
            error("Only one of valgrind, operf, or opgprof can be selected")

        script += "cd %s\n" % (job.wdir)

	if options.serial:
	    cmd = job.cmd
	    args = job.args
	    if options.debug:
                #cmd = options.gcmd + " "  + " ".join(options.gargs) + " " + cmd
		cmd = "%s --args %s" % (options.gcmd, cmd)
	    if options.valgrind:
		cmd = "valgrind -v " + cmd
        else:
            initCmd = options.mpi.getInitCmd(job)
            if initCmd: script += "%s\n\n" % (initCmd)
            cmd = options.mpi.getCmd()
            args = options.mpi.getArgs(job, options)

	if options.operf:
	    cmd = "operf --callgraph " + cmd

        if stdout and (stdout == stderr):
            stderr = "&1"

        if args: cmd += " %s" % (args)
        if stdout: cmd += " >%s" % stdout
        if stderr: cmd += " 2>%s" % stderr
        script += cmd + "\n"

	epilog = job.getEpilog("\\\\\\")
	if epilog:
	    if hosts and options.rsh:
		script += "for host in %s; do\n" % (hosts)
		script += "%s $host \"%s\";\n" % (options.rsh, epilog)
		script += "done\n"
	    else:
		script += epilog

	script += "\n"

        return script


# Shell scheduler
class Shell(Scheduler):
    cmd = "/bin/sh"
    key = "sh"
    jobid = "$$"
    name = "Shell"

    def get_script(self, job, options):
        
        script = "#/bin/sh\n"
        if options.stack: script += "ulimit -s %i\n" % (options.stack)

        for (k,v) in job.env.items():
            script += "export %s=\"%s\"\n" %(k,v)
        script += "\n"

        script += Scheduler.get_script(self, job, options)
        
        return (self.cmd, script)    


# SGE scheduler
class SGE(Scheduler):
    key = "sge"
    name = "SGE"
    cmd = "qsub"
    jobid = "$JOB_ID"
    name = "SGE"

    def __init__(self):
        pass
     	#if options.qarg: self.cmd += (" " + options.qarg)
        #if options.interactive: self.cmd += " -I"
        
    @staticmethod
    def add_options(parser):
        parser.add_option("", "--sge", action="append",
                         help="SGE directives", metavar="X")

    def get_script(self, job, options):

        script = "#$ -S /bin/sh\n"
        
        job.hostfile = "/$TMPDIR/machines"
        
        for a in options.sge or []:
            script += "#$ %s\n" % a
        if options.account:
            script += "#$ -A %s\n" % (options.account)
        if options.queue:
            script += "#$ -q %s\n" % (options.queue)

	if options.mpi:
	    n = options.np
	    pe = "mpi"
	    pe = { MPICH.key : "mpich",
		   OpenMPI.key : "orte" }.get(options.mpi.key)
	    if options.nodes: (pe, n) = ("mpichthr", 8*options.nodes)
	    script += "#$ -pe %s %i\n" % (pe, n)
	#script += "#$ -pe %s %i\n" % ("mpi", options.nodes or options.np or 1)

        #if not options.interactive:
        script += "#$ -o %s\n" % (job.output)
        script += "#$ -e %s\n" % (job.error)
        if job.output == job.error: script += "#$ -j yes\n"

        # if options.vars: script += "#$ -V\n"
	# script += "\n"
        
	#script += "\n".join(["#$ %s" % directive for directive in options.sge or []])
	script += "\n"
	script += "\n"

        # job environment
        for (name, value) in job.env.items():
            script += "export %s=\"%s\"\n" % (name, value)
	script += "\n"

        script += Scheduler.get_script(self, job, options)
        
        return (self.cmd, script)


# PBS scheduler
class PBS(Scheduler):
    key = "pbs"
    name = "PBS"
    cmd = "qsub"
    jobid = "$PBS_JOBID"
    name = "PBS"

    def __init__(self):
        pass
     	#if options.qarg: self.cmd += (" " + options.qarg)
        #if options.interactive: self.cmd += " -I"
        
    @staticmethod
    def add_options(parser):
        parser.add_option("", "--pbs", action="append",
                         help="PBS directives", metavar="X")

    def get_script(self, job, options):
        
        job.hostfile = "/$PBS_NODEFILE"

        script = "#PBS -S /bin/sh\n"

        script = "#PBS -N %s\n" % job.name
        
        for directive in options.pbs or []:
            script += "#PBS %s\n" % directive
        if options.account:
            script += "#PBS -A %s\n" % (options.account)
        if options.queue:
            script += "#PBS -q %s\n" % (options.queue)

	if options.nodes:
	    script += "#PBS -l nodes=%i\n" % options.nodes

        if options.wtime:
            script += "#PBS -l walltime=%s\n" % options.wtime

	script += "#PBS -o %s.oe%s\n" % (job.name, self.jobid)
	script += "#PBS -joe\n"
        # #if not options.interactive:
        # script += "#PBS -o %s\n" % (job.output)
        # if job.output == job.error:
        #     script += "#PBS -joe\n"
        # else:
        #     script += "#PBS -e %s\n" % (job.error)
        
	script += "\n"
	script += "\n"

        # job environment
        for (name, value) in job.env.items():
            script += "export %s=\"%s\"\n" % (name, value)
	script += "\n"

        script += Scheduler.get_script(self, job, options,
                                       stdout=job.output,
                                       stderr=job.error)
        
        return (self.cmd, script)
    

class MPI:
    name = "MPI"
    cmd = "mpiexec"
    classes = ()

    def getOptionGroup(parser):
        return None
    getOptionGroup = staticmethod(getOptionGroup)
    
    def matches(cmd, regex = None):
	if not regex: return False
	p = Popen(cmd, stdout = subprocess.PIPE, 
        stderr = subprocess.STDOUT, shell = True)
	p.wait()
	while True:
	    line = p.stdout.read()
	    if not line: break
	    if re.search(regex, line): return True
	return False
    matches = staticmethod(matches)

    def __init__(self, options):
        if options.mpi: self.cmd = options.mpi
        pass

    def getInitCmd(self, job):
        return None

    def getCmd(self):
        return self.cmd

    def getArgs(self):
        return None

    def add_options(mpi_classes, parser):
	help = "MPI: "
	help += ", ".join(["%s - %s" % (c.key, c.name) for c in mpi_classes])
	parser.add_option("-m", "--mpi", dest="mpi", help=help)
	parser.add_option("", "--mpicmd", help="mpi command")
	parser.add_option("-f", "--hostfile", help="hostfile")
	parser.add_option("-l", "--hostlist", help="hostlist")
    add_options = staticmethod(add_options)

    
# BlueGene job
class BlueGene(MPI):
    key = "bg"
    cmd = "mpirun"
    name = "Blue Gene"
    
    def getOptionGroup(parser):
        group = optparse.OptionGroup(parser, "Blue Gene options")
        group.add_option("", "--bg-mode", choices=["co", "vn", "smp"],
                         help="Blue Gene mode: co, vn, or smp", metavar="X")
        group.add_option("", "--bg-partition", type="string",
                          help="Blue Gene partition", metavar="X")
        group.add_option("", "--bg-host", type="string",
                          help="Blue Gene host", metavar="X")
        group.add_option("", "--bg-timeout", type="string",
                          help="Blue Gene timeout", metavar="X")
        group.add_option("", "--bg-flat-profile", action="store_true",
                         default=False, help ="FLAT_PROFILE")
        return group
    getOptionGroup = staticmethod(getOptionGroup)

    def __init__(self):
        global options
        
        if options.mpicmd: self.cmd = options.mpicmd
        
        if not options.ppn:
            options.ppn = 1
            if options.bg_mode == "smp":
                options.ppn = 4
            if options.bg_mode == "vn":
                options.ppn = 2
        if not options.nodes:
            options.nodes = 1
            if options.np:
                options.nodes = math.ceil(float(options.np)/options.ppn)
        if not options.np: options.np = options.nodes * options.ppn

#         if mpi == BlueGene: script += "#@job_type = bluegene\n"
        
#         if mpi == BlueGene:
#             if mpi.partition: script += "#@bg_partition = %s\n" % (mpi.partition)
#             if options.np: script += "#@bg_size = %s\n" % (options.np)

                                                
    def getArgs(self, job):
        global options
        args = ""
        
        if options.bg_mode:
            args += " -mode %s" % (options.bg_mode.upper())
        if options.bg_partition:
            args += " -partition %s" % (options.bg_partition)
        if options.bg_host:
            args += " -host %s" % (options.bg_host)
        if options.bg_timeout:
            args += " -timeout %s" % (options.bg_timeout)
       
        args += " -np %s" % (options.np)
        
        if job.wdir: args += " -cwd %s" % (job.wdir)
        for (k, v) in job.env.items():
            args += " -env \"%s=%s\"" % (k, v)
        if job.cmd: args += " %s" % (job.cmd)
        if job.args: args += " %s" % (job.args)

        return args.strip()


# MPICH job
class MPICH(MPI):
    key = "mpich"
    name = "MPICH"
    
    def getOptionGroup(parser):
        group = optparse.OptionGroup(parser, "MPICH options")
        group.add_option("", "--mpich-mpdboot", type="string", default="mpdboot",
                         help="mpdboot", metavar="X")
        return group
    getOptionGroup = staticmethod(getOptionGroup)

    def matches(cmd):
	pattern = "mpiexec \[-h or -help or --help\]"
        if MPI.matches(cmd, pattern): return True
        pattern = "HYDRA"
	return MPI.matches("%s -info" % cmd, pattern)
    matches = staticmethod(matches)

    def getInitCmd(self, job):
        return ""
	# global options
	# cmd = ""
	# if options.mpich_mpdboot and job.hostfile:
        #     mpdboot = "%s -f %s" % (options.mpich_mpdboot, job.hostfile)
        #     cmd += "hosts=$(sort %s | uniq | wc -l)\n" % job.hostfile
        #     cmd += "%s -n $(expr $hosts + 1) || " % mpdboot
        #     cmd += "%s -n $(expr $hosts + 0)\n" % mpdboot
	# return cmd

    def getArgs(self, job, options):
        args = ""

        # global args

	#args += " -demux select"
        #args += "-iface ib0"
	if job.hostlist: args += " -hosts %s" %(",".join(job.hostlist))
        if job.hostfile: args += " -machinefile %s" %(job.hostfile)
	#if options.ppn: args += " -perhost %s" %(options.ppn)
        if options.debug and options.debugger == "gdb": args += " -gdb"

	np = options.np
	if not np and job.hostfile:
	    np = 0
	    f = open(job.hostfile)
	    for line in f: np += bool(line.strip())
	    f.close()

        # local args
        if np: args += " -np %s" % np
        #if job.hostlist: args += " -host %s" % job.hostlist[0]
        if job.wdir: args += " -wdir %s" % (job.wdir)
        envlist = ",".join(job.env.keys())
        if envlist: args +=  " -envlist " + envlist
	if options.valgrind: args += " valgrind "
        if job.cmd: args += " %s" % (job.cmd)
        if job.args: args += " %s" % (job.args)

        return args.strip()


# OpenMPI job
class OpenMPI(MPI):
    key = "ompi"
    name = "OpenMPI"

    def matches(cmd):
	return MPI.matches(cmd, "(OpenRTE|Open MPI)")
    matches = staticmethod(matches)

    def getArgs(self, job, options):
        args = ""
        if options.hostfile: args += " -hostfile %s" % (options.hostfile)
	if options.hostlist: args += " -host %s" % (",".join(options.hostlist))
        if options.np: args += " -np %s" % options.np
	if options.debug: args += " --debug"
	#if options.debugger: args += " --debugger %s" % (options.debugger)
	#if options.verbose: args += " -v -d"
        if job.wdir: args += " -wdir %s" % (job.wdir)
        for (k, v) in job.env.items():
            args += " -x %s=\"%s\"" % (k, v)
        if job.cmd: args += " %s" % (job.cmd)
        if job.args: args += " %s" % (job.args)

        return args.strip()


# XT job
class XT(MPI):
    key = "xt"
    mpicmd = "aprun"
    name = "XT"
    def __init__(self):
	pass
         # if options.np: options.pbs += ["-l size=%i" % (options.np)]

    def getCmd(self):
        global options
        cmd = self.mpicmd
        if options.debug and options.debugger == "tv":
            cmd = options.gcmd + " " + cmd + " -a"
        return cmd
            
    def getOptionGroup(parser):
        group = optparse.OptionGroup(parser, "XT options")
        group.add_option("", "--xt-symm-heap", type="string",
                         help="XT_SYMMETRIC_HEAP", metavar="X")
        return group
    getOptionGroup = staticmethod(getOptionGroup)

    def getArgs(self, job):
        global options
        args = ""
        if options.np: args += " -n %s" % (options.np)
        if options.ppn: args += " -N %s" % (options.ppn)
        if options.threads: args += " -d %s" % (options.threads)
        if job.cmd: args += " %s" % (job.cmd)
        if job.args: args += " %s" % (job.args)
        return args.strip()

class Job:
    hostfile = None
    hostlist = None
    name = None
    wdir = None
    output = None
    error = None
    env = {}
    cmd = None
    args = None

class MPQC(Job):
    cmd = "mpqc"

    def getOptionGroup(parser):
        return None
    getOptionGroup = staticmethod(getOptionGroup)

    def __init__(self, input, options):

        if not os.path.isabs(input):
            input = os.path.abspath(input)

        self.env = options.env

        if options.cmd: self.cmd = options.cmd
        self.cmd += " %s " % input

        self.name = os.path.basename(input.strip())
        # strip off suffix
        if self.name.lower().endswith(".in"):
            self.name = self.name[:len(self.name)-3]

        # working directory
        jobid = options.scheduler.jobid
	if not options.wdir:
	    self.wdir = os.path.dirname(os.path.abspath(input))
        else:
            self.wdir = options.wdir # + ".%s" % (jobid)
        #if jobid: self.wdir = self.wdir + ".%s" % (jobid)
        # if options.wdir:
        #     self.wdir = os.path.join(options.wdir, self.wdir)

        # output
        #if not options.interactive:
        #if options.output: self.output = options.output
        self.output = os.path.join(os.getcwd(), self.name + ".out")
        #if options.error: self.error = options.error
        self.error = os.path.join(os.getcwd(), self.name + ".out")
                                    
        self.cmd = os.path.expandvars(self.cmd)
        self.args = ""

    # def files(self):

    #     class Input:
    #         name = ""
    #         string = ""

    #     f = open(self.input, "r")
    #     inp = Input
    #     inp.name = self.name
    #     inp.string = f.read()
    #     # for line in f:
    #     #     if line.strip() and line.strip()[0] == "!": continue
    #     #     inp.string += "%s" % line #(line[:-1].replace("'", ""))
    #     f.close()
    #     return [ inp ]


    # write working dir
    def getProlog(self, escape="\\"):
        global options
        wdir = self.wdir
        prolog = ""

        #prolog += "rm -fr %s/%s.*\n" % (self.wdir, self.name)
        #prolog += "rm -fr %s/\n" % (self.wdir)
        prolog += "mkdir -p %s\n" % (self.wdir)
        prolog += "\n"

        return prolog

    # write ending of the script
    def getEpilog(self, escape="\\"):
        global options
        script = ""
        # script += "rm -f %s/*.F*\n" % (self.wdir)
        # script += "rm -f %s/*.h5\n" % (self.wdir)
        # if not options.keep: script += "rm -fr %s/\n" % (self.wdir)
        return script


def warning(message):
    print "Warning: %s" % (message)

    
def error(message):
    print "Error: %s" % (message)
    sys.exit(1)


def opentemp(prefix=None, suffix=None):
    (fd, fn) = tempfile.mkstemp(prefix=prefix, suffix=suffix)
    os.close(fd)
    return open(fn, "w")


if __name__ == "__main__":
    main()
