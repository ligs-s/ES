#!/usr/bin/env python

from ROOT import gROOT
import os, sys

def get_list(listname):
    f = open(listname, 'r')
    filelist = []
    for line in f:
        if line.startswith('#') or line=='\n':
            continue
        line = line.strip()
        filelist.append(line)
    return filelist

def make_job(jobdir, jobname, cmd):

    if not os.path.exists(jobdir):
        os.system('mkdir -p %s' % jobdir)
    jobfile = jobdir+'/'+jobname
    f = open(jobfile, 'w')
    f.write(cmd)
    f.close()
    os.system("chmod 755 %s" % jobfile)
    sub_cmd = 'bsub -W 60 -R rhel60 -o %s/%s.out -e %s/%s.err %s' % (jobdir, jobname, jobdir, jobname, jobfile)

    #sub_cmd = 'bsub -W 60 -R rhel60 "%s"' % cmd

    print sub_cmd
    os.system(sub_cmd)
    return
    
        
if __name__=='__main__':
    """ """
    filelist = get_list('list.txt')
    jobdir = './job'
    workdir = os.getcwd()
    job_idx = 0
    for line in filelist:   
        cmd = "#!/bin/bash\n"
        cmd += "cd %s\n" % workdir
        #cmd += """root -l -b -q 'dat2root.C("%s", "wave0.dat")'\n""" % line
        cmd += "./compiled/dat2root %s\n" % line
        cmd += """root -l -b -q 'process_wf.C("root/%s/wave0.dat.root", "root/%s/charge.root")'\n""" % (line, line)

        #cmd = "cd %s;" % workdir
        #cmd += "./compiled/dat2root %s;" % line
        #cmd += "root -l -b -q 'process_wf.C(\"root/%s/wave0.dat.root\", \"root/%s/charge.root\")';" % (line, line)

        make_job(jobdir, "job%d.sh" % job_idx, cmd)
        job_idx += 1
