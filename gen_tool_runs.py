import argparse
import sys
import os
from configparser import ConfigParser
from collections import defaultdict


#def gen_leafcutter_run(outfile, name, grp1_name, grp2_name, grp1list, grp2list, outDir):
def gen_leafcutter_run(comp_dict, cmdline, args): 

    outDir = args.output
    outname = args.outname
    nthreads = args.nthreads
    bam_dir = args.bam_dir
    annot_db = args.annot_db
    nreads = args.nreads 
    intronlen = args.intronlen
    minsamplesintron = args.min_samples_x_intron
    minsamplesgroup = args.min_samples_x_group
    mincov = args.min_coverage


    if not os.path.exists('%s/conf_juncfiles' % outDir):
        os.makedirs('%s/conf_juncfiles' % outDir)

    if not os.path.exists('%s/diff_files' % outDir):
        os.makedirs('%s/diff_files' % outDir)

    with open('%s/run_%s.sh' % (outDir, outname), 'w+') as outfile:

        outfile.write('#!/bin/bash\n')
        outfile.write('#\n')
        outfile.write('#This file has been generated using:\n# \t %s \n \n' % cmdline)
        outfile.write('\n')
        outfile.write('\n')
        outfile.write('leafcutter=%s\n' % args.bindir)
        outfile.write('\n')
        outfile.write('mkdir -p %s/logs\n' % outDir)
        outfile.write('mkdir -p %s/out_gz\n' % outDir)
        outfile.write('mkdir -p %s/output_ds\n' % outDir)
        outfile.write('\n')
        outfile.write('\n')
       
        for name, dd in comp_dict.items():
            cond1 = dd['grp1']
            cond2 = dd['grp2']
            ln =    dd['cmd']            

            minsizegroup = min (len(cond1[1]), len(cond2[1]))
            minsamplesgroup = min(minsizegroup, args.min_samples_x_group)
            minsamplesintron = min(minsizegroup, args.min_samples_x_intron) 

            jfile_name = '%s/conf_juncfiles/%s_juncfile.txt' % (outDir, name)
            diff_intron_file = '%s/diff_files/%s.txt' %(outDir, name)

            with open(jfile_name, 'w+') as jfile, open(diff_intron_file, 'w+') as di_file:
                for nn in cond1[1]:
                    jfile.write('%s/%s.bam.junc\n' % (bam_dir, nn) )
                    di_file.write('%s.bam\t%s\n' %(nn, cond1[0]))

                for nn in cond2[1]:
                    jfile.write('%s/%s.bam.junc\n' % (bam_dir, nn) )
                    di_file.write('%s.bam\t%s\n' %(nn, cond2[0]))

            outgz = 'out_gz/%s' % (name)
            logfile = '%s/logs/%s' %(outDir, name)

            outfile.write('#%s\n' % ln)                       
            cmd = "python $leafcutter/clustering/leafcutter_cluster.py -j %s -m %s -o %s -l %s > %s" % (jfile_name, nreads, outgz, intronlen, logfile)
            outfile.write('%s\n' % cmd)
        
#            cmd = "$leafcutter/scripts/leafcutter_ds.R -p %s -e %s -o %s/output_ds/%s  %s_perind_numers.counts.gz %s >> %s" %(nthreads, annot_db, outDir, name, outgz, diff_intron_file, logfile)
            cmd = "$leafcutter/scripts/leafcutter_ds.R -p %s -e %s -o %s/output_ds/%s --min_samples_per_intron %s --min_samples_per_group %s --min_coverage %s  %s_perind_numers.counts.gz %s >> %s" %(nthreads, annot_db,
                                                                                                                                 outDir, name, minsamplesintron, minsamplesgroup, mincov, outgz, diff_intron_file, logfile)
            outfile.write('%s\n\n' % cmd)


def gen_rMATS_run(comp_dict, cmdline, args): 
    
    outDir = args.output
    outname = args.outname
    nthreads = args.nthreads
    bam_dir = args.bam_dir
    threshold = args.threshold
    readtype = args.readtype
    libtype = args.libtype
    readlen = args.readlen
    annot_db = args.annot_db
    
    if not os.path.exists('%s/run_files' % outDir):
        os.makedirs('%s/run_files' % outDir)

    with open('%s/run_%s.sh' % (outDir, outname), 'w+') as outfile:
   
        outfile.write('#!/bin/bash\n')
        outfile.write('#\n')
        outfile.write('#This file has been generated using:\n#\t %s \n \n' % cmdline)
        outfile.write('\n')
        outfile.write('\n')
        outfile.write('rmatsdir=%s\n' % args.bindir)
        outfile.write('\n')
        outfile.write('mkdir -p %s/output_ds\n' % outDir)
        outfile.write('mkdir -p %s/logs\n' % outDir)
        outfile.write('\n')
        outfile.write('\n')

          
        for name, dd in comp_dict.items():
            cond1 = dd['grp1']
            cond2 = dd['grp2']
            ln =    dd['cmd'] 
            
            outfile.write('#%s\n' % ln)
            jfile_name1 = '%s/run_files/%s_%s.txt' % (outDir, name, cond1[0])
            jfile_name2 = '%s/run_files/%s_%s.txt' % (outDir, name, cond2[0])
   
            with open(jfile_name1, 'w+') as jfile:
                print (cond1[1])
                kk = ','.join(['%s/%s.bam' % (bam_dir, nn) for nn in cond1[1]])
                jfile.write('%s\n' %kk)
           
            with open(jfile_name2, 'w+') as jfile:
                kk = ','.join(['%s/%s.bam' % (bam_dir, nn) for nn in cond2[1]])
                jfile.write('%s\n' %kk)

            cmd = "python $rmatsdir/rmats.py --b1 %s --b2 %s --gtf %s --od %s/output_ds/%s -t %s --readLength %s --cstat %s --libType %s --nthread %s  --tstat %s --tmp /tmp/rmats_%s --variable-read-length > %s/logs/%s"  % (jfile_name1, jfile_name2, annot_db, outDir, name, readtype, readlen, threshold, libtype, nthreads, nthreads, name, outDir, name)
  
            outfile.write('%s\n' % cmd)


def gen_whippet_run(comp_dict, cmdline, args):

    outDir = args.output
    outname = args.outname
    nthreads = args.nthreads
    fa_dir = args.bam_dir
    threshold = args.threshold
    annot_db = args.annot_db
    annot_fa = args.annot_fa
    outPsi = "%s/psi" % outDir
    index = args.index
    fa_ext = args.fasta_ext

    if not os.path.exists(outPsi):
        os.makedirs(outPsi)

    samples_d = {}
    deltas = []


    with open('%s/run_%s.sh' % (outDir, outname), 'w+') as outfile:

        outfile.write('#!/bin/bash\n')
        outfile.write('#\n')
        outfile.write('#This file has been generated using:\n#\t %s \n \n' % cmdline)
        outfile.write('\n')
        outfile.write('\n')
        outfile.write('whippetdir=%s\n' % args.bindir)
        outfile.write('\n')
        outfile.write('mkdir -p %s/output_ds\n' % outDir)
        outfile.write('mkdir -p %s/logs\n' % outDir)
        outfile.write('\n')
        outfile.write('\n')

        if index is None:
            print("Not index provided, generating index")
            index = "%s/index.%s.jls" % (outDir, os.path.basename(annot_fa))
            cmd = " julia $whippetdir/whippet-index.jl --fasta %s --gtf %s --index %s " %(annot_fa, annot_db, index)
            outfile.write('%s\n' % cmd)

        for name, dd in comp_dict.items():
            cond1 = dd['grp1']
            cond2 = dd['grp2']
            ln = dd['cmd']

            outfile.write('\n# %s\n' % ln)

            for fa in cond1[1]:
                sfa = "%s/%s" %(fa_dir, fa)
                cmd = "julia $whippetdir/whippet-quant.jl %s.1.%s %s.2.%s -o %s/%s -x %s &> %s/logs/%s.log" % (sfa, fa_ext, sfa, fa_ext, outPsi, fa, index, outDir, fa)
                if not args.separate_run:
                    if fa not in samples_d:
                        samples_d[fa] = cmd
                else:
                    outfile.write('%s\n' % cmd)

            for fa in cond2[1]:
                sfa = "%s/%s" %(fa_dir, fa)
                cmd = "julia $whippetdir/whippet-quant.jl %s.1.%s %s.2.%s -o %s/%s -x %s &> %s/logs/%s.log" % (sfa, fa_ext, sfa, fa_ext, outPsi, fa, index, outDir, fa)
                if not args.separate_run:
                    if fa not in samples_d:
                        samples_d[fa] = cmd
                else:
                    outfile.write('%s\n' % cmd)



            
            grp1 = ','.join(['%s/%s.psi.gz' % (outPsi, nn) for nn in cond1[1]])
            grp2 = ','.join(['%s/%s.psi.gz' % (outPsi, nn) for nn in cond2[1]])
            if len(cond1[1])==1:
                grp1+=','
            if len(cond2[1])==1:
                grp2+=','

            psi_list = ' '.join(['%s/%s.psi.gz' % (outPsi, nn) for nn in cond1[1]])
            psi_list += ' ' + ' '.join(['%s/%s.psi.gz' % (outPsi, nn) for nn in cond2[1]])
            cmd = "mkdir -p %s/%s\n" % (outDir, name)
            cmd += "julia $whippetdir/whippet-delta.jl -a %s -b %s -o %s/%s/%s &> %s/logs/%s.log\n "  % (grp1, grp2, outDir, name, name, outDir, name)
            cmd += "ln -s %s %s/.\n" %(psi_list, name)
            if not args.separate_run:
                deltas.append(cmd)
            else:
                outfile.write('%s' % cmd)

        if not args.separate_run:
            
            for cmd in samples_d.values():
                outfile.write('%s\n' % cmd)

            for cmd in deltas:
                outfile.write('%s\n' % cmd)


def gen_SUPPA_run(comps, cmdline, args):


    # generate salmon runs
    outdir = args.output
    outname = args.outname
    nthreads = args.nthreads
    faformat = args.fa_format

    quantdir = "%s/quantification" % (outdir)

    with open('%s/run_salmon_%s.sh' % (outdir, outname), 'w+') as outsalmon:
        outsalmon.write('#!/bin/bash\n')
        outsalmon.write('#\n')
        outsalmon.write('#This file has been generated using:\n#\t %s \n \n' % cmdline)
        # outsalmon.write('. %s/activate\n' % args.bindir)
        outsalmon.write('\n')
        outsalmon.write('mkdir -p %s\n' % quantdir)

        file_dict = {}
        for name, samples in comps.items():
            for grpid in ['grp1', 'grp2']:
                for sample in samples[grpid][1]:
                    if args.paired:
                        smps = "-1 %s/%s_1.%s -2 %s/%s_2.%s" % (args.bam_dir, sample, faformat, 
                                                                args.bam_dir, sample, faformat)
                    else:
                        smps = "-r %s/%s.%s" % (args.bam_dir, sample, faformat) 

                    cmd = "%s/salmon quant -i %s -l %s %s --validateMappings -p %s -o %s/%s" % (args.salmondirbin, args.salmon_index, 
                                                                       args.libtype, smps, args.nthreads, quantdir, sample)
                    if sample not in file_dict :
                        file_dict[sample] = cmd

        # lst_quants = ""
        for sample, cmd in file_dict.items():
            outsalmon.write("%s\n" % cmd)
            # lst_quants += "%s/%s/quant.sf " % (quantdir, sample)

    scrpt_dir = os.path.dirname(os.path.realpath(__file__))
    for name, samples in comps.items():
        odir = "%s/%s" % (outdir, name)
        if not os.path.exists(odir):
            os.makedirs(odir)
        with open('%s/run_%s.sh' % (odir, outname), 'w+') as outfile :

            s_events = "%s/suppa_events/" % odir
            event_dir = "%s/events" % odir
            exec_quant = "%s/quant_for_suppa" % odir

            outfile.write('mkdir -p %s\n' % s_events)
            outfile.write('mkdir -p %s\n' % exec_quant)
            outfile.write('mkdir -p %s\n' % event_dir)
            outfile.write('mkdir -p %s/output_ds/%s/\n' % (outdir, name))

            outfile.write("##1. Extract the TPM values from the Salmon output\n")
            lst_quants  = ' '.join(["%s/%s/quant.sf" % (quantdir, x) for x in samples['grp1'][1]])
            lst_quants += ' '
            lst_quants += ' '.join(["%s/%s/quant.sf" % (quantdir, x) for x in samples['grp2'][1]])

            outfile.write("python %s/multipleFieldSelection.py -i %s -k 1 -f 4 -o %s/iso_tpm.txt\n" %(args.bindir,
                                                                                                      lst_quants,
                                                                                                      exec_quant))

            cmd = "#2. Before running SUPPA, we need to calculate the AS events on the hg19 annotation\n" \
                  "#2.1: Generate the events:\n" \
                  "python %s/suppa.py generateEvents -i %s -f ioe -o %s -e SE SS MX RI FL \n" %(args.bindir,
                                                                                         args.annot_db, s_events)
            outfile.write("%s\n" % cmd)

            cmd = "for file in $(ls %s | grep .ioe);do\n" \
                  "\techo \"Processing $file...\"\n" \
                  "\tcat %s/$file >> %s/summary.events_formatted.ioe\n done" %(s_events, s_events, s_events)
            outfile.write("%s\n" % cmd)

            outfile.write("#2.2: Put all the gtf events in the same file:\n")
            cmd = "for file in $(ls %s | grep .gtf);do\n" \
                  "\t echo \"Processing $file...\" \n" \
                  "\t cat %s/$file >> %s/summary.events_formatted.gtf\n done " %(s_events, s_events, s_events)
            outfile.write("%s\n" % cmd)

            outfile.write("#3: Run SUPPA for getting the psi values of the events:\n")
            cmd = "python %s/suppa.py psiPerEvent -i %s/summary.events_formatted.ioe " \
                  "-e %s/iso_tpm.txt -o %s/events" % (args.bindir, s_events, exec_quant, event_dir)

            outfile.write("%s\n\n" % cmd)

        # for name, samples in comps.items():
            grp1 = samples['grp1'][0]
            grp2 = samples['grp2'][0]
            smplst1 = ' '.join([x for x in samples['grp1'][1]])
            smplst2 = ' '.join([x for x in samples['grp2'][1]])

            cmd = "python %s/split_file.py %s/iso_tpm.txt -c %s -o %s/%s_iso.tpm " %(scrpt_dir, exec_quant,
                                                                                     smplst1, exec_quant, grp1)
            outfile.write("%s\n" % cmd)
            cmd = "python %s/split_file.py %s/events.psi -c %s -o %s/%s_events.psi" %(scrpt_dir, event_dir, smplst1,
                                                                                      exec_quant, grp1)
            outfile.write("%s\n" % cmd)

            cmd = "python %s/split_file.py %s/iso_tpm.txt -c %s -o %s/%s_iso.tpm " %(scrpt_dir, exec_quant,
                                                                                     smplst2, exec_quant, grp2)
            outfile.write("%s\n" % cmd)
            cmd = "python %s/split_file.py %s/events.psi -c %s -o %s/%s_events.psi" %(scrpt_dir, event_dir, smplst2,
                                                                                      exec_quant, grp2)
            outfile.write("%s\n" % cmd)

            m = max(len(samples['grp1'][1]), len(samples['grp2'][1]))
            if m>=10:
                runtype = "classical"
            else:
                runtype = "empirical"

            cmd = "python %s/suppa.py diffSplice -m %s -i %s/summary.events_formatted.ioe " \
                  "-e %s/%s_iso.tpm %s/%s_iso.tpm -p %s/%s_events.psi %s/%s_events.psi " \
                  "-o %s/output_ds/%s" % (args.bindir, runtype, s_events, exec_quant, grp1, exec_quant, grp2, exec_quant,
                                grp1, exec_quant, grp2, outdir, name)
            outfile.write("%s\n\n" % cmd)

    with open('%s/run_%s.sh' % (outdir, outname), 'w+') as outfile:
        outfile.write('#!/bin/bash\n')
        outfile.write('#\n')
        outfile.write('#This file has been generated using:\n#\t %s \n \n' % cmdline)

        for name, samples in comps.items():
            odir = "%s/%s" % (outdir, name)
            fname = '%s/run_%s.sh' % (odir, outname)

            outfile.write("#%s\n" %samples['cmd'])
            outfile.write("sh %s\n\n" % fname)





    # generate suppa runs


            
def gen_majiq_run(comps, cmdline, args):
    def print_run(outfile, method, directory, grp1, grp2, name1, name2, extra_flags=''):
        deltapsi_dir = os.path.join(outdir, '%s_%s' % (method, directory))
        cmd = 'majiq %s -o %s -grp1 %s -grp2 %s -j %d --names %s %s %s %s' % (method, deltapsi_dir, grp1, grp2, nthreads, name1, name2, ' '.join(args.quant_flags), extra_flags)
        outfile.write('%s\n' % cmd)
        voila_file_name = os.path.join(deltapsi_dir, '%s-%s.%s.voila' % (name1, name2, 'het' if method == 'heterogen' else method))
        cmd = 'voila tsv -f %s.tsv %s %s %s' % (voila_file_name, voila_file_name, build_dir, ' '.join(args.voila_flags))
        outfile.write('%s\n' % cmd)
    
    outdir = args.output
    outname = args.outname
    build_dir = os.path.join(outdir, 'build')
    ini_fname = os.path.join(outdir, 'settings_%s.ini' % outname)
    nthreads = args.nthreads

    majiq_ini = ConfigParser()
    majiq_ini.optionxform = str

    majiq_ini.add_section('info')
    majiq_ini.add_section('experiments')

    majiq_ini.set('info', 'genome', args.genome)
    majiq_ini.set('info', 'bamdirs', args.bam_dir)
    majiq_ini.set('info', 'readlen', str(args.readlen))

#    groups = defaultdict(set)
    groups = {}
    with open('%s/run_%s.sh' % (outdir, outname), 'w+') as outfile:
        outfile.write('#!/bin/bash\n')
        outfile.write('#\n')
        outfile.write('#This file has been generated using:\n#\t %s \n \n' % cmdline)
        outfile.write('. %s/activate\n' % args.bindir)
        outfile.write('\n')
        if not args.separate_builder:
            outfile.write('#Build samples\n')
            outfile.write('majiq build %s -o %s -c %s -j %d %s\n' % (args.annot_db, build_dir, ini_fname, nthreads, ' '.join(args.build_flags)))
            outfile.write('\n')
            outfile.write('\n')

        print(comps)
        for directory, samples in comps.items():
            outfile.write('#%s\n' % samples['cmd'])
            if args.separate_builder:
                build_dir = os.path.join(outdir, 'build_%s' % directory)
                ini_fname = '%s/settings_%s.ini' % (outdir, directory)
                outfile.write('majiq build %s -o %s -c %s -j %d %s\n' % (args.annot_db, build_dir, ini_fname, nthreads, ' '.join(args.build_flags)))
            name1, samps1 = samples['grp1']
            name2, samps2 = samples['grp2']
            grp1 = ' '.join([os.path.join(build_dir, '%s.majiq' % sample) for sample in samps1])
            grp2 = ' '.join([os.path.join(build_dir, '%s.majiq' % sample) for sample in samps2])

            for method in args.quant_methods:
                print_run(outfile, method, directory, grp1, grp2, name1, name2)
                if args.rr_pattern is not None and args.rr_pattern in directory:
                    print_run(outfile, method, directory + '_relaxed', grp1, grp2, name1, name2, extra_flags='--minreads 2 --minpos 2')

            # These lines assume the group names in the TSV files are of the format (tissue_name)a or (tissue_name)b.
            if args.separate_builder:
                majiq_ini.set('experiments', name1, ','.join(samps1))
                majiq_ini.set('experiments', name2, ','.join(samps2))
                with open(ini_fname, 'w') as ini_file:
                    majiq_ini.write(ini_file)
                majiq_ini.remove_option('experiments', name1)
                majiq_ini.remove_option('experiments', name2)
            else:
                for xx in samps1:
                    groups[xx] = xx
                for xx in samps2:
                    groups[xx] = xx
#                groups[name1].update(samps1)
#                groups[name2].update(samps2)

    if not args.separate_builder:
        for key, value in groups.items():
#            majiq_ini.set('experiments', key, ','.join(value))
            majiq_ini.set('experiments', key, value)
        with open(ini_fname, 'w') as ini_file:
            majiq_ini.write(ini_file)


def parse_settings(inFile):

    comparison = {}

    with open(inFile) as fp:
        for ln in fp.readlines():
            if ln.startswith('#') : continue
            tab = ln.strip().split()
            print (tab)
            comparison[tab[2]] = {"grp1": (tab[0], [xx for xx in tab[3].split(',')]),
                                  "grp2": (tab[1], [xx for xx in tab[4].split(',')]),
                                  "cmd" : ln.strip()}
    return comparison


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    common = argparse.ArgumentParser(add_help=False)
   
    common.add_argument('outname', action='store') 
    common.add_argument('-i', '--input_file', action='store')
    common.add_argument('-b', '--bam_dir', action='store')
    common.add_argument('-o', '--output', action='store')
    common.add_argument('-j', '--nthreads', action='store', type=int, default=4)
    common.add_argument('-db', '--annot_db', action='store')
    common.add_argument('--time', action='store_true', default=False)
    common.add_argument('--bindir', default='/usr/local/bin/')
    
    rmats = argparse.ArgumentParser(add_help=False)
    rmats.add_argument('--readtype', action='store', default='paired')
    rmats.add_argument('--libtype', action='store', default='fr-unstranded')
    rmats.add_argument('--readlen', action='store', type=int, default=100)
    rmats.add_argument('--threshold', action='store', type=float, default=0.2)


    whippet = argparse.ArgumentParser(add_help=False)
    whippet.add_argument('--threshold', action='store', type=float, default=0.2)
    whippet.add_argument('--index', action='store', required=True)
    whippet.add_argument('--annot-fa', action='store', default='fa.gz')
    whippet.add_argument('--fasta-ext', action='store', default='fasta')
    whippet.add_argument('--separate-run', action='store_true', default=False)

    suppa = argparse.ArgumentParser(add_help=False)
    suppa.add_argument('--transcript_fa', dest="transcriptfa", action='store', type=str)
    suppa.add_argument('--threshold', action='store', type=float, default=0.2)
    suppa.add_argument('--libtype', dest="libtype", action='store', type=str, default="A")
    suppa.add_argument('--salmon-dir', dest="salmondirbin", default='/usr/local/bin/')
    suppa.add_argument('--salmon-index', dest="salmon_index", required=True)
    suppa.add_argument('--fa-format', dest="fa_format", default='fastq')
    suppa.add_argument('--single_end', dest="paired", default=True, action='store_false')

    leafcutter = argparse.ArgumentParser(add_help=False)
    leafcutter.add_argument('--nreads', action='store', type=int, default=50)
    leafcutter.add_argument('--intronlen', action='store', type=int, default=500000)
    leafcutter.add_argument('--min_samples_x_intron', action='store', type=int, default=5)
    leafcutter.add_argument('--min_samples_x_group', action='store', type=int, default=3)
    leafcutter.add_argument('--min_coverage', action='store', type=int, default=20)

    majiq = argparse.ArgumentParser(add_help=False)
    majiq.add_argument('--genome', default='hg19')
    majiq.add_argument('--readlen', action='store', type=int, default=100)
    majiq.add_argument('--build-flags', nargs='*', default=[])
    majiq.add_argument('--quant-flags', nargs='*', default=[])
    majiq.add_argument('--voila-flags', nargs='*', default=['--show-all'])
    majiq.add_argument('--separate-builder', action='store_true')
    majiq.add_argument('--rr-pattern')
    majiq.add_argument('--quant-methods', nargs='+', choices=['deltapsi', 'heterogen'], default=['deltapsi'])

    subparsers = parser.add_subparsers(help='')
    
    parser_suppa = subparsers.add_parser('suppa', help='Generate Whippet call from input file', parents=[common, suppa])
    parser_suppa.set_defaults(func=gen_SUPPA_run, bindir='/home/jordi/software/suppa')

    parser_rmats = subparsers.add_parser('rmats', help='Generate rMATS call from input file', parents=[common, rmats])
    parser_rmats.set_defaults(func=gen_rMATS_run, bindir='/home/jordi/software/rMATS.4.0.2/rMATS-turbo-Linux-UCS4')

    parser_leafcutter = subparsers.add_parser('leafcutter', help='Generate leafcutter call from input file', parents=[common, leafcutter])
    parser_leafcutter.set_defaults(func=gen_leafcutter_run, bindir='/home/jordi/software/leafcutter')

    parser_majiq = subparsers.add_parser('majiq', help='Generate MAJIQ call from input file', parents=[common, majiq])
    parser_majiq.set_defaults(func=gen_majiq_run, bindir='/opt/venvs/py3/bin')

    parser_whippet = subparsers.add_parser('whippet', help='Generate Whippet call from input file', parents=[common, whippet])
    parser_whippet.set_defaults(func=gen_whippet_run, bindir='/home/jordi/.julia/v0.6/Whippet/bin')


    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)


    args = parser.parse_args()
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    print(args.bindir)
    comps = parse_settings(args.input_file)
    cmdline = ' ' .join(sys.argv)
    print(args.bindir)
    args.func(comps, cmdline, args)








