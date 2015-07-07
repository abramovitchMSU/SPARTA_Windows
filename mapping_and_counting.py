__author__ = 'benkjohnson'


import os
import subprocess
import glob
from shutil import copy
import check_dependencies_windows
import qc_analysis

class Mapping_and_Counting(object):
    def __init__(self):

        return

    def bowtie(self, datalocation, analysislocation, options):
        """Run Bowtie for SE reads less than 50 bp in length.
        Will add the ability to run Bowtie2 for PE and SE with
        reads greater than 50 bp."""

        cd = check_dependencies_windows.CheckDependencies()
        qc = qc_analysis.QC_analysis()
        gff, genref = qc.findreferencefiles(datalocation)
        copy(genref, os.path.join(analysislocation, 'Bowtie'))
        genrefname = genref.split("/")[-1]
        copy(gff, os.path.join(analysislocation, 'HTSeq'))
        # subprocess.Popen("cp " + genref + " " + analysislocation + "/Bowtie", shell=True).wait()
        # subprocess.Popen("cp " + gff + " " + analysislocation + "/HTSeq", shell=True).wait()
        os.chdir(os.path.join(cd.getSPARTAdir(options), "Mapping_and_counting"))
        # os.chdir(cd.getSPARTAdir() + "/Mapping_and_counting")
        # if not os.path.lexists(os.path.join(cd.getSPARTAdir(), "Mapping_and_counting", "bowtie-1.1.1")):
        #     #This will be a problem for Windows users. Distribute with unzipped binaries?
        #     subprocess.call(["unzip", "bowtie-1.1.1-linux-x86_64.zip"], stdout=open(os.devnull, 'wb'))
        os.chdir(os.path.join(cd.getpwd(), "bowtie-1.1.1"))
        for file in os.listdir(os.path.join(analysislocation, "QC")):
            extension = file.split(".")[-1]
            if extension == "gz":
                subprocess.Popen("gzip -dc " + os.path.join(analysislocation, "QC", file) + " > " + os.path.join(analysislocation, "Bowtie", os.path.splitext(file)[0]), shell=True).wait()
            else:
                copy(os.path.join(analysislocation, "QC", file), os.path.join(analysislocation, "Bowtie"))

        if options.cleanup:
            for file in os.listdir(os.path.join(analysislocation, "QC")):
                extension = file.split(".")[-1]
                if extension == "gz" or extension in ["fq, fastq"]:
                    subprocess.Popen("rm " + os.path.join(analysislocation, "QC", "{file}".format(file=file)), shell=True).wait()

        print "Building the Bowtie index from the reference genome"
        #This could be a problem... Have the user specify if the index is large or small
        if options.verbose:
            subprocess.Popen(".\\bowtie-build-s.exe -f " + os.path.join(analysislocation, "Bowtie", genrefname) + " " + glob.glob(os.path.join(analysislocation, "Bowtie") + "/*.fa*")[0].split(".")[0], shell=True).wait()
        else:
            subprocess.Popen(".\\bowtie-build-s.exe -q -f " + os.path.join(analysislocation, "Bowtie", genrefname) + " " + glob.glob(os.path.join(analysislocation, "Bowtie") + "/*.fa*")[0].split(".")[0], shell=True).wait()

        allebwtfiles = glob.glob("*.ebwt")[:]
        for ebwtfile in allebwtfiles:
            copy(ebwtfile, os.path.join(analysislocation, "Bowtie"))
            # subprocess.Popen("cp " + ebwtfile + " " + analysislocation + "/Bowtie/", shell=True).wait()
        print "Mapping reads to the reference genome with Bowtie"
        for file in os.listdir(os.path.join(analysislocation, "Bowtie")):
            extension = os.path.splitext(file)[1]
            if extension == ".fq" or extension == ".fastq":
                fname = os.path.splitext(file)[0]
                strippedfile = fname[len('trimmed'):]
                #Again this could be a problem, have the user define whether the index is large or small
                if options.verbose:
                    subprocess.Popen(".\\bowtie-align-s.exe -S --threads {threads} ".format(threads=options.threads) + glob.glob(os.path.join(analysislocation, "Bowtie") + "/*.fa*")[0].split(".")[0] + " " + os.path.join(analysislocation, "Bowtie", file) + " > " + os.path.join(analysislocation, "Bowtie", "align" + strippedfile + ".sam"), shell=True).wait()
                elif options.mismatch:
                    subprocess.Popen(".\\bowtie-align-s.exe -S --threads {threads} -v {mismatch} ".format(threads=options.threads, mismatch=options.mismatch) + glob.glob(os.path.join(analysislocation, "Bowtie") + "/*.fa*")[0].split(".")[0] + " " + os.path.join(analysislocation, "Bowtie", file) + " > " + os.path.join(analysislocation, "Bowtie", "align" + strippedfile + ".sam"), shell=True).wait()
                elif options.otherbowtieoptions:
                    subprocess.Popen(".\\bowtie-align-s.exe -S --threads {threads} {otherbowtieoptions}".format(threads=options.threads, otherbowtieoptions=options.otherbowtieoptions) + glob.glob(os.path.join(analysislocation, "Bowtie") + "/*.fa*")[0].split(".")[0] + " " + os.path.join(analysislocation, "Bowtie", file) + " > " + os.path.join(analysislocation, "Bowtie", "align" + strippedfile + ".sam"), shell=True).wait()
                elif options.mismatch and options.otherbowtieoptions:
                    subprocess.Popen(".\\bowtie-align-s.exe -S --threads {threads} -v {mismatch} {otherbowtieoptions} ".format(threads=options.threads, mismatch=options.mismatch, otherbowtieoptions=options.otherbowtieoptions) + glob.glob(os.path.join(analysislocation, "Bowtie") + "/*.fa*")[0].split(".")[0] + " " + os.path.join(analysislocation, "Bowtie", file) + " > " + os.path.join(analysislocation, "Bowtie", "align" + strippedfile + ".sam"), shell=True).wait()
                else:
                    subprocess.Popen(".\\bowtie-align-s.exe -S --quiet --threads {threads} ".format(threads=options.threads) + glob.glob(os.path.join(analysislocation, "Bowtie") + "/*.fa*")[0].split(".")[0] + " " + os.path.join(analysislocation, "Bowtie", file) + " > " + os.path.join(analysislocation, "Bowtie", "align" + strippedfile + ".sam"), shell=True).wait()

        return

    def htseq(self, analysislocation, options):
        """Run htseq-count to count gene features post-Bowtie mapping"""

        cd = check_dependencies_windows.CheckDependencies()
        os.chdir(os.path.join(cd.getSPARTAdir(options), "Mapping_and_counting"))
        gff = glob.glob(os.path.join(analysislocation, "HTSeq") + "/*.g*")[0]
        print "Counting gene features with HTSeq"
        for mapfile in os.listdir(os.path.join(analysislocation, "Bowtie")):
            extension = os.path.splitext(mapfile)[1]
            if extension == ".sam":
                fname = os.path.splitext(mapfile)[0]
                strippedmapfile = fname[len('align'):]
                if options.verbose:
                    subprocess.Popen("python -m HTSeq.scripts.count --mode={mode} --stranded={stranded} --order={order} --type={type} -a {minqual} --idattr={idattr} ".format(mode=options.mode, stranded=options.stranded, order=options.order, type=options.type, minqual=options.minqual, idattr=options.idattr) + os.path.join(analysislocation, "Bowtie", mapfile) + " " + gff + " > " + os.path.join(analysislocation, "HTSeq", "map" + strippedmapfile + ".sam"), shell=True).wait()
                else:
                    subprocess.Popen("python -m HTSeq.scripts.count --quiet --mode={mode} --stranded={stranded} --order={order} --type={type} -a {minqual} --idattr={idattr} ".format(mode=options.mode, stranded=options.stranded, order=options.order, type=options.type, minqual=options.minqual, idattr=options.idattr) + os.path.join(analysislocation, "Bowtie", mapfile) + " " + gff + " > " + os.path.join(analysislocation, "HTSeq", "map" + strippedmapfile + ".sam"), shell=True).wait()

        if options.cleanup:
            for file in os.listdir(os.path.join(analysislocation, "Bowtie")):
                subprocess.Popen("rm " + os.path.join(analysislocation, "Bowtie", "{file}".format(file=file)), shell=True).wait()

        return