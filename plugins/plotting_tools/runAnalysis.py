#!/usr/bin/python

import sys,os

rootcmd  = "root -b -x -q "
create   = "runData.C"
create2  = "runData2.C"
plot     = "analysisPlots.C"
plot2    = "analysisPlots2.C"



cmd = rootcmd+create
print cmd
os.system(cmd)

#vbadd  = "hadd VectorBoson_Background.root W*.root Z*.root"
vbadd  = "hadd SC5_MET_350/VectorBoson_Background.root SC5_MET_350/WJets.root SC5_MET_350/Z*Jets.root"
print vbadd
os.system(vbadd)

qcdadd = "hadd SC5_MET_350/QCD_Background.root SC5_MET_350/QCD_pt*.root"
print qcdadd
os.system(qcdadd)

smadd = "hadd SC5_MET_350/SM_Background.root SC5_MET_350/QCD_Background.root SC5_MET_350/VectorBoson_Background.root SC5_MET_350/TTJets.root"
print smadd
os.system(smadd)

#cmd = "mv *.root SC5_MET_350/"
#print cmd
#os.system(cmd)

#################################

#cmd = rootcmd+plot
#print cmd
#os.system(cmd)

#makepdfs = """gs -q -dSAFER -dNOPAUSE -dBATCH -sOutputFile=pre_cuts_plots.pdf -sDEVICE=pdfwrite -c .setpdfwrite -f ./*pre_cuts*.ps"""
#print makepdfs
#os.system(makepdfs)
#
#makepdfs = """gs -q -dSAFER -dNOPAUSE -dBATCH -sOutputFile=N1_cuts_plots.pdf -sDEVICE=pdfwrite -c .setpdfwrite -f ./h_N1*.ps"""
#print makepdfs
#os.system(makepdfs)
#
#makepdfs = """gs -q -dSAFER -dNOPAUSE -dBATCH -sOutputFile=individual_cuts_plots.pdf -sDEVICE=pdfwrite -c .setpdfwrite -f ./*individual_cuts*.ps"""
#print makepdfs
#os.system(makepdfs)
#
#makepdfs = """gs -q -dSAFER -dNOPAUSE -dBATCH -sOutputFile=post_cuts_plots.pdf -sDEVICE=pdfwrite -c .setpdfwrite -f ./*post_cuts*.ps h_selections.ps"""
#print makepdfs
#os.system(makepdfs)
#
#cmd = rootcmd+plot2
#print cmd
#os.system(cmd)
#
#makepdfs = """gs -q -dSAFER -dNOPAUSE -dBATCH -sOutputFile=pre_cuts_plots_backgrounds.pdf -sDEVICE=pdfwrite -c .setpdfwrite -f ./*pre_cuts*.ps"""
#print makepdfs
#os.system(makepdfs)
#
#makepdfs = """gs -q -dSAFER -dNOPAUSE -dBATCH -sOutputFile=N1_cuts_plots_backgrounds.pdf -sDEVICE=pdfwrite -c .setpdfwrite -f ./h_N1*.ps"""
#print makepdfs
#os.system(makepdfs)
#
#makepdfs = """gs -q -dSAFER -dNOPAUSE -dBATCH -sOutputFile=individual_cuts_plots_backgrounds.pdf -sDEVICE=pdfwrite -c .setpdfwrite -f ./*individual_cuts*.ps"""
#print makepdfs
#os.system(makepdfs)
#
#makepdfs = """gs -q -dSAFER -dNOPAUSE -dBATCH -sOutputFile=post_cuts_plots_backgrounds.pdf -sDEVICE=pdfwrite -c .setpdfwrite -f ./*post_cuts*.ps h_selections.ps"""
#print makepdfs
#os.system(makepdfs)
#
##copyfiles = """scp *.pdf sturdy@lxplus.cern.ch:~/public/html/.susyplots/sixth_iteration"""
##print copyfiles
##os.system(copyfiles)
#

##################################

cmd = rootcmd+create2
print cmd
os.system(cmd)

#vbadd  = "hadd VectorBoson_Background.root W*.root Z*.root"
vbadd  = "hadd SC5_MET_200/VectorBoson_Background.root SC5_MET_200/WJets.root SC5_MET_200/Z*Jets.root"
print vbadd
os.system(vbadd)

qcdadd = "hadd SC5_MET_200/QCD_Background.root SC5_MET_200/QCD_pt*.root"
print qcdadd
os.system(qcdadd)

smadd = "hadd SC5_MET_200/SM_Background.root SC5_MET_200/QCD_Background.root SC5_MET_200/VectorBoson_Background.root SC5_MET_200/TTJets.root"
print smadd
os.system(smadd)

#cmd = "mv *.root SC5_MET_200/"
#print cmd
#os.system(cmd)

#cmd = rootcmd+plot
#print cmd
#os.system(cmd)
#
#makepdfs = """gs -q -dSAFER -dNOPAUSE -dBATCH -sOutputFile=pre_cuts_plots.pdf -sDEVICE=pdfwrite -c .setpdfwrite -f ./*pre_cuts*.ps"""
#print makepdfs
#os.system(makepdfs)
#
#makepdfs = """gs -q -dSAFER -dNOPAUSE -dBATCH -sOutputFile=N1_cuts_plots.pdf -sDEVICE=pdfwrite -c .setpdfwrite -f ./h_N1*.ps"""
#print makepdfs
#os.system(makepdfs)
#
#makepdfs = """gs -q -dSAFER -dNOPAUSE -dBATCH -sOutputFile=individual_cuts_plots.pdf -sDEVICE=pdfwrite -c .setpdfwrite -f ./*individual_cuts*.ps"""
#print makepdfs
#os.system(makepdfs)
#
#makepdfs = """gs -q -dSAFER -dNOPAUSE -dBATCH -sOutputFile=post_cuts_plots.pdf -sDEVICE=pdfwrite -c .setpdfwrite -f ./*post_cuts*.ps h_selections.ps"""
#print makepdfs
#os.system(makepdfs)
#
##################################
#cmd = rootcmd+plot2
#print cmd
#os.system(cmd)
#
#makepdfs = """gs -q -dSAFER -dNOPAUSE -dBATCH -sOutputFile=pre_cuts_plots_backgrounds.pdf -sDEVICE=pdfwrite -c .setpdfwrite -f ./*pre_cuts*.ps"""
#print makepdfs
#os.system(makepdfs)
#
#makepdfs = """gs -q -dSAFER -dNOPAUSE -dBATCH -sOutputFile=N1_cuts_plots_backgrounds.pdf -sDEVICE=pdfwrite -c .setpdfwrite -f ./h_N1*.ps"""
#print makepdfs
#os.system(makepdfs)
#
#makepdfs = """gs -q -dSAFER -dNOPAUSE -dBATCH -sOutputFile=individual_cuts_plots_backgrounds.pdf -sDEVICE=pdfwrite -c .setpdfwrite -f ./*individual_cuts*.ps"""
#print makepdfs
#os.system(makepdfs)
#
#makepdfs = """gs -q -dSAFER -dNOPAUSE -dBATCH -sOutputFile=post_cuts_plots_backgrounds.pdf -sDEVICE=pdfwrite -c .setpdfwrite -f ./*post_cuts*.ps h_selections.ps"""
#print makepdfs
#os.system(makepdfs)
#
##copyfiles = """scp *.pdf sturdy@lxplus.cern.ch:~/public/html/.susyplots/sixth_iteration"""
##print copyfiles
##os.system(copyfiles)

#cmd = ".."
#print cmd
#os.chdir(cmd)
