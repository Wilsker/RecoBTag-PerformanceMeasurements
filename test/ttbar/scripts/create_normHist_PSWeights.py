import ROOT, json
from ROOT import TFile, TTree, gDirectory,gROOT, TH1D, vector
from array import array
"""
Takes a directory on eos (starting from /store/...) and returns a list of all files with 'prepend' prepended
"""
def getEOSlslist(directory, mask='', prepend='root://eoscms//eos/cms'):
    from subprocess import Popen, PIPE
    print 'looking into: '+directory+'...'
    eos_cmd = '/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select'
    data = Popen([eos_cmd, 'ls', '/eos/cms/'+directory],stdout=PIPE)
    out,err = data.communicate()

    full_list = []

    ## if input file was single root file:
    if directory.endswith('.root'):
        if len(out.split('\n')[0]) > 0:
            return [prepend + directory]

    ## instead of only the file name append the string to open the file in ROOT
    for line in out.split('\n'):
        if len(line.split()) == 0: continue
        full_list.append(prepend + directory + '/' + line)

    ## strip the list of files if required
    if mask != '':
        stripped_list = [x for x in full_list if mask in x]
        return stripped_list

    ## return
    return full_list

def produceNormalisationHistogram(samplesList):

    norigEvents = ROOT.TH1D("","",1095,0,1095)
    for tag,sample in samplesList:
        xsec=sample[0]
        print 'tag: ', tag
        input_list=getEOSlslist(directory='/store/group/phys_btag/Commissioning/TTbar/KinMethod_StructuredDir_Moriond2018-02-20_ttbar/'+tag)
        #input_list=['/eos/cms/store/group/phys_btag/Commissioning/TTbar/KinMethod_StructuredDir_Moriond2018-02-20_ttbar/MC13TeV_TTJets_AH_training/JetTree_mc_FatJets_1.root']
        #input_list=['/eos/cms/store/group/phys_btag/Commissioning/TTbar/KinMethod_StructuredDir_Moriond2018-02-20_ttbar/MC13TeV_TTJets_AH_training/JetTree_mc_FatJets_1.root']

        print 'Number of files to process: ' , len(input_list)
        count_files = 0
        for f in input_list:
            count_files = count_files+1
            if count_files%200 == 0:
                print 'Processing file number: ', count_files
                print 'File: ', f
            fIn = TFile(f)
            # Declare ttree.
            myttree = gDirectory.Get("btagana/ttree")
            # Get entries in ttree.
            entries = myttree.GetEntriesFast()
            print '# entries: ', entries
            # Define histogram to hold sum of weights.
            temp_hist = ROOT.TH1D("","",1095,0,1095)
            # If first file, create new histogram.
            for entry in xrange(entries):
                # Load ttree to memory.
                ientry = myttree.LoadTree(entry)
                if ientry < 0:
                    break
                # Get values of branches for entry in ttree
                nb = myttree.GetEntry(entry)
                if nb <=0:
                    continue
                #for weight_index in xrange(1082, myttree.ttbar_nw):
                for weight_index in xrange(0, myttree.ttbar_nw):
                    temp_hist.Fill(weight_index,myttree.ttbar_w[weight_index])
            norigEvents.Add(temp_hist)
            fIn.Close()
    print "norigEvents.GetBinContent(0): " , norigEvents.GetBinContent(0)
    print "norigEvents.GetBinContent(1): " , norigEvents.GetBinContent(1)
    print "norigEvents.GetBinContent(2): " , norigEvents.GetBinContent(2)
    print "norigEvents.GetBinContent(3): " , norigEvents.GetBinContent(3)
    try:
        #Loop over the bins in norigEvents histogram and calculate ratio w.r.t. xsec in sample .json file.
        for xbin in xrange(1,norigEvents.GetNbinsX()+1):
            if norigEvents.GetBinContent(xbin)>0:
                print xsec/norigEvents.GetBinContent(xbin)
                norigEvents.SetBinContent(xbin,xsec/norigEvents.GetBinContent(xbin))
                norigEvents.SetBinError(xbin,0.)
    except:
        print 'No normalization histogram for ',tag
    return norigEvents

def main():
    outfile_name = 'norm_histos_PSWeights_allsamples.root'
    #Only run on ttbar else pointless.
    jsonFile = open('../data/samples_Run2017_MC_ttbarAH.json','r')
    samplesList=json.load(jsonFile,encoding='utf-8').items()

    fOut = TFile(outfile_name,'RECREATE')

    norm_hist = ROOT.TH1D("norma_hist","norma_hist",1095,0,1095)
    norm_hist = produceNormalisationHistogram(samplesList)
    norm_hist.SetDirectory(fOut)
    #norm_hist.Write()
    print 'fOut: ', fOut.ls()
    fOut.Write()
    fOut.Close()

main()
