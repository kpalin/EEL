# -*- coding: UTF-8 -*- 

##
## $Log$
## Revision 1.4.2.2  2005/05/10 13:13:00  kpalin
## Fixed the default directory thing. Now it should be in working order.
##
## Revision 1.4.2.1  2005/05/09 07:16:46  kpalin
## Remember browsed directory and fixed button labeling.
##
## Revision 1.4  2005/03/22 13:17:13  kpalin
## Merged some fixes surfacing from testing the public version.
##
## Revision 1.3.2.1  2005/03/22 12:19:26  kpalin
## Fixed the problems that came up with testing.
##
## Revision 1.3  2005/01/13 13:16:42  kpalin
## Moved the requesting of sequences to be aligned to Python side
## of stuff. Much better.
##
## Revision 1.2  2005/01/12 13:55:09  kpalin
## Path handling fixed.
##
## Revision 1.1  2005/01/12 13:34:55  kpalin
## Added Tkinter/Tix Graphical user interface and command -no-gui to
## avoid it.
##
##



from eellib.Interface import Interface
from eellib.tixgui import tixgui,fontStr


import tkMessageBox
from Tix import *
#from Tkinter import *

import os.path,os


class alignUI(Frame):
    "Frame for starting the alignment"
    def __init__(self, *args, **params):
    ## Standard heading: initialization
        apply(Frame.__init__, (self,) + args, params)
        self._parent = None
        if len(args) != 0: self._parent = args[0]

        self._widgets={}

        self._widgets["params"]=alignParams(self)
        self._widgets["params"].grid(column=1,row=1,sticky="w",padx=20)

        self._widgets["data"]=alignData(self)
        self._widgets["data"].grid(column=1,row=2)

        self._widgets["data"].setDefaultRadio(self._parent.master.haveMatches())
        
        self._widgets["commands"]=alignCommand(self)
        self._widgets["commands"].grid(column=1,row=3)

    def cancel_button(self):
        self.master.destroy()
        pass


    def align(self,seq1=None,seq2=None):
        if int(self.loadFile.get())>0:
            fname=self.siteFileName.get()
            try:
                if not os.access(fname,os.R_OK):
                    return 0
            except(ImportError,AttributeError):
                pass
        else:
            fname='.'

        lambdaVal=float(self.lambdaStr.get())
        xi=float(self.xiStr.get())
        mu=float(self.muStr.get())
        nu=float(self.nuStr.get())
        nuc_per_rot=float(self.nprStr.get())
        resultReq=int(self.resultCount.get())

        self._parent.master.align(fname,resultReq,lambdaVal,xi,mu,nu,nuc_per_rot,seq1,seq2)
        if hasattr(self._parent.master,"alignment"):
            self._parent.master._widgets["saveAlign_button"].config(state="normal")
            self.master.destroy()
            pass
        
    def align_button(self):
        #print self.loadFile.get()
        try:
            self.align()
        except AttributeError,val:
            if len(val.args)==1 or (len(val.args)>1 and len(val.args[1])<2):
                tkMessageBox.showerror("Error..",val.args[0])
            else:
                root=Toplevel(self)
                root.title("Select sequences..")
                self._widgets["selectSeqWindow"]=selectSequences(root,seqsToChoose=val.args[1])
                self._widgets["selectSeqWindow"].pack()
                #tkMessageBox.showerror("Error..",val.args[0])
                msg,seqs=val
          


class selectSequences(Frame):
    def __init__(self, *args, **params):
    ## Standard heading: initialization
        self.seqsToChoose=params["seqsToChoose"]
        del(params["seqsToChoose"])
        
        apply(Frame.__init__, (self,) + args, params)
        self._parent = None
        if len(args) != 0: self._parent = args[0]

        self._widgets={}
        self.firstSeqName=StringVar()
        self.secondSeqName=StringVar()
        
        self._widgets["firstSeqCombo"]=ComboBox(self,label="Select first sequence:",variable=self.firstSeqName,listcmd=lambda x="firstSeqCombo":self.addComboList(x),dropdown=TRUE,command=lambda x:self.setValidSeqs("firstSeqCombo",x))
        self._widgets["firstSeqCombo"].grid(row=1,column=1,columnspan=2,sticky="w")

        self._widgets["secondSeqCombo"]=ComboBox(self,label="Select second sequence:",variable=self.secondSeqName,listcmd=lambda x="secondSeqCombo":self.addComboList(x),dropdown=TRUE,command=lambda x:self.setValidSeqs("secondSeqCombo",x))
        self._widgets["secondSeqCombo"].grid(row=2,column=1,columnspan=2,sticky="w")


        self._widgets["ok_button"]=Button(self,name="button_OK",font=fontStr,text="OK",command=self.button_OK)
        self._widgets["ok_button"].grid(column=2,row=3)

        self._widgets["cancel_button"]=Button(self,name="button_cancel",font=fontStr,text="Cancel",command=self.button_cancel)
        self._widgets["cancel_button"].grid(column=1,row=3)

        self.firstSeqName.set(self.seqsToChoose[0])
        self.secondSeqName.set(self.seqsToChoose[1])

    def setValidSeqs(self,fixed,newvalue):
        for x in self.seqsToChoose:
            if x!=newvalue:
                break
        if fixed=="firstSeqCombo" and newvalue==self.secondSeqName.get():
            self.secondSeqName.set(x)
        if fixed=="secondSeqCombo" and newvalue==self.firstSeqName.get():
            self.firstSeqName.set(x)
            

    def button_cancel(self):
        self._parent.destroy()

    def button_OK(self):
        
        self._parent.master.align(self.firstSeqName.get(),self.secondSeqName.get())

    
    def addComboList(self,comboStr):
        for i in self.seqsToChoose:
            self._widgets[comboStr].subwidget("listbox").insert(END,i)

class alignData(Frame):
    def __init__(self, *args, **params):
    ## Standard heading: initialization
        apply(Frame.__init__, (self,) + args, params)
        self._parent = None
        if len(args) != 0: self._parent = args[0]

        self._widgets={}
        self._parent.loadFile=StringVar()
        self._widgets["computed_radio"]=Radiobutton(self,name="computed_radio",font=fontStr,text="Computed binding sites",value="0",variable=self._parent.loadFile)

        self._widgets["computed_radio"].grid(column=1,row=1,sticky="w")

        self._widgets["load_radio"]=Radiobutton(self,name="load_radio",font=fontStr,text="Saved binding sites:",value="1",variable=self._parent.loadFile)

        self._widgets["load_radio"].grid(column=1,row=2,sticky="w")

        self._parent.siteFileName=StringVar()
        self._widgets["load_entry"]=FileEntry(self,name="file_entry",activatecmd=self.customEntryDial,variable=self._parent.siteFileName)
        self._widgets["load_entry"].grid(column=2,row=2)


    def setDefaultRadio(self,hascomputed):
        if hascomputed:
            self._widgets["computed_radio"].config(state="normal")
            self._widgets["computed_radio"].select()
        else:
            self._widgets["computed_radio"].config(state="disabled")
            self._widgets["load_radio"].select()

    def customEntryDial(self):
        _self=self._widgets["load_entry"]
        fdialStr=_self.tk.call(_self._w,"filedialog")
        fsboxStr=_self.tk.call(fdialStr,"subwidget","fsbox")
        
        _self.tk.call(fsboxStr,"configure","-pattern","*.gff")
        _self.tk.call(fsboxStr,"configure","-directory",self._parent._parent.master.defaultDir)
        pass

class alignCommand(Frame):
    def __init__(self, *args, **params):
    ## Standard heading: initialization
        apply(Frame.__init__, (self,) + args, params)
        self._parent = None
        if len(args) != 0: self._parent = args[0]

        self._widgets={}

        self._widgets["cancel"]=Button(self,name="cancel_button",font=fontStr,text="Cancel",command=self._parent.cancel_button)
        self._widgets["cancel"].grid(column=1,row=1)
        
        self._widgets["align"]=Button(self,name="align_button",font=fontStr,text="Align",command=self._parent.align_button)
        self._widgets["align"].grid(column=2,row=1)
        
class alignParams(Frame):
    def __init__(self, *args, **params):

    ## Standard heading: initialization
        apply(Frame.__init__, (self,) + args, params)
        self._parent = None
        if len(args) != 0: self._parent = args[0]

        self._widgets={}

        self._parent.lambdaStr=StringVar()
        self._parent.lambdaStr.set("2.0")
        self._widgets["lambda"]=LabelEntry(self,label="Lambda:")
        self._widgets["lambda"].entry.config(textvariable=self._parent.lambdaStr)
        self._widgets["lambda"].grid(column=1,row=1,sticky="w")



        self._parent.xiStr=StringVar()
        self._parent.xiStr.set("200.0")
        self._widgets["xi"]=LabelEntry(self,label="Xi:")
        self._widgets["xi"].entry.config(textvariable=self._parent.xiStr)
        self._widgets["xi"].grid(column=1,row=2,sticky="w")

        self._parent.muStr=StringVar()
        self._parent.muStr.set("0.5")

        self._widgets["mu"]=LabelEntry(self,label="Mu:")
        self._widgets["mu"].entry.config(textvariable=self._parent.muStr)
        self._widgets["mu"].grid(column=1,row=3,sticky="w")

        self._parent.nuStr=StringVar()
        self._parent.nuStr.set("200.0")
        self._widgets["nu"]=LabelEntry(self,label="Nu:")
        self._widgets["nu"].entry.config(textvariable=self._parent.nuStr)
        self._widgets["nu"].grid(column=1,row=4,sticky="w")


        self._parent.nprStr=StringVar()
        self._parent.nprStr.set("10.4")
        self._widgets["nuc_per_rot"]=LabelEntry(self,label="Nucleotides Per Rotation:")
        self._widgets["nuc_per_rot"].entry.config(textvariable=self._parent.nprStr)
        self._widgets["nuc_per_rot"].grid(column=1,row=5,sticky="w")

        self._parent.resultCount=StringVar()
        self._parent.resultCount.set("100")
        self._widgets["res_count"]=LabelEntry(self,label="Number of sub-optimal alignments:")
        self._widgets["res_count"].entry.config(textvariable=self._parent.resultCount)
        self._widgets["res_count"].grid(column=1,row=6,sticky="w")


class getTFBS(Frame):
    "Frame for getting Transcription Factor Binding Sites"
    def __init__(self, *args, **params):
    ## Standard heading: initialization
        apply(Frame.__init__, (self,) + args, params)
        self._parent = None
        if len(args) != 0: self._parent = args[0]

        self._widgets={}

        self._widgets["background"]=_bgSelect(self)
        self._widgets["background"].grid(column=1,row=1)

        self._widgets["other"]=_tfbsOther(self)
        self._widgets["other"].grid(column=1,row=2)


    def getTFBS(self,cutoff=9.0):
        bgStr=self._widgets["background"].bg_variable.get()
        if bgStr=="radio_single":
            if not self._widgets["background"]._widgets["single"].hasNormalFreq():
                return
            self._widgets["background"].selectSingle()
            
        #self._widgets["background"]._widgets[bgStr].select()
        isAbsolute=int(self._widgets["other"].bg_mode.get())

        self._parent.master.getTFBS(cutoff,isAbsolute)
        self._parent.master._widgets["saveSites_button"]["state"]="normal"
        self._parent.destroy()
        


class _tfbsOther(Frame):
    "Frame for giving the cutoff and start the tfbs scan"
    def __init__(self, *args, **params):
    ## Standard heading: initialization
        apply(Frame.__init__, (self,) + args, params)
        self._parent = None
        if len(args) != 0: self._parent = args[0]


        self._widgets={}
        self._widgets["label"]=Label(self,name="label",font=fontStr,text="Cutoff: ")
        self._widgets["label"].grid(column=1,row=1,sticky="w")

        self.cutoff=StringVar()

        self._widgets['cutoff']=Entry(self,name='cutoff', font=fontStr, textvariable=self.cutoff)
        self.cutoff.set("9")
        self._widgets["cutoff"].grid(column=2,row=1,sticky="w")

        self.bg_mode=StringVar()

        self._widgets["absolute_radio"]=Radiobutton(self,name="radio_abs",font=fontStr,text="Absolute",value="1",variable=self.bg_mode)
        self._widgets["absolute_radio"].grid(column=3,row=1)
        self._widgets["absolute_radio"].select()
        self._widgets["relative_radio"]=Radiobutton(self,name="radio_rel",font=fontStr,text="Relative",value="0",variable=self.bg_mode)
        self._widgets["relative_radio"].grid(column=4,row=1)


        self._widgets["getTFBS"]=Button(self,name="gettfbs_button",font=fontStr,text="Get TFBS",command=self.getTFBS_button)
        self._widgets["getTFBS"].grid(column=5,row=1)


    def getTFBS_button(self):
        self._parent.getTFBS(float(self.cutoff.get()))


class _bgSelect(Frame):
    "Frame for selecting parameters for background distribution"
    def __init__(self, *args, **params):
    ## Standard heading: initialization
        apply(Frame.__init__, (self,) + args, params)
        self._parent = None
        if len(args) != 0: self._parent = args[0]

        self.bg_variable=StringVar()

        self._widgets={}
        self._widgets["single"]=_bgFreq(self)
        self._widgets["single"].grid(column=2,row=1)

        self.markovFile=StringVar()
        self._widgets["markov"]=FileEntry(self,name="markov",activatecmd=self.newFileEntry,command=self.setMarkov,variable=self.markovFile)

        self._widgets["markov"].grid(column=2,row=2,sticky="w")

        self._widgets["radio_single"]=Radiobutton(self,name="radio_single",text="Single background: ",variable=self.bg_variable,value="radio_single",command=self.selectSingle)
        self._widgets["radio_single"].grid(column=1,row=1)
        
        self._widgets["radio_markov"]=Radiobutton(self,name="radio_markov",text="Markov background: ",variable=self.bg_variable,value="radio_markov",command=self.selectMarkov)
        self._widgets["radio_markov"].grid(column=1,row=2)
        self._widgets["radio_single"].select()



        self.grid_rowconfigure(1, weight=0, minsize=20)
        self.grid_rowconfigure(2, weight=1, minsize=20)
        self.grid_columnconfigure(1, weight=0, minsize=10)
        self.grid_columnconfigure(2, weight=0, minsize=100)



    def newFileEntry(self):
        _self=self._widgets["markov"]
        fdialStr=_self.tk.call(_self._w,"filedialog")
        fsboxStr=_self.tk.call(fdialStr,"subwidget","fsbox")
        
        _self.tk.call(fsboxStr,"configure","-pattern","*.bg")
        _self.tk.call(fsboxStr,"configure","-directory",self._parent._parent.master.defaultDir)
        #print repr(self._widgets["markov"].file_dialog())
        pass

    def selectSingle(self):
        self._parent._parent.master.setBGFreq(self._widgets["single"].getNormalFreq())
        
    def selectMarkov(self):
        self._widgets["markov"].button.invoke()
        fname=self.markovFile.get()
        
        self.setMarkov(fname)

    def setMarkov(self,fname):
        try:
            self._parent._parent.master.setMarkovBG(fname)
            self._widgets["radio_markov"].select()
        except (SyntaxError,IOError,IndexError):
            self._widgets['radio_single'].select()
        


class _bgFreq(Frame):
    "Frame for single nucleotide background frequencies"
    def __init__(self, *args, **params):
    ## Standard heading: initialization
        apply(Frame.__init__, (self,) + args, params)
        self._parent = None
        if len(args) != 0: self._parent = args[0]



        self._widgets={}

        self.freqA=StringVar()
        self._widgets["labelA"]=Label(self,name="labelA",font=fontStr,text="'A' Frequency")
        self._widgets["labelA"].grid(row=1,column=1)
        self._widgets['freqA']=Entry(self,name='freqA', font=fontStr, textvariable=self.freqA)
        self._widgets["freqA"].grid(row=2,column=1)


        self.freqC=StringVar()
        self._widgets["labelC"]=Label(self,name="labelC",font=fontStr,text="'C' Frequency")
        self._widgets["labelC"].grid(row=1,column=2)
        self._widgets['freqC']=Entry(self,name='freqC', font=fontStr, textvariable=self.freqC)
        self._widgets["freqC"].grid(row=2,column=2)

        self.freqG=StringVar()
        self._widgets["labelG"]=Label(self,name="labelG",font=fontStr,text="'G' Frequency")
        self._widgets["labelG"].grid(row=1,column=3)
        self._widgets['freqG_entry']=Entry(self,name='freqG_entry',textvariable=self.freqG)
        self._widgets["freqG_entry"].grid(row=2,column=3)

        self.freqT=StringVar()
        self._widgets["labelT"]=Label(self,name="labelT",font=fontStr,text="'T' Frequency")
        self._widgets["labelT"].grid(row=1,column=4)
        self._widgets['freqT']=Entry(self,name='freqT', font=fontStr, textvariable=self.freqT)
        self._widgets["freqT"].grid(row=2,column=4)


    ## Resize behavior(s)
        self.grid_rowconfigure(1, weight=0, minsize=10)
        self.grid_rowconfigure(2, weight=1, minsize=10)
        self.grid_columnconfigure(1, weight=1, minsize=10)
        self.grid_columnconfigure(2, weight=1, minsize=10)
        self.grid_columnconfigure(3, weight=1, minsize=10)
        self.grid_columnconfigure(4, weight=1, minsize=10)

        self.freqA.set("0.25")
        self.freqC.set("0.25")
        self.freqG.set("0.25")
        self.freqT.set("0.25")


    def hasNormalFreq(self):
        a,c,g,t=self.freqA.get(),self.freqC.get(),self.freqG.get(),self.freqT.get()
        self.getNormalFreq()
        A,C,G,T=self.freqA.get(),self.freqC.get(),self.freqG.get(),self.freqT.get()
        if A==a and C==c and G==g and T==t:
            return 1
        else:
            return 0

    def getNormalFreq(self):
        
        a,c,g,t=float(self.freqA.get()),float(self.freqC.get()),float(self.freqG.get()),float(self.freqT.get())
        tot=1.0*sum([a+c+g+t])
        a,c,g,t=a/tot,c/tot,g/tot,t/tot
        self.freqA.set(str(a))
        self.freqC.set(str(c))
        self.freqG.set(str(g))
        self.freqT.set(str(t))
        return (a,c,g,t)

        
class eelGUI(Interface,tixgui):
    def __init__(self,*args,**keywd):
        Interface.__init__(self)
        tixgui.__init__(self,*args,**keywd)
        # Give consistent default directory. Changes when user changes dir.
        self.defaultDir=os.getcwd()

        self.a,self.c,self.g,self.t=0.25,0.25,0.25,0.25

  ## METHOD saveAlign_button:
  ## ~~~~~~
    def saveAlign_button(self):
        root=Toplevel(self)
        self._widgets["showAlignWindow"]=showAligns(root)
        self._widgets["showAlignWindow"].pack(fill=BOTH)
        root.title("Conserved Enhancer Elements..")
        
            
        pass

  ## METHOD align_button:
  ## ~~~~~~
    def align_button(self):
        root=Toplevel(self)
        self._widgets["alignWindow"]=alignUI(root)
        self._widgets["alignWindow"].pack(side=TOP,fill=BOTH,expand=1)
        root.title("Align Binding Sites..")
    
        pass


  ## METHOD saveSites_button:
  ## ~~~~~~
    def saveSites_button(self):
        root=Toplevel(self)
        self._widgets["saveWindow"]=showSites(root)
        self._widgets["saveWindow"].pack(side=TOP,fill=BOTH,expand=1)
        root.title("Putative TF Binding sites")
        pass

  ## METHOD getSites_button:
  ## ~~~~~~
    def getSites_button(self):
        root=Toplevel(self)
        self._widgets["bgFreq"]=getTFBS(root)
        self._widgets["bgFreq"].pack(side=TOP, fill=BOTH, expand=1)
        root.title("Get TFBS..")
        
        pass



    def removeSeq(self):
        listbox=self._widgets["sequence_listbox"].subwidget("listbox")
        for i in listbox.curselection():
            seqName=listbox.get(i)
            self.removeSequence(seqName)
            listbox.delete(i)
        self.resetButtonStatus()

    def removeMat(self):
        listbox=self._widgets["matrices_listbox"].subwidget("listbox")
        for i in listbox.curselection():
            self.removeMatrix(int(i))
            #listbox.delete(i)
        self.updateMatrixList()
        #print mat
        

    def updateSequenceList(self):
        listbox=self._widgets["sequence_listbox"].subwidget("listbox")
        listbox.delete(0,END)
        for i in self.seq.getNames():
            listbox.insert(END,i)

        self.resetButtonStatus()
    
    def resetButtonStatus(self):
        matCount=len(self._widgets["matrices_listbox"].subwidget("listbox").get(0,END))
        seqCount=len(self._widgets["sequence_listbox"].subwidget("listbox").get(0,END))
        
        if matCount>0 and seqCount>0:
            self._widgets["getSites_button"]["state"]="normal"
        else:
            self._widgets["getSites_button"]["state"]="disabled"
            
    def addAllSequences(self):
        filelist=self._widgets["sequences_pick"].subwidget("filelist").subwidget("listbox")

        self.defaultDir=self._widgets["sequences_pick"]["directory"]
                
        files=filelist.get(0,END)
        self.addSequence([os.path.join(self.defaultDir,x) for x in files])
        self.updateSequenceList()

    def addSequences(self,filename):
        self.defaultDir=self._widgets["sequences_pick"]["directory"]
        self.addSequence([filename])
        self.updateSequenceList()


    def updateMatrixList(self):

        self._widgets["matrices_listbox"].subwidget("listbox").delete(0,END)
        for i in self.matlist:
            #i.setBGfreq(self.a,self.c,self.g,self.t)
            self._widgets["matrices_listbox"].subwidget("listbox").insert(END,"%s (I=%g)"%(i.name,i.InfoContent))
        self.resetButtonStatus()


    def addAllMatrices(self):
        filelist=self._widgets["mat_pick"].subwidget("filelist").subwidget("listbox")


        self.defaultDir=self._widgets["mat_pick"]["directory"]
        
        files=filelist.get(0,END)
        self.addMatrix([os.path.join(self.defaultDir,x) for x in files])
        self.updateMatrixList()
        
    def addMatrices(self,filename):
        
        self.defaultDir=self._widgets["mat_pick"]["directory"]
        self.addMatrix([filename])
        self.updateMatrixList()


class showAligns(Frame):
    "Frame for displaying and saving putative binding sites"
    def __init__(self, *args, **params):
    ## Standard heading: initialization
        apply(Frame.__init__, (self,) + args, params)
        self._parent = None
        if len(args) != 0: self._parent = args[0]

        self._widgets={}

        self._widgets["scrollText"]=ScrolledText(self,scrollbar="auto")
        #self._widgets["scrollText"].config(textfont="Courier")
        self._widgets["scrollText"].text["font"]="courier"

        from eellib import Output
        siteStr=Output.formatalign(self._parent.master.alignment,self._parent.master.seq)

        
        
        self._widgets["scrollText"].text.insert(END,siteStr)
        self._widgets["scrollText"].text.config(state="disabled")
        
        self._widgets["scrollText"].grid(column=1,row=1,sticky="nesw")
        self._widgets["scrollText"].text.config(wrap="none")
        
        
        self._widgets["buttonFrame"]=Frame(self)

        self._widgets["button_OK"]=Button(self._widgets["buttonFrame"],name="button_OK",font=fontStr,text="Done",command=self.button_OK)
        self._widgets["button_OK"].grid(column=3,row=1)

        self._widgets["button_Save"]=Button(self._widgets["buttonFrame"],name="button_Save",font=fontStr,text="Save human readable..",command=self.button_Save)
        self._widgets["button_Save"].grid(column=2,row=1)


        self._widgets["buttonFrame"].grid(column=1,row=2)

        self._widgets["button_SaveGFF"]=Button(self._widgets["buttonFrame"],name="button_SaveGFF",font=fontStr,text="Save computer readable..",command=self.button_SaveGFF)
        self._widgets["button_SaveGFF"].grid(column=1,row=1)


        self._widgets["buttonFrame"].grid(column=1,row=2,sticky="S")

        self.grid_rowconfigure(1,weight=1)
        self.grid_columnconfigure(1,weight=1)




    def button_OK(self):
        self._parent.destroy()

    def button_Save(self):
        #self.save_w=Toplevel(self)
        #self.save_w.title("Save sites..")
        self._widgets["saveBox"]=FileSelectDialog(self,command=self.savealign_click)
        self._widgets["saveBox"].fsbox.config(pattern="*.aln",directory=self._parent.master.defaultDir)

        self._widgets["saveBox"].popup()


    def savealign_click(self,*others):
        self._parent.master.defaultDir=self._widgets["saveBox"].fsbox["directory"]
        self._parent.master.savealign(*others)

    def savealignGFF_click(self,*others):
        self._parent.master.defaultDir=self._widgets["saveBox"].fsbox["directory"]
        self._parent.master.savealignGFF(*others)

    def button_SaveGFF(self):
        #self.save_w=Toplevel(self)
        #self.save_w.title("Save sites..")
        self._widgets["saveBox"]=FileSelectDialog(self,command=self.savealignGFF_click)
        self._widgets["saveBox"].fsbox.config(pattern="*.gff",directory=self._parent.master.defaultDir)

        self._widgets["saveBox"].popup()
        



class showSites(Frame):
    "Frame for displaying and saving putative binding sites"
    def __init__(self, *args, **params):
    ## Standard heading: initialization
        apply(Frame.__init__, (self,) + args, params)
        self._parent = None
        if len(args) != 0: self._parent = args[0]

        self._widgets={}

        self._widgets["scrollText"]=ScrolledText(self,scrollbar="auto")
        siteStr=self._parent.master.getmatchStr()
        
        self._widgets["scrollText"].text.insert(END,siteStr)
        self._widgets["scrollText"].text.config(state="disabled")
        
        self._widgets["scrollText"].grid(column=1,row=1,sticky="nswe")
        self._widgets["scrollText"].text.config(wrap="none")

        self._widgets["buttonFrame"]=Frame(self)

        self._widgets["button_OK"]=Button(self._widgets["buttonFrame"],name="button_OK",font=fontStr,text="Done",command=self.button_OK)
        self._widgets["button_OK"].grid(column=2,row=1)

        self._widgets["button_Save"]=Button(self._widgets["buttonFrame"],name="button_Save",font=fontStr,text="Save..",command=self.button_Save)
        self._widgets["button_Save"].grid(column=1,row=1)


        self._widgets["buttonFrame"].grid(column=1,row=2)
        self.grid_rowconfigure(1,weight=1)
        self.grid_columnconfigure(1,weight=1)

    def button_OK(self):
        self._parent.destroy()

    def button_Save(self):
        #self.save_w=Toplevel(self)
        #self.save_w.title("Save sites..")
        self._widgets["saveBox"]=FileSelectDialog(self,command=self.saveSelect)
        self._widgets["saveBox"].fsbox.config(pattern="*.gff")
        self._widgets["saveBox"].fsbox.config(directory=self._parent.master.defaultDir)

        self._widgets["saveBox"].popup()
        

        

    def saveSelect(self,fname):
        self._parent.master.defaultDir=self._widgets["saveBox"].fsbox["directory"]
        self._parent.master.savematch(fname)
        



if __name__ == '__main__':
  root = Tk()
  o = eelGUI(root)
  o.pack(side=TOP, fill=BOTH, expand=1)
  root.mainloop()
