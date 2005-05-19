#!/usr/bin/env python
# -*- coding: UTF-8 -*- 

##
## $Log$
## Revision 1.1.4.2  2005/05/10 13:13:00  kpalin
## Fixed the default directory thing. Now it should be in working order.
##
## Revision 1.1.4.1  2005/05/09 07:17:08  kpalin
## Fixed labeling etc.
##
## Revision 1.1  2005/01/12 13:34:55  kpalin
## Added Tkinter/Tix Graphical user interface and command -no-gui to
## avoid it.
##
##


import eelgui

from Tkinter import *

from Tix import *

fontStr='-*-helvetica-Medium-R-Normal-*-*-140-*-*-*-*-*-*'

class tixgui(Frame):
  def __init__(self, *args, **params):
    ## Standard heading: initialization
    apply(Frame.__init__, (self,) + args, params)
    self._parent = None
    if len(args) != 0: self._parent = args[0]

    self._parent.title("Enhancer Element Locator")
    
    self._init_before()
    ## Widget creation
    self._widgets = {}
    self._widgets['sequence_listbox'] = ScrolledListBox(self, name='sequence_listbox',scrollbar="auto", height=0, width=0,command=self.removeSeq)
    self._widgets['sequence_listbox'].grid(column=1, row=1, sticky='nesw')
    self._widgets['matrices_listbox'] = ScrolledListBox(self, name='matrices_listbox', scrollbar="auto",height=0, width=0,command=self.removeMat)
    self._widgets['matrices_listbox'].grid(column=2, row=1, sticky='nesw')
    self._widgets['sequences_button'] = Button(self, name='sequences_button', command=self.sequences_button, font=fontStr, text='Add Sequences')
    self._widgets['sequences_button'].grid(column=1, row=2)
    self._widgets['matrices_button'] = Button(self, name='matrices_button', command=self.matrices_button, font=fontStr, text='Add matrices')
    self._widgets['matrices_button'].grid(column=2, row=2)
    self._widgets['getSites_button'] = Button(self, name='getSites_button', command=self.getSites_button, font=fontStr, text='Get Binding Sites',state="disabled")
    self._widgets['getSites_button'].grid(column=1, row=3)
    self._widgets['saveSites_button'] = Button(self, name='saveSites_button', command=self.saveSites_button, font=fontStr, text='Show/save Sites',state="disabled")
    self._widgets['saveSites_button'].grid(column=2, row=3)
    self._widgets['align_button'] = Button(self, name='align_button', command=self.align_button, font=fontStr, text='Align Sites')
    self._widgets['align_button'].grid(column=1, row=4)
    self._widgets['saveAlign_button'] = Button(self, name='saveAlign_button', command=self.saveAlign_button, font=fontStr, text='Show/save Alignments',state="disabled")
    self._widgets['saveAlign_button'].grid(column=2, row=4)
    ## Scroll commands

    ## Resize behavior(s)
    self.grid_rowconfigure(1, weight=0, minsize=162)
    self.grid_rowconfigure(2, weight=1, minsize=50)
    self.grid_rowconfigure(3, weight=1, minsize=50)
    self.grid_rowconfigure(4, weight=1, minsize=50)
    self.grid_columnconfigure(1, weight=0, minsize=50)
    self.grid_columnconfigure(2, weight=0, minsize=50)
    ## Call to post-init method
    self._init_after()


## vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv ##
## vvvv BEGINNING OF CODE WORTH READING; THIS IS STILL AUTO-GENERATED CODE - DO NOT MODIFY YET! vvvv ##
## vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv ##

  ## METHOD _init_before:
  ## ~~~~~~
  def _init_before(self):
    ## The following are the variables used in the interface
    ## Feel free to copy/move them to your _init_specificBefore method if they must be defined otherwise
    ## (End of variables)
    self._init_specificBefore()

  ## METHOD _init_after:
  ## ~~~~~~
  def _init_after(self):
    ## The following are the menus used in the interface
    ## Feel free to copy/move them to your _init_specificAfter method if they must be defined otherwise
    ## (End of menus)
    self._init_specificAfter()


## vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv ##
## vvvv BEGINNING OF USER-CUSTOMIZABLE CODE vvvv ##
## vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv ##
#~CUST_BEGIN (marker; do not delete or edit!)
  ## METHOD _init_specificBefore:
  ## ~~~~~~
  def _init_specificBefore(self):
    pass


  ## METHOD _init_specificAfter:
  ## ~~~~~~
  def _init_specificAfter(self):
    pass




  ## METHOD sequences_button:
  ## ~~~~~~
  def sequences_button(self):
    
    self.seq_w=Toplevel()
    frm=Frame(self.seq_w)
    self._widgets["sequences_pick"]=FileSelectBox(self.seq_w,command=self.addSequences,pattern="*.fa",directory=self.defaultDir)
    self._widgets["sequences_pick"].pack()

    frm2=Frame(frm)

    Button(frm2,text="All Sequences",command=self.addAllSequences).pack(side=LEFT,padx=5,pady=3)
    Button(frm2,text="Done",command=self.seq_w.destroy).pack(side=LEFT,padx=5,pady=3)
    frm2.pack()
    frm.pack()
    pass


  ## METHOD matrices_button:
  ## ~~~~~~
  def matrices_button(self):
    self.mat_w=Toplevel()
    print "self.defaultDir",self.defaultDir
    frm=Frame(self.mat_w)
    self._widgets["mat_pick"]=FileSelectBox(frm,command=self.addMatrices,pattern="*.pfm",directory=self.defaultDir)
    self._widgets["mat_pick"].pack()

    frm2=Frame(frm)

    Button(frm2,text="All matrices",command=self.addAllMatrices).pack(side=LEFT)
    Button(frm2,text="Done",command=self.mat_w.destroy).pack(side=LEFT)
    frm2.pack()
    frm.pack()
    pass









#~CUST_END (marker; do not delete or edit!)
## ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ##
## ^^^^ END OF CUSTOMIZABLE CODE; DO NOT MODIFY BELOW... ^^^^ ##
## ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ##


if __name__ == '__main__':
  root = Tk()
  o = tixgui(root)
  o.pack(side=TOP, fill=BOTH, expand=1)
  root.mainloop()
