import sys
from cStringIO import StringIO


#
# $Log$
#
# For making even sized sequence "logos" which are not sequence logos

# Legaleze for the postscript

##WebLogo License

##Copyright (c) 2002 Gavin E. Crooks, Gary Hon, 
##John-Marc Chandonia and Steven E. Brenner

##Permission is hereby granted, free of charge, to
##any person obtaining a copy of this software and
##associated documentation files (the "Software"),
##to deal in the Software without restriction,
##including without limitation the rights to use,
##copy, modify, merge, publish, distribute,
##sublicense, and/or sell copies of the Software,
##and to permit persons to whom the Software is
##furnished to do so, subject to the following
##conditions:
                               
##The above copyright notice and this permission
##notice shall be included in all copies or
##substantial portions of the Software.
                               
##THE SOFTWARE IS PROVIDED "AS IS", WITHOUT
##WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
##INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
##MERCHANTABILITY, FITNESS FOR A PARTICULAR
##PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
##THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
##ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
##IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
##ARISING FROM, OUT OF OR IN CONNECTION WITH THE
##SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
##SOFTWARE.

##(This is the MIT Open Source License, 
##http://www.opensource.org/licenses/mit-license.html)



class PfmFormat:
    psProlog1="""%%!PS-Adobe-3.0 EPSF-3.0
%%%%Title: Sequence Logo:
%%%%Creator: pfm2ps
%%%%CreationDate: 
%%%%BoundingBox:   0  0  %d  100 
%%%%Pages: 0
%%%%DocumentFonts: 
%%%%EndComments

"""
    psProlog2="""
%% ---- CONSTANTS ----
/cmfactor 72 2.54 div def %% defines points -> cm conversion
/cm {cmfactor mul} bind def %% defines centimeters


%% ---- VARIABLES ----

/black [0 0 0] def
/red [0.8 0 0] def
/green [0 0.8 0] def
/blue [0 0 0.8] def
/yellow [1 0.71 .0] def
/purple [0.8 0 0.8] def
/orange [1 0.7 0] def

/charsPerLine %d def

/logoWidth 18 cm def
/logoHeight 5 cm def
/logoTitle () def

/yaxis false def
/yaxisLabel (bits) def
/yaxisBits  2 def %% bits
/yaxisTicBits 1 def


/xaxis false def
/xaxisLabel () def
/showEnds (d) def %% d: DNA, p: PROTEIN, -: none

/showFineprint false def
/fineprint (U. of Helsinki) def

/logoLines 1 def

/showingBox (n) def    %%n s f
/shrinking false def
/shrink  0.5 def
/outline false def

/IbeamFraction  1 def
/IbeamGray      0.50 def
/IbeamLineWidth 0.5 def

/fontsize       12 def
/titleFontsize  14 def
/smallFontsize   6 def

/defaultColor black def 

/charWidth %d def

"""
    psProlog3="""
% Standard DNA/RNA color scheme
/colorDict << 
   (G)  orange
   (T)  blue   
   (C)  green  
   (A)  red 
   (U)  black   

>> def

% Standard DNA/RNA color scheme
% /colorDict << 
%   (G)  orange
%   (T)  red   
%   (C)  blue  
%   (A)  green 
%   (U)  red   
% >> def

% Standard Amino Acid colors
%/colorDict << 
%  (G)  green  
%  (S)  green  
%  (T)  green  
%  (Y)  green  
%  (C)  green  
%  (N)  purple 
%  (Q)  purple 
%  (K)  blue   
%  (R)  blue   
%  (H)  blue   
%  (D)  red    
%  (E)  red    
%  (P)  black  
%  (A)  black  
%  (W)  black  
%  (F)  black  
%  (L)  black  
%  (I)  black  
%  (M)  black  
%  (V)  black  
%>> def



% ---- DERIVED PARAMETERS ----

/leftMargin
  fontsize 1.5 mul

def 

/bottomMargin
  fontsize 0.75 mul

  % Add extra room for axis
  xaxis {fontsize 1.75 mul add } if
  xaxisLabel () eq {} {fontsize 0.75 mul add} ifelse
def


/topMargin 
  logoTitle () eq { 10 }{titleFontsize 4 add} ifelse
def

/rightMargin 
  %Add extra room if showing ends
  showEnds (-) eq { fontsize}{fontsize 1.5 mul} ifelse
def

/yaxisHeight 
  logoHeight 
  bottomMargin sub  
  topMargin sub
def

/ticWidth fontsize 2 div def

/pointsPerBit yaxisHeight yaxisBits div  def

/isBoxed 
  showingBox (s) eq
  showingBox (f) eq or { 
    true
  } {
    false
  } ifelse
def

/stackMargin 1 def

% Do not add space aroung characters if characters are boxed
/charRightMargin 
  isBoxed { 0.0 } {stackMargin} ifelse
def

/charTopMargin 
  isBoxed { 0.0 } {stackMargin} ifelse
def

%/charWidth
%  logoWidth
%  leftMargin sub
%  rightMargin sub
%  charsPerLine div
%  charRightMargin sub
%def

/charWidth4 charWidth 4 div def
/charWidth2 charWidth 2 div def

/stackWidth 
  charWidth charRightMargin add
def
 
/numberFontsize 
  fontsize charWidth lt {fontsize}{charWidth} ifelse
def

% movements to place 5'/N and 3'/C symbols
/leftEndDeltaX  fontsize neg         def
/leftEndDeltaY  fontsize 1.5 mul neg def
/rightEndDeltaX fontsize 0.25 mul     def
/rightEndDeltaY leftEndDeltaY        def

% Outline width is proporional to charWidth, 
% but no less that 1 point
/outlinewidth 
  charWidth 32 div dup 1 gt  {}{pop 1} ifelse
def


% ---- PROCEDURES ----

/StartLogo { 
  % Save state
  save 
  gsave 

  % Print Logo Title, top center 
  gsave 
    SetTitleFont

    logoWidth 2 div
    logoTitle
    stringwidth pop 2 div sub
    logoHeight logoLines mul  
    titleFontsize sub
    moveto

    logoTitle
    show
  grestore

  % Print X-axis label, bottom center
  gsave
    SetStringFont

    logoWidth 2 div
    xaxisLabel stringwidth pop 2 div sub
    fontsize 3 div
    moveto

    xaxisLabel
    show
  grestore

  % Show Fine Print
  showFineprint {
    gsave
      SetSmallFont
      logoWidth
        fineprint stringwidth pop sub
        smallFontsize sub
          smallFontsize 3 div
      moveto
    
      fineprint show
    grestore
  } if

  % Move to lower left corner of last line, first stack
  leftMargin bottomMargin translate

  % Move above first line ready for StartLine 
  0 logoLines logoHeight mul translate

  SetLogoFont
} bind def

/EndLogo { 
  grestore 
  showpage 
  restore 
} bind def


/StartLine{ 
  % move down to the bottom of the line:
  0 logoHeight neg translate
  
  gsave 
    yaxis { MakeYaxis } if
    xaxis { ShowLeftEnd } if
} bind def

/EndLine{ 
    xaxis { ShowRightEnd } if
  grestore 
} bind def


/MakeYaxis {
  gsave    
    stackMargin neg 0 translate
    ShowYaxisBar
    ShowYaxisLabel
  grestore
} bind def


/ShowYaxisBar { 
  gsave  
    SetStringFont

    /str 10 string def % string to hold number  
    /smallgap stackMargin 2 div def

    % Draw first tic and bar
    gsave    
      ticWidth neg 0 moveto 
      ticWidth 0 rlineto 
      0 yaxisHeight rlineto
      stroke
    grestore

   
    % Draw the tics
    % initial increment limit proc for
    0 yaxisTicBits yaxisBits abs %cvi
    {/loopnumber exch def

      % convert the number coming from the loop to a string
      % and find its width
      loopnumber 10 str cvrs
      /stringnumber exch def % string representing the number

      stringnumber stringwidth pop
      /numberwidth exch def % width of number to show

      /halfnumberheight
         stringnumber CharBoxHeight 2 div
      def

      numberwidth % move back width of number
      neg loopnumber pointsPerBit mul % shift on y axis
      halfnumberheight sub % down half the digit

      moveto % move back the width of the string

      ticWidth neg smallgap sub % Move back a bit more  
      0 rmoveto % move back the width of the tic  

      stringnumber show
      smallgap 0 rmoveto % Make a small gap  

      % now show the tic mark
      0 halfnumberheight rmoveto % shift up again
      ticWidth 0 rlineto
      stroke
    } for
  grestore
} bind def

/ShowYaxisLabel {
  gsave
    SetStringFont

    % How far we move left depends on the size of
    % the tic labels.
    /str 10 string def % string to hold number  
    yaxisBits yaxisTicBits div cvi yaxisTicBits mul 
    str cvs stringwidth pop
    ticWidth 1.5 mul  add neg  


    yaxisHeight
    yaxisLabel stringwidth pop
    sub 2 div

    translate
    90 rotate
    0 0 moveto
    yaxisLabel show
  grestore
} bind def


/StartStack {  % <stackNumber> startstack
  xaxis {MakeNumber}{pop} ifelse
  gsave
} bind def

/EndStack {
  grestore
  stackWidth 0 translate
} bind def


% Draw a character whose height is proportional to symbol bits
/MakeSymbol{ % charbits character MakeSymbol
  gsave
    /char exch def
    /bits exch def

    /bitsHeight 
       bits pointsPerBit mul 
    def

    /charHeight 
       bitsHeight charTopMargin sub
       dup 
       0.0 gt {}{pop 0.0} ifelse % if neg replace with zero 
    def 
 
    charHeight 0.0 gt {
      char SetColor
      charWidth charHeight char ShowChar

      showingBox (s) eq { % Unfilled box
        0 0 charWidth charHeight false ShowBox
      } if

      showingBox (f) eq { % Filled box
        0 0 charWidth charHeight true ShowBox
      } if

    } if

  grestore

  0 bitsHeight translate 
} bind def


/ShowChar { % <width> <height> <char> ShowChar
  gsave
    /tc exch def    % The character
    /ysize exch def % the y size of the character
    /xsize exch def % the x size of the character

    /xmulfactor 1 def 
    /ymulfactor 1 def

  
    % if ysize is negative, make everything upside down!
    ysize 0 lt {
      % put ysize normal in this orientation
      /ysize ysize abs def
      xsize ysize translate
      180 rotate
    } if

    shrinking {
      xsize 1 shrink sub 2 div mul
        ysize 1 shrink sub 2 div mul translate 

      shrink shrink scale
    } if

    % Calculate the font scaling factors
    % Loop twice to catch small correction due to first scaling
    2 {
      gsave
        xmulfactor ymulfactor scale
      
        ysize % desired size of character in points
        tc CharBoxHeight 
        dup 0.0 ne {
          div % factor by which to scale up the character
          /ymulfactor exch def
        } % end if
        {pop pop}
        ifelse

        xsize % desired size of character in points
        tc CharBoxWidth  
        dup 0.0 ne {
          div % factor by which to scale up the character
          /xmulfactor exch def
        } % end if
        {pop pop}
        ifelse
      grestore
    } repeat

    % Adjust horizontal position if the symbol is an I
    tc (I) eq {
      charWidth 2 div % half of requested character width
      tc CharBoxWidth 2 div % half of the actual character
      sub 0 translate
      % Avoid x scaling for I 
      /xmulfactor 1 def 
    } if


    % ---- Finally, draw the character
  
    newpath
    xmulfactor ymulfactor scale

    % Move lower left corner of character to start point
    tc CharBox pop pop % llx lly : Lower left corner
    exch neg exch neg
    moveto

    outline {  % outline characters:
      outlinewidth setlinewidth
      tc true charpath
      gsave 1 setgray fill grestore
      clip stroke
    } { % regular characters
      tc show
    } ifelse

  grestore
} bind def


/ShowBox { % x1 y1 x2 y2 filled ShowBox
  gsave
    /filled exch def 
    /y2 exch def
    /x2 exch def
    /y1 exch def
    /x1 exch def
    newpath
    x1 y1 moveto
    x2 y1 lineto
    x2 y2 lineto
    x1 y2 lineto
    closepath

    clip
    
    filled {
      fill
    }{ 
      0 setgray stroke   
    } ifelse

  grestore
} bind def


/MakeNumber { % number MakeNumber
  gsave
    SetNumberFont
    stackWidth 0 translate
    90 rotate % rotate so the number fits
    dup stringwidth pop % find the length of the number
    neg % prepare for move
    stackMargin sub % Move back a bit
    charWidth (0) CharBoxHeight % height of numbers
    sub 2 div %
    moveto % move back to provide space
    show
  grestore
} bind def


/Ibeam{ % heightInBits Ibeam
  gsave
    % Make an Ibeam of twice the given height in bits
    /height exch  pointsPerBit mul def 
    /heightDRAW height IbeamFraction mul def

    IbeamLineWidth setlinewidth
    IbeamGray setgray 

    charWidth2 height neg translate
    ShowIbar
    newpath
      0 0 moveto
      0 heightDRAW rlineto
    stroke
    newpath
      0 height moveto
      0 height rmoveto
      currentpoint translate
    ShowIbar
    newpath
    0 0 moveto
    0 heightDRAW neg rlineto
    currentpoint translate
    stroke
  grestore
} bind def


/ShowIbar { % make a horizontal bar
  gsave
    newpath
      charWidth4 neg 0 moveto
      charWidth4 0 lineto
    stroke
  grestore
} bind def


/ShowLeftEnd {
  gsave
    SetStringFont
    leftEndDeltaX leftEndDeltaY moveto
    showEnds (d) eq {(5) show ShowPrime} if
    showEnds (p) eq {(N) show} if
  grestore
} bind def


/ShowRightEnd { 
  gsave
    SetStringFont
    rightEndDeltaX rightEndDeltaY moveto
    showEnds (d) eq {(3) show ShowPrime} if
    showEnds (p) eq {(C) show} if
  grestore
} bind def


/ShowPrime {
  gsave
    SetPrimeFont
    (\242) show 
  grestore
} bind def

 
/SetColor{ % <char> SetColor
  dup colorDict exch known {
    colorDict exch get aload pop setrgbcolor
  } {
    pop
    defaultColor aload pop setrgbcolor
  } ifelse 
} bind def

% define fonts
/SetTitleFont {/Times-Bold findfont titleFontsize scalefont setfont} bind def
%/SetLogoFont  {/Helvetica-Narrow-Bold findfont charWidth  scalefont setfont} bind def
/SetLogoFont  {/Times-Bold findfont charWidth  scalefont setfont} bind def
/SetStringFont{/Helvetica-Bold findfont fontsize scalefont setfont} bind def
/SetPrimeFont {/Symbol findfont fontsize scalefont setfont} bind def
/SetSmallFont {/Helvetica findfont smallFontsize scalefont setfont} bind def

/SetNumberFont {
    /Helvetica-Bold findfont 
    numberFontsize
    scalefont
    setfont
} bind def

%Take a single character and return the bounding box
/CharBox { % <char> CharBox <lx> <ly> <ux> <uy>
  gsave
    newpath
    0 0 moveto
    % take the character off the stack and use it here:
    true charpath 
    flattenpath 
    pathbbox % compute bounding box of 1 pt. char => lx ly ux uy
    % the path is here, but toss it away ...
  grestore
} bind def


% The height of a characters bounding box
/CharBoxHeight { % <char> CharBoxHeight <num>
  CharBox
  exch pop sub neg exch pop
} bind def


% The width of a characters bounding box
/CharBoxWidth { % <char> CharBoxHeight <num>
  CharBox
  pop exch pop sub neg 
} bind def


% Deprecated names
/startstack {StartStack} bind  def
/endstack {EndStack}     bind def
/makenumber {MakeNumber} bind def
/numchar { MakeSymbol }  bind def

%%%%EndProlog

%%%%Page: 1 1
StartLogo
StartLine %% line number 1
"""

    psEpilog="""
EndLine
EndLogo

%%EOF
"""
    def __init__(self,fname=None,sortOrder=1):
        self.fname=fname
        self.charwidth=40
        if self.fname:
            self.pfm=self.readpfm(self.fname)

        # From top
        #  0-A,C,G,T
        #  1-Larger on Top
        #  2-Smaller on Top
        self.sortOrder=sortOrder
        

    def readpfm(self,fname):
        return [map(int,x.split()) for x in open(fname).readlines()]

    def makeStacks(self,pfm=None):
        if pfm==None:
            pfm=self.pfm
        self.stacks=len(pfm[0])

        stackCodes=["A","C","G","T"]

        self.stackStrs=[]
        for i in range(self.stacks):
            tot=0
            for j in range(len(pfm)):
                tot=tot+pfm[j][i]
            sStr=[]
            for j in range(len(pfm)):
                if pfm[j][i]>0:
                    height=pfm[j][i]*1.0/tot
                    char=stackCodes[j]
                    sStr.append((height,char," %f (%s) numchar\n"%(height,char)))
            self.stackStrs.append(sStr)
        return self.stackStrs

    def setSort(self,sortOrder):
        self.sortOrder=sortOrder


    def sort(self):
        for i in range(self.stacks):
            if self.sortOrder==0:
                self.stackStrs[i].sort(lambda x,y:cmp(x[1],y[1]))
            else:
                self.stackStrs[i].sort(lambda x,y:cmp((x[0],-ord(x[1])),(y[0],-ord(y[1]))))
                if self.sortOrder==2:
                    self.stackStrs[i].reverse()

    def stacksStr(self,s=StringIO()):
        self.sort()
        for (i,stack) in zip(range(self.stacks),self.stackStrs):
            s.write("(%d) startstack\n"%(i+1))
            for height,char,psString in stack:
                s.write(psString)
            s.write("endstack\n\n")
        return s
                    
    def __str__(self):
        s=StringIO()
        if not hasattr(self,"stacks"):
            self.makeStacks()
        s.write(self.psProlog1%((self.stacks+1)*self.charwidth))
        s.write(self.psProlog2%(self.stacks,self.charwidth))
        s.write(self.psProlog3)
        self.stacksStr(s)
        s.write(self.psEpilog)
        return s.getvalue()
    
                                
if __name__=="__main__":
    if len(sys.argv)<2:
        print "usage: python pfm2ps.py TF/joku.pfm"
        sys.exit(1)
    pfm=PfmFormat(sys.argv[1])
    
    print str(pfm)
    #print "\n".join(makeStacks(pfm))
