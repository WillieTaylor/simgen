%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ws-rv9x6.cls   :   27 January 2003
%% Class file for Review Volume Trim Size [9in x 6in] written in Latex2E.
%%
%% The content, structure, format and layout of this style file is the 
%% property of World Scientific Publishing Co. Pte. Ltd. 
%% Copyright 1995, 2002 by World Scientific Publishing Co. 
%% All rights are reserved.
%%
%% Review Volume Trim Size: 9in x 6in
%% Text Area: 7.35in (include runningheads) x 4.5in
%% Main Text is 10/13pt					  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{ws-rv9x6}[1994/07/13 v1.2u Standard LaTeX document 
class]
\@twosidetrue\@mparswitchtrue\@twocolumnfalse
\newif\if@openright
\@openrighttrue
\newif\if@restonecol
%
\newif\if@draft
%
\newdimen\draftrule
\draftrule=0pt
%
%
\DeclareOption{splissue}{\@splissuetrue}
\DeclareOption{draft}{\@drafttrue\setlength\overfullrule{5pt}\draftrule0.01pt}%

%
\hoffset=0.1in  
\voffset=0.65in

%
\DeclareOption{final}{\setlength\overfullrule{0pt}}
%
\ExecuteOptions{a4paper,10pt,twoside,onecolumn,final}
\ProcessOptions
%
\newcommand\chaptername{Chapter}
\newcommand\contentsname{Contents} 
\newcommand\indexname{Index}  
\newcommand\appendixname{Appendix}
\newcommand\abstractname{Abstract}
\newcommand\bibname{References}%
%

\newcommand\today{\ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
	\space\number\day, \number\year}
\newcount\minute
\newcount\hour
\def\currenttime{%
	\minute\time
	\hour\minute
	\divide\hour60
	\the\hour:\multiply\hour60\advance\minute-\hour\the\minute}
%
\newdimen\trimheight
\newdimen\trimwidth
\newdimen\typeheight
\newdimen\typewidth
\newdimen\tempdimen
\newdimen\tablewidth
\newdimen\normaltextheight
\newbox\tempbox
\newdimen\tablewd
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fonts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\foliofont{\fontsize{8}{10}\selectfont}
\def\bibfont{\fontsize{9}{11}\selectfont}
\def\rhfont{\footnotesize\itshape{}}
\def\CNfont{\fontsize{10}{13}\bfseries\boldmath\selectfont{}}
\def\CTfont{\fontsize{10}{13}\bfseries\boldmath\selectfont{}}
\def\authorfont{\small\centering{}}
\def\abstractfont{\small}
\def\affiliationfont{\small\itshape\centering{}}
\def\datefont{\footnotesize}
\def\sectionfont{\bfseries\boldmath\rightskip2pc plus1fill{}}
\def\subsectionfont{\bfseries\boldmath\itshape\rightskip2pc plus1fill{}}
\def\subsubsectionfont{\itshape\rightskip2pc plus1fill{}}
\def\paragraphfont{\bfseries\boldmath}
\def\subparagraphfont{\itshape}
\def\figurecaptionfont{\footnotesize}
\def\tablefont{\footnotesize}
\def\tablecaptionfont{\footnotesize}

%
\renewcommand\normalsize{%
   \@setfontsize\normalsize\@xpt{13}%\@xiiipt
   \abovedisplayskip6\p@\@plus4\p@\@minus2\p@
   \abovedisplayshortskip\abovedisplayskip
   \belowdisplayshortskip\abovedisplayshortskip
   \belowdisplayskip\abovedisplayskip
   \let\@listi\@listI}
\normalsize
%
\newcommand\small{%
   \@setfontsize\small\@ixpt{11}%
   \abovedisplayskip 8\p@ \@plus3\p@ \@minus4\p@
   \abovedisplayshortskip \z@ \@plus2\p@
   \belowdisplayshortskip 4\p@ \@plus2\p@ \@minus2\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 4\p@ \@plus2\p@ \@minus2\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip}
%
\newcommand\footnotesize{%
   \@setfontsize\footnotesize\@viiipt{10}%
   \abovedisplayskip 6\p@ \@plus2\p@ \@minus4\p@
   \abovedisplayshortskip \z@ \@plus\p@
   \belowdisplayshortskip 3\p@ \@plus\p@ \@minus2\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 3\p@ \@plus\p@ \@minus\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip}
%
\newcommand\scriptsize{\@setfontsize\scriptsize\@viipt\@viiipt}
\newcommand\tiny{\@setfontsize\tiny\@vpt\@vipt}
\newcommand\large{\@setfontsize\large\@xiipt{14}}%
\newcommand\Large{\@setfontsize\Large\@xivpt{16}}%%
\newcommand\LARGE{\@setfontsize\LARGE\@xviipt{22}}
\newcommand\huge{\@setfontsize\huge\@xxpt{25}}
\newcommand\Huge{\@setfontsize\Huge\@xxvpt{30}}
\setlength\smallskipamount{3\p@ \@plus 1\p@ \@minus 1\p@}
\setlength\medskipamount{6\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\bigskipamount{12\p@ \@plus 4\p@ \@minus 4\p@}
%
\DeclareOldFontCommand{\rm}{\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\scshape}{\@nomath\sc}
\DeclareOldFontCommand{\bi}{\bfseries\itshape}{\bfseries\itshape}
\renewcommand{\rmdefault}{cmr}%
%
\DeclareMathSizes{28}{28}{25}{20}
\DeclareMathSizes{24}{24}{21}{18}
\DeclareMathSizes{22}{22}{19}{17}
\DeclareMathSizes{20}{20}{17}{15}
\DeclareMathSizes{18}{18}{15}{10}
\DeclareMathSizes{12}{11}{8}{5}
\DeclareMathSizes{11}{11}{8}{5}
\DeclareMathSizes{10}{10}{7}{5}
\DeclareMathSizes{9.5}{9}{7}{5}
\DeclareMathSizes{9}{9}{7}{5}
\DeclareMathSizes{8}{8}{6}{5}
\DeclareMathSizes{7}{7}{5}{5}
\DeclareMathSizes{6}{6}{5}{5}
\DeclareMathSizes{5}{5}{5}{5}
\def\bm#1{\mathchoice
	{\mbox{\boldmath$#1$}}%
          {\mbox{\boldmath$#1$}}%
          {\mbox{\boldmath$\scriptstyle#1$}}%
          {\mbox{\boldmath$\scriptscriptstyle#1$}}}%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Initializations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

\trimheight8.97truein	%measure is 9in
\trimwidth5.97truein
%
\typeheight7.35in
\typewidth4.5truein
%
\setlength\parindent{1.5em}
\setlength\headheight{5\p@}  
\setlength\headsep   {22pt}%
\setlength\topskip   {7\p@}%
%\setlength\footskip{0.17in}%
\setlength\footskip{0.36in}%
\setlength\maxdepth{4\p@}
\setlength\@maxdepth\maxdepth
\setlength\textwidth{4.5truein}  
\setlength\textheight\typeheight
\advance\textheight-\headheight
\advance\textheight-\headsep
\setlength\normaltextheight{\textheight}
\setlength\oddsidemargin{0.73in}%
\setlength\evensidemargin{\trimwidth}
\addtolength\evensidemargin{-\textwidth}
\addtolength\evensidemargin{-\oddsidemargin}
\setlength\marginparwidth{.75in}
\setlength\marginparsep{7\p@}
\setlength\marginparpush{5\p@}
\setlength\topmargin{0.5in}%
\setlength\footnotesep{6\p@}%6pt
\setlength{\skip\footins}{18\p@ \@plus 5\p@ \@minus 2\p@}
\setlength\floatsep    {18\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\textfloatsep{18\p@ \@plus 2\p@ \@minus 4\p@}
\setlength\intextsep   {18\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\dblfloatsep    {18\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\dbltextfloatsep{18\p@ \@plus 2\p@ \@minus 4\p@}
\setlength\@fptop{0\p@}
\setlength\@fpsep{8\p@ \@plus 2fil}
\setlength\@fpbot{0\p@ \@plus 1fil}
\setlength\@dblfptop{0\p@ \@plus 1fil}
\setlength\@dblfpsep{8\p@ \@plus 2fil}
\setlength\@dblfpbot{0\p@ \@plus 1fil}
\setlength\partopsep{0\p@}
%
\setlength\lineskip{1\p@}
\setlength\normallineskip{1\p@}%
\renewcommand\baselinestretch{}%{}
\setlength\parskip{0\p@}
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Local Definitions %%%%%%%%%%%%%%%%%%%%%%
%

\def\uc#1{\uppercase{#1}}
%\def\cdash{$^{\mbox{\sevenrm --}}$}
\def\cdash{$^{\mbox{\scriptsize --}}$}

\newcommand{\cal}{\protect\pcal}
\newcommand{\pcal}{\@fontswitch{\relax}{\mathcal}}
%
\def\mtb{\fboxsep=0pt\fbox{\phantom{e}}}
\def\breakhere{\par\vskip-\lastskip
			\pagebreak}
\def\forcebreak{\vadjust{\pagebreak}}
%
\def\thebepage{\thepage}
\def\theeepage{\pageref{eof}}
\def\pub#1#2{\par%
	\addvspace{-2pt plus-2pt minus-2pt}
	{\centering{
	\authorfont #1\par
	#2\par}\par}
	\addvspace{15pt plus2pt minus2pt}}%
%

\def\and{\hbox{\lowercase{and }}}
\def\aand{\\[10pt]{\upshape{and}}\\[-3pt]}
\def\@aabuffer{}

%
\def\author#1{%\expandafter\def\expandafter\@aabuffer\expandafter
	{\authorfont{#1}\par}\par
	\vspace*{3pt}}

\def\address#1{%\expandafter\def\expandafter\@aabuffer\expandafter
	{\affiliationfont{#1}\par}\par
	\vspace*{20pt}}

%
\def\keywords#1{\par
	\vspace*{8pt}
	{\authorfont{\leftskip18pt\rightskip\leftskip
	\noindent{\it Keywords}\/:\ #1\par}}\par}
%
\def\category#1#2{\par
	\vspace*{8pt}
	{\authorfont{\leftskip18pt\rightskip\leftskip
	\noindent{\it Categories}\/:\ #1:\ #2\par}}\par}
%
\newenvironment{abstract}{\par%
	\abstractfont
  	\noindent\leftskip18pt\rightskip18pt
}{%
  \par}
%
\def\@citex[#1]#2{\if@filesw\immediate\write\@auxout
	{\string\citation{#2}}\fi	
\def\@citea{}\@cite{\@for\@citeb:=#2\do		
	{\@citea\def\@citea{,}\@ifundefined	
	{b@\@citeb}{{\bf ?}\@warning 
	{Citation `\@citeb' on page \thepage \space undefined}} 
	{\csname b@\@citeb\endcsname}}}{#1}} 
 \newif\if@cghi 
\def\cite{\@cghitrue\@ifnextchar [{\@tempswatrue 
	\@citex}{\@tempswafalse\@citex[]}} 
\def\citelow{\@cghifalse\@ifnextchar [{\@tempswatrue 
	\@citex}{\@tempswafalse\@citex[]}} 
\def\@cite#1#2{{$\null^{#1}$\if@tempswa\typeout 
	{WSPC warning: optional citation argument  
	ignored: `#2'} \fi}} 
\newcommand{\citeup}{\cite} 
% 
\def\@refcitex[#1]#2{\if@filesw\immediate\write\@auxout 
	{\string\citation{#2}}\fi 
\def\@citea{}\@refcite{\@for\@citeb:=#2\do 
	{\@citea\def\@citea{, }\@ifundefined 
	{b@\@citeb}{{\bf ?}\@warning 
	{Citation `\@citeb' on page \thepage \space undefined}} 
	\hbox{\csname b@\@citeb\endcsname}}}{#1}} 
 \def\@refcite#1#2{{#1\if@tempswa\typeout 
        {WSPC warning: optional citation argument 
	ignored: `#2'} \fi}} 
 \def\refcite{\@ifnextchar[{\@tempswatrue
	\@refcitex}{\@tempswafalse\@refcitex[]}}
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Lists %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\@listI{\leftmargin\leftmargini
            \listparindent\itemindent
            \parsep \z@\labelsep.5em
            \topsep 6.5\p@ \@plus3\p@% \@minus2\p@
            \itemsep0\p@}
%
\let\@listi\@listI
\@listi
%
\def\@listii {\leftmargin\leftmarginii
              \labelwidth\leftmarginii
              \advance\labelwidth-\labelsep
              \topsep    3.5\p@ \@plus2\p@ \@minus\p@
              \parsep    \z@
              \itemsep\z@}
\def\@listiii{\leftmargin\leftmarginiii
              \labelwidth\leftmarginiii
              \advance\labelwidth-\labelsep
              \topsep 3.5\p@ \@plus\p@\@minus\p@
              \parsep \z@
             % \partopsep \p@ \@plus\z@ \@minus\p@
              \itemsep0\p@}%
\def\@listiv {\leftmargin\leftmarginiv
              \labelwidth\leftmarginiv
              \advance\labelwidth-\labelsep}
\def\@listv  {\leftmargin\leftmarginv
              \labelwidth\leftmarginv
              \advance\labelwidth-\labelsep}
\def\@listvi {\leftmargin\leftmarginvi
              \labelwidth\leftmarginvi
              \advance\labelwidth-\labelsep}
%
\setlength\leftmargini{3pc}
\setlength\leftmarginii{2.2em}
\setlength\leftmarginiii{1.87em}
\setlength\leftmarginiv{1.7em}
\setlength\leftmarginv{1em}
\setlength\leftmarginvi{1em}
\setlength\leftmargin{\leftmargini}
\setlength\listparindent{\parindent}
\setlength\labelsep{.5em}
\setlength\labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}
\renewcommand\theenumi{\arabic{enumi}}
\renewcommand\theenumii{\alph{enumii}}
\renewcommand\theenumiii{\roman{enumiii}}
\renewcommand\theenumiv{\Alph{enumiv}}
\newcommand\labelenumi{(\theenumi)}
\newcommand\labelenumii{(\theenumii)}
\newcommand\labelenumiii{\theenumiii.}
\newcommand\labelenumiv{\theenumiv.}
\renewcommand\p@enumii{\theenumi}
\renewcommand\p@enumiii{\theenumi(\theenumii)}
\renewcommand\p@enumiv{\p@enumiii\theenumiii}
\newcommand\labelitemi{$\m@th\bullet$}
\newcommand\labelitemii{\normalfont\bfseries --}
\newcommand\labelitemiii{$\m@th\ast$}
\newcommand\labelitemiv{$\m@th\cdot$}
%
\def\enummax#1{\setbox\tempbox=\hbox{#1\hskip\labelsep}%
  \expandafter\global\csname leftmargin\romannumeral\the\@enumdepth\endcsname\wd\tempbox}
%
\def\enumerate{\@ifnextchar[{\@enumerate}{\@enumerate[\csname label\@enumctr\endcsname]}}
%
\def\@enumerate[#1]{\ifnum \@enumdepth >3 \@toodeep\else
      \advance\@enumdepth \@ne\edef\@enumctr{enum\romannumeral\the\@enumdepth}%
      \enummax{#1}\list
      {\csname label\@enumctr\endcsname}{\usecounter
        {\@enumctr}\def\makelabel##1{\hss\llap{##1}}}\fi}
%
\let\Item\item
\newenvironment{enumeroman}{%
  \def\theenumi{\roman{enumi}}\def\theenumii{\alph{enumii}}%
  \def\labelenumi{(\theenumi)}\def\labelenumii{(\theenumii)}%
		\let\item\Item
  \begin{enumerate}%
}{%
  \end{enumerate}}
%
\newenvironment{alphlist}{%
  \def\theenumi{\alph{enumi}}\def\theenumii{\roman{enumii}}%
  \def\labelenumi{(\theenumi)}\def\labelenumii{(\theenumii)}%
		\let\item\Item
  \begin{enumerate}%
}{%
  \end{enumerate}}
%
\newenvironment{arabiclist}{%
	\def\theenumi{\arabic{enumi}}\def\theenumii{\alph{enumii}}
	\def\labelenumi{(\theenumi)}\def\labelenumii{(\theenumii)}%
	\let\item\Item
	\begin{enumerate}
	}{%
	\end{enumerate}}
%
\newenvironment{romanlist}{%
	\def\theenumi{\roman{enumi}}\def\theenumii{\alph{enumii}}
	\def\labelenumi{(\theenumi)}\def\labelenumii{(\theenumii)}%
	\let\item\Item
	\begin{enumerate}
	}{%
	\end{enumerate}}
%
\newenvironment{itemlist}{%
	\def\labelenumi{\labelitemi}
	\let\item\Item
	\begin{enumerate}
	}{%
	\end{enumerate}}
%
\newenvironment{description}
               {\list{}{\labelwidth\z@ \itemindent-\leftmargin
                        \let\makelabel\descriptionlabel}}
               {\endlist}
\newcommand*\descriptionlabel[1]{\hspace\labelsep
                                \normalfont\bfseries #1}
%
\newenvironment{unnumlist}{%
	\let\item\Item
  	\leftmargini2pc
   	\ifnum \@enumdepth >3 \@toodeep\else
      	\advance\@enumdepth \@ne
      	\list{}{\itemindent-2pc\topsep6pt
              \def\makelabel##1{\hss\llap{##1}}}%
   \fi
}{%
  \endlist}
%
\newenvironment{quote}
               {\list{}{\rightmargin18pt\leftmargin18pt}%
                \item[]}
               {\endlist}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\setcounter{secnumdepth}{3}
\newcounter {chapter}
\newcounter {section}
\newcounter {subsection}[section]
\newcounter {subsubsection}[subsection]
\newcounter {paragraph}[subsubsection]
\newcounter {subparagraph}[paragraph]
\renewcommand\thesection      {\arabic{section}}
\renewcommand\thesubsection   {\thesection.\arabic{subsection}}
\renewcommand\thesubsubsection{\thesubsection .\arabic{subsubsection}}
\renewcommand\theparagraph    {\thesubsubsection.\arabic{paragraph}}
\renewcommand\thesubparagraph {\theparagraph.\arabic{subparagraph}}
%
\newcommand\section{\@startsection{section}{1}{\z@}{-18pt \@plus-3pt minus-6pt}{5pt}{\sectionfont}}
\newcommand\subsection{\@startsection{subsection}{2}{\z@}{-18pt \@plus-3pt minus-6pt}{5pt}{\subsectionfont}}
\newcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}{-18pt \@plus-3pt minus-3pt}{5pt}{\subsubsectionfont}}
\newcommand\paragraph{\@startsection{paragraph}{4}{\z@}{12pt \@plus -3pt}{-.5em}{\paragraphfont}}
\newcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}{3.25ex \@plus1ex}{-1em}{\subparagraphfont}}
%
\def\@seccntformat#1{{\upshape{\csname the#1\endcsname.}}\hskip .5em}
%
\gdef\@author{}
%
%
\newcommand*\l@section{\@dottedtocline{1}{\z@}{1.5em}}%1.5em
\newcommand*\l@subsection{\@dottedtocline{2}{1.5em}{2.5em}}%3.8em
\newcommand*\l@subsubsection{\@dottedtocline{3}{4em}{3.3em}}%7.0
\newcommand*\l@paragraph{\@dottedtocline{4}{8.7em}{5em}}%10em
\newcommand*\l@subparagraph{\@dottedtocline{5}{12em}{6em}}
%
%
\newcommand\@pnumwidth{1.55em}
\newcommand\@tocrmarg{2.55em}
\newcommand\@dotsep{4.5}
\setcounter{tocdepth}{3}
\newcommand\tableofcontents{%
{\global
\@topnum\z@ 
\@afterindentfalse 
\if@twocolumn
\@restonecoltrue
\onecolumn 
\else 
\@restonecolfalse 
\fi
\vspace*{24pt}
\noindent 
{\small\bf Contents}\par 
\vskip1em 
\nobreak}
{\small
\@starttoc{toc}%
}\if@restonecol
\twocolumn
\fi}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% chapter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newcommand\chapter{\cleardoublepage
                    \thispagestyle{plain}%
                    \global\@topnum\z@
                    \@afterindentfalse
                    \secdef\@chapter\@schapter}
%
\def\@chapter[#1]#2{\ifnum \c@secnumdepth >\m@ne
			\renewcommand\thepage{\arabic{page}}
                        \else \fi
                    \chaptermark{#1}%
                    \addtocontents{lof}{\protect\addvspace{10\p@}}%
                    \addtocontents{lot}{\protect\addvspace{10\p@}}%
                    \if@twocolumn
                      \@topnewpage[\@makechapterhead{#2}]%
                    \else
                      \@makechapterhead{#2}%
                      \@afterheading
                    \fi}
%
\def\@makechapterhead#1{%
	{\vbox to 115pt{%  
	\refstepcounter{chapter}%
%	\addcontentsline{toc}{chapter}{\chapname#1}
	\vspace*{37pt}
	\parindent\z@\raggedright\reset@font
        {\centering{{\CNfont CHAPTER \thechapter\par}
        \vskip 18\p@
    	\vbox{\CTfont #1}\par}\par}\nobreak\vfill}}}
%
\def\@schapter#1{
	\markboth{#1}{#1}
        \@makeschapterhead{#1}%
        \@afterheading}
%
\def\@makeschapterhead#1{%
%	{\addcontentsline{toc}{schapter}{#1}
	{\vbox to 110pt{%
	\vspace*{33pt}
	\parindent\z@\raggedright\reset@font
        {\centering{{\CNfont \phantom{Chapter \thechapter}\par}
        \vskip 18\p@
    	\vbox{\CTfont #1}\par}\par}\nobreak\vfill}}}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Floats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\setcounter{topnumber}{3}
\def\topfraction{.98}
\setcounter{bottomnumber}{3}
\def\bottomfraction{.98}
\setcounter{totalnumber}{6}%4
\def\textfraction{.01}
\def\floatpagefraction{.98100}
\setcounter{dbltopnumber}{2}
\def\dbltopfraction{.981}
\def\dblfloatpagefraction{.99}
%
%
\newlength\abovecaptionskip
\newlength\belowcaptionskip
\setlength\abovecaptionskip{6\p@}% plus3pt minus3pt}
\setlength\belowcaptionskip{-1\p@}
%
\def\FigName{figure}%
%
\InputIfFileExists{psfig.sty}{\typeout{^^Jpsfig.sty inputed...ok}}{\typeout{^^JWarning: psfig.sty could be be found.^^J}}
\InputIfFileExists{epsfsafe.tex}{\typeout{^^Jepsfsafe.tex inputed...ok}}
			{\typeout{^^JWarning: epsfsafe.tex could not be found.^^J}}
\InputIfFileExists{epsfig.sty}{\typeout{^^Jepsfig.sty inputed...ok}}{\typeout{^^JWarning: epsfig.sty could not be found.^^J}}
\InputIfFileExists{epsf.sty}{\typeout{^^Jepsf.sty inputed...ok}}{\typeout{^^JWarning: epsf.sty could not be found.^^J}}%
\InputIfFileExists{psfig.sty}{}{}
%\InputIfFileExists{epsf.sty}{}{}
%
\newbox\captionbox
%
\long\def\@makecaption#1#2{%
  \ifx\FigName\@captype
    	\vskip\abovecaptionskip
 	\setbox\captionbox\hbox{%
	{\figurecaptionfont#1.\unskip\ignorespaces\hskip1em#2\par}}
	\ifdim\wd\captionbox>\hsize
	{\figurecaptionfont#1.\unskip\ignorespaces\hskip1em#2\par}\par
	\else
	\centerline{\figurecaptionfont#1.\unskip\ignorespaces\hskip1em#2}\par
	\fi\else\global\setbox\captionbox\hbox{%
	{\tablecaptionfont#1.\ignorespaces\hskip1em#2\par}}
	\fi
}
%
\def\figurename{Fig.}
\newcounter{figure}
\renewcommand\thefigure{\@arabic\c@figure}
\def\fps@figure{htbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\figurename~\thefigure}
\newenvironment{figure}
               {\@float{figure}}
               {\end@float}
\newenvironment{figure*}
               {\@dblfloat{figure}}
               {\end@dblfloat}
%
\def\ArtDir{\jobname/}
\def\figurebox#1#2#3{%
        \def\arg{#3}%
        \ifx\arg\empty
        {\hfill\vbox{\hsize#2\hrule\hbox to #2{%
                \vrule\hfill\vbox to #1{\hsize#2\vfill}\vrule}\hrule}\hfill}%
        \else
        \hbox to \hsize{\hfill\epsfbox{\ArtDir#3.eps}\hfill}
        \fi}
%
\def\PSFIG#1#2#3{\par
  {\hbox to \hsize{%
	\hfill{\psfig{file=\ArtDir#1.eps,width=#2,angle=#3}}\hfill}\par}}%
%
\newcounter{table}
\def\tablename{Table}
\renewcommand\thetable{\@arabic\c@table}
\def\fps@table{htbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename~\thetable}
\newenvironment{table}
               {\@float{table}}
               {\end@float}
\newenvironment{table*}
               {\@dblfloat{table}}
               {\end@dblfloat}
%
\long\def\tbl#1#2{%
	\parindent\z@\ignorespaces\noindent\tablecaptionfont
	\caption{#1}%
  	\par\setbox\tempbox\hbox{\tablefont #2}%
  	\tablewd\hsize\advance\tablewd-\wd\tempbox\global\divide\tablewd\tw@
	\ifdim\wd\captionbox<\wd\tempbox\centerline{\unhbox\captionbox}
	\else\leftskip\tablewd\rightskip\leftskip{\unhbox\captionbox}\par
	\fi\vskip5pt\centerline{\box\tempbox}
}%
\long\def\ttbl#1#2#3{%
	\parindent\z@
	\noindent
  	\tablecaptionfont
	\caption{#2}%
	\tablewd\hsize
	\advance\tablewd-#1
  	\global\divide\tablewd\tw@
	\ifdim \wd\captionbox<#1
	\centerline{\unhbox\captionbox}
	\else
  	{\leftskip\tablewd% plus1fill{}
	\rightskip\tablewd% plus1fill{}
  	\unhbox\captionbox\par}
	\fi
 	\vskip-8pt
	\centerline{#3}
}%
%
\newenvironment{InTextTable}{%
 	\par\addvspace{16pt plus6pt minus6pt}%
 	\begin{center}
	\tablecaptionfont\noindent\ignorespaces
}{%
   	\end{center}
	\par\addvspace{17pt plus6pt minus6pt}}
%
\newenvironment{tabnote}{%
	\par\vskip5pt\footnotesize%
	\@ifnextchar[{\@tabnote}{\@tabnote[]}}{% 
	\par\vskip-5pt}
\def\@tabnote[#1]{%
	\leftskip\tablewd\rightskip\leftskip {\it #1} \ignorespaces}
%
\def\x{@{\extracolsep{\fill}}}
\def\toprule{\\[-6pt]\hline\\[-5.5pt]}
\def\colrule{\\[-7.5pt]\hline\\[-5.5pt]}
\def\botrule{\\[-7pt]\hline\\[-8.5pt]}
\def\crule#1{\\[-6pt]\cline{#1}\\[-6.5pt]}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Maths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\setlength\arraycolsep{1.5\p@}
\setlength\tabcolsep{6\p@}
\setlength\arrayrulewidth{.4\p@}
\setlength\doublerulesep{2\p@}
\setlength\tabbingsep{\labelsep}
\skip\@mpfootins = \skip\footins
\setlength\fboxsep{3\p@}
\setlength\fboxrule{.4\p@}
\renewcommand\theequation{\arabic{equation}}
%
\setlength\columnsep{10\p@}
\setlength\columnseprule{0\p@}
%
\def\bstrut{\vrule width0pt depth5pt}
\def\tstrut{\vrule width0pt height10pt}
\jot=3pt
%
%
\renewcommand{\[}{\relax\ifmmode\@badmath\else\bgroup\let\eqnnum\relax\let\refstepeqcnt\relax\@equation\fi}
\renewcommand{\]}{\relax\ifmmode\end@equation\egroup\@endparenv\@doendpe\else\@badmath\fi}
%
%\@namedef{equation*}{\let\eqnnum\relax\let\refstepeqncnt\relax\@equation}
%\@namedef{endequation*}{\end@equation}
\renewenvironment{equation}{\@equation}{\end@equation}
%
\def\refstepeqcnt{\refstepcounter{equation}}
\def\eqnnum{\displaywidth\linewidth\hbox{\@eqnnum}}
%
\def\@equation{\let\topsep\abovedisplayskip\@beginparpenalty\predisplaypenalty\@endparpenalty\postdisplaypenalty
     \refstepeqcnt
     \trivlist \item[]\leavevmode
       \hbox to\linewidth\bgroup\hfill $\m@th% $
       \displaystyle}%
\def\end@equation{$\hfill\eqnnum\egroup\endtrivlist}
%
\renewenvironment{eqnarray}{%
   \stepcounter{equation}%
   \def\@currentlabel{\p@equation\theequation}%
   \global\@eqnswtrue
   \m@th
   \global\@eqcnt\z@
   \tabskip\@centering
   \let\\\@eqncr
   $$\everycr{}\halign to\displaywidth\bgroup
       \hskip\@centering$\displaystyle\tabskip\z@skip{##}$\@eqnsel
      &\global\@eqcnt\@ne\hskip \tw@\arraycolsep \hfil${##}$\hfil
      &\global\@eqcnt\tw@ \hskip \tw@\arraycolsep
         $\displaystyle{##}$\hfil\tabskip\@centering
      &\global\@eqcnt\thr@@ \hb@xt@\z@\bgroup\hss##\egroup
         \tabskip\z@skip
      \cr}
        {%
      \@@eqncr
      \egroup
      \global\advance\c@equation\m@ne
   $$\global\@ignoretrue
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Footnote %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\renewcommand\footnoterule{%
  \kern-3\p@
%  \hrule depth0pt height0pt \@width.2\columnwidth
  \hrule width .4\columnwidth
  \kern2.6\p@}
\newcommand\@makefntext[1]{%
    \parindent 1em%
    \noindent
     \@makefnmark#1}
%
\def\@fnsymbol#1{\ensuremath{\ifcase#1\or *\or \dagger\or \ddagger\or
   	\mathchar "278\or \mathchar "27B\or \|\or **\or \dagger\dagger
   	\or \ddagger\ddagger \else\@ctrerr\fi}}
%
\def\thefootnote{\alph{footnote}}
\def\@makefnmark{{$^{\rm \@thefnmark}$}}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Page Style %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\mark{{}{}}
%
\def\ps@plain{%
     \let\@mkboth\@gobbletwo
     \let\@oddhead\@empty
     %\def\@oddfoot{\hfill\mbox{\lower12pt\hbox{\foliofont\thepage}}\hfill}%
      \def\@oddfoot{\hfill\foliofont\thepage\hfill}%
    \let\@evenhead\@empty
     \let\@evenfoot\@oddfoot}
%
\def\RightMark{}
%
\def\ps@headings{%
    \let\@oddfoot\@empty
    \let\@evenfoot\@empty
%    \def\@evenhead{\rhfont{\upshape{\thepage}}\hskip1em\leftmark\hfill}%
%    \def\@oddhead{\hfill\rhfont\rightmark\hskip1em{\upshape{\thepage}}}%
    \def\@evenhead{\rhfont{\upshape{\thepage}}\hfill\leftmark\hfill}%
    \def\@oddhead{\hfill\rhfont\rightmark\hfill{\upshape{\thepage}}}%
    \let\@mkboth\markboth
    \def\authormark##1{\markboth{##1}{\RightMark}}%
    \def\chaptermark##1{\gdef\RightMark{##1}\markright{##1}}}%
    %\def\sectionmark##1{\gdef\RightMark{##1}\markright{##1}}}%
%
\def\ps@myheadings{%
    \let\@oddfoot\@empty\let\@evenfoot\@empty
    \def\@evenhead{\rlap{\foliofont \thepage}\hfill{\rhfont\leftmark}\hfill}%
    \def\@oddhead{\hfill{\rhfont\rightmark}\hfill{\foliofont\thepage}}%
    \let\@mkboth\@gobbletwo
    \let\chaptermark\@gobble
    \let\sectionmark\@gobble}
%
\def\ps@pblind{%
	\let\@mkboth\@gobbletwo
	\def\@evenhead{}%
	\def\@oddhead{}
	\def\@oddfoot{\hfill\mbox{\lower8pc\hbox{}}\hfill}
   \let\@evenfoot\@oddfoot}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% theorem environments %%%%%%%%%%%%%%%%%%%%
%
\def\@xthm#1#2{\@begintheorem{#2}{\csname the#1\endcsname}{}\ignorespaces}
\def\@ythm#1#2[#3]{\@begintheorem{}{\csname the#1\endcsname}{#3}\ignorespaces}
%
\let\Conventionfont\upshape
\let\Conditionfont\upshape 
\let\Observationfont\upshape
\let\Assumptionfont\itshape
\let\Conjecturefont\itshape
\let\Definitionfont\upshape%
\let\Theoremfont\upshape
\let\Corollaryfont\itshape
\let\Propositionfont\itshape
\let\Lemmafont\itshape
\let\Examplefont\upshape
\let\Theoremheadfont\bfseries
\let\Corollaryheadfont\bfseries
\let\Propositionheadfont\bfseries
\let\Lemmaheadfont\bfseries
\let\Definitionheadfont\bfseries
\let\Remarkfont\upshape
\let\Assumptionheadfont\bfseries
\let\Conditionheadfont\bfseries
\let\Conjectureheadfont\bfseries
\let\Conventionsheadfont\bfseries
\let\Observationheadfont\bfseries
\let\Stepfont\upshape
%
\def\Exampleheadfont{\bfseries}
\def\Remarkheadfont{\bfseries}
\def\Stepheadfont{\bfseries}

%
\def\@begintheorem#1#2#3{\par\addvspace{8pt plus3pt minus2pt}%
              \noindent{\csname#1headfont\endcsname#1\ \ignorespaces#3 #2:}%
              \csname#1font\endcsname\hskip.5em\ignorespaces}
\def\@endtheorem{\par\addvspace{8pt plus3pt minus2pt}\@endparenv}
%
\newcounter{lemma}
\newcounter{corollary}
\newcounter{proposition}
\newcounter{definition}
\newcounter{example}
\newcounter{remark}
\newcounter{question}
\newcounter{notation}

%
%\newtheorem{theorem}{Theorem}%[section]
\newtheorem{theorem}{Theorem}
\newtheorem{corollary}[theorem]{Corollary}%[section]
\newtheorem{lemma}[theorem]{Lemma}
\newtheorem{proposition}[theorem]{Proposition}%[section]
\newtheorem{definition}[theorem]{Definition}%[section]
%\newtheorem{example}[example]{Example}%[section]
\newtheorem{example}[theorem]{Example}%[section]
\newtheorem{remark}[theorem]{Remark}%[section]
\newtheorem{notation}[theorem]{Notation}%[section]
\newtheorem{question}[theorem]{Question}%[section]
\newtheorem{step}{Step}
\newtheorem{assumption}{Assumption}%[section]
\newtheorem{condition}{Condition}%[section]
\newtheorem{conjecture}{Conjecture}%[section]
\newtheorem{conventions}{Conventions}%[section]
\newtheorem{observation}{Observation}%[section]
%

%
\newif\iflogo
\def\prbox{\hbox to \hsize{\hfill\fboxsep0pt\fbox{\phantom{\vrule width5pt height5pt depth0pt}}}\global\logofalse}
%
\newenvironment{proof}{%
\par\addvspace{12pt plus3pt minus3pt}\global\logotrue%
\noindent{\bf Proof:\hskip.5em}\ignorespaces}{%
	\par\iflogo\vskip-\lastskip
	\vskip-\baselineskip\prbox\par
	\addvspace{12pt plus3pt minus3pt}\fi}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Bibliography %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newdimen\bibindent
\bibindent=2em
\def\@biblabel#1{#1.}
\newenvironment{thebibliography}[1]
     {\section*{\bibname}\bibfont%
	\addcontentsline{toc}{section}{\bibname}
      \list{\@biblabel{\arabic{enumiv}}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
	    \labelsep6pt%.33em
            \advance\leftmargin\labelsep
	    \itemsep0pt\parsep0pt
            \usecounter{enumiv}%
            \renewcommand\theenumiv{\arabic{enumiv}}}%
      \sloppy\clubpenalty4000\widowpenalty4000%
      \sfcode`\.=\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}
\newcommand\newblock{}
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% APPENDIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newcounter{appendix}
\newcommand\appendix{\par
  	\refstepcounter{appendix}
	\setcounter{section}{0}%
        \setcounter{lemma}{0}
        \setcounter{theorem}{0}
	\setcounter{definition}{0}
        \setcounter{corollary}{0}
	\setcounter{assumption}{0}
	\setcounter{condition}{0}
	\setcounter{conjecture}{0}
	\setcounter{conventions}{0}
	\setcounter{example}{0}
	\setcounter{observation}{0}
	\setcounter{proposition}{0}
	\setcounter{equation}{0}
	\@addtoreset{equation}{section}
	\renewcommand\thelemma{\Alph{section}.\arabic{lemma}}
	\renewcommand\theremark{\Alph{section}.\arabic{remark}}
	\renewcommand\thetheorem{\Alph{section}.\arabic{theorem}}
	\renewcommand\thecorollary{\Alph{section}.\arabic{corollary}}
	\renewcommand\theassumption{\Alph{section}.\arabic{assumption}}
	\renewcommand\thecondition{\Alph{section}.\arabic{condition}}
	\renewcommand\theconjecture{\Alph{section}.\arabic{conjecture}}
	\renewcommand\theconventions{\Alph{section}.\arabic{conventions}}
	\renewcommand\theexample{\Alph{section}.\arabic{example}}
	\renewcommand\theobservation{\Alph{section}.\arabic{observation}}
	\renewcommand\theproposition{\Alph{section}.\arabic{proposition}}
	\renewcommand\thedefinition{\Alph{section}.\arabic{definition}}
	\renewcommand\thesection{\appendixname\ \Alph{section}}%
	\renewcommand\theequation{\Alph{section}.\arabic{equation}}}%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% APPENDIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (Index) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newenvironment{theindex}
	{\raggedright\small
	\if@twocolumn
         	\@restonecolfalse
         \else
         	\@restonecoltrue
         \fi
         \columnseprule \z@
         \columnsep 18\p@
                \twocolumn[\@makeschapterhead{\MakeUppercase{\indexname}}\vspace*{4pt}]%
                \@mkboth{\indexname}{\indexname}%
                \thispagestyle{plain}\parindent\z@
                \parskip\z@ \@plus .3\p@\relax
                \let\item\@idxitem}
               {\if@restonecol\onecolumn\else\clearpage\fi}
\newcommand\@idxitem  {\par\hangindent 1em}
\newcommand\subitem   {\par\hangindent 40\p@ \hspace*{10\p@}}
\newcommand\subsubitem{\par\hangindent 40\p@ \hspace*{30\p@}}
\newcommand\indexspace{\par\vskip\baselineskip}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (Index) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Trim Marks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\draftnote{\hfill\hbox to \trimwidth{\footnotesize\today\quad\currenttime
\qquad\quad WSPC/Trim Size: 9in x 6in for Review 
Volume\hfill\jobname\hfill}}%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Crop Marks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\cropmarks{\nointerlineskip%
     \vbox to 0pt{\vskip-8.6pt
	\parindent0pt\hspace*{6pt}\draftnote
     \par\advance\trimheight by 2pc
     \vbox to \trimheight{\offinterlineskip\parindent0pt\vskip-2pc
     \hbox to \trimwidth{\raise2.6pc\hbox{\vrule height2pc width0.3pt}\raisebox{2pc}
   {\hspace*{-2.6pc}\vrule width2pc height0.3pt}\hfill\raisebox{2pc}
   {\vrule width2pc height0.3pt}\hspace*{-2.6pc}\raise2.6pc\hbox{\vrule height2pc width0.3pt}}\vfill
     \hbox to \trimwidth{\lower1.2pc\hbox{\vrule height2pc width0.3pt}\raisebox{1.4pc}
        {\hspace*{-2.6pc}\vrule width2pc height0.3pt}\hfill 
    \lower1.2pc\hbox{\vrule height2pc width0.3pt}\rlap{\hspace*{0.6pc}\raisebox{1.4pc}{\vrule width2pc height0.3pt}}}}}
     \vbox to 0pt{\offinterlineskip\parindent0pt\vskip\topmargin%\vskip1.3pt
     \tempdimen\normaltextheight
     \advance\tempdimen27pt  
     \moveright\@themargin
     \vbox to \tempdimen{\hrule height\draftrule 
     \vbox to 0pt{\vskip\headheight\vskip\headsep%\vskip3pt
                \vrule height\draftrule width\textwidth}%
      \hbox to \textwidth{\vrule width\draftrule \hfill\vrule height\tempdimen width\draftrule}%
           \hrule height\draftrule}}}%%
\def\xcropmarks{\nointerlineskip%
     \vbox to 0pt{\vskip-8.6pt
	\parindent0pt\hspace*{6pt}\draftnote
     \par\advance\trimheight by 2pc
     \vbox to \trimheight{\offinterlineskip\parindent0pt\vskip-2pc
     \hbox to \trimwidth{\raise2.6pc\hbox{\vrule height2pc width0.3pt}\raisebox{2pc}
   {\hspace*{-2.6pc}\vrule width2pc height0.3pt}\hfill\raisebox{2pc}
   {\vrule width2pc height0.3pt}\hspace*{-2.6pc}\raise2.6pc\hbox{\vrule height2pc width0.3pt}}\vfill
     \hbox to \trimwidth{\lower1.2pc\hbox{\vrule height2pc width0.3pt}\raisebox{1.4pc}
        {\hspace*{-2.6pc}\vrule width2pc height0.3pt}\hfill 
    \lower1.2pc\hbox{\vrule height2pc width0.3pt}\rlap{\hspace*{0.6pc}\raisebox{1.4pc}{\vrule width2pc height0.3pt}}}}}
     \vbox to 0pt{\offinterlineskip\parindent0pt\vskip\topmargin%\vskip1.3pt
     \tempdimen\normaltextheight
     \advance\tempdimen27pt
     \moveright\@themargin
     \vbox to \tempdimen{\hrule height\draftrule 
     \vbox to 0pt{\vskip\headheight\vskip\headsep%\vskip3pt
                \vrule height\draftrule width\textwidth}%
      \hbox to \textwidth{\vrule width\draftrule \hfill\vrule height\tempdimen width\draftrule}%
           \hrule height\draftrule}}}%%
%
\let\trimmarks\cropmarks
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Crop Marks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output Routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
\def\@outputpage{%
\begingroup           % the \endgroup is put in by \aftergroup
  \let \protect \noexpand
  \@resetactivechars
  \let\-\@dischyph
  \let\'\@acci\let\`\@accii\let\=\@acciii
  \let\\\@normalcr
  \let\par\@@par          %% 15 Sep 87
  \shipout \vbox{%
    \set@typeset@protect
    \aftergroup \endgroup
    \aftergroup \set@typeset@protect
                                % correct? or just restore by ending
                                % the group?
  \if@specialpage
    \global\@specialpagefalse\@nameuse{ps@\@specialstyle}%
  \fi
  \if@twoside
    \ifodd\count\z@ \let\@thehead\@oddhead \let\@thefoot\@oddfoot
         \let\@themargin\oddsidemargin
    \else \let\@thehead\@evenhead
       \let\@thefoot\@evenfoot \let\@themargin\evensidemargin
    \fi
  \fi
  \reset@font
  \normalsize
  \baselineskip\z@skip \lineskip\z@skip \lineskiplimit\z@
    \@begindvi\trimmarks
    \vskip \topmargin
    \moveright\@themargin \vbox {%
      \setbox\@tempboxa \vbox to\headheight{%
        \vfil
        \color@hbox
          \normalcolor
          \hb@xt@\textwidth {%
            \let \label \@gobble
            \let \index \@gobble
            \let \glossary \@gobble %% 21 Jun 91
            \@thehead
            }%
        \color@endbox
        }%                        %% 22 Feb 87
      \dp\@tempboxa \z@
      \box\@tempboxa
      \vskip \headsep
      \box\@outputbox
      \baselineskip \footskip
      \color@hbox
        \normalcolor
        \hb@xt@\textwidth{%
          \let \label \@gobble
          \let \index \@gobble      %% 22 Feb 87
          \let \glossary \@gobble   %% 21 Jun 91
          \@thefoot
          }%
      \color@endbox
      }%
    }%
\global \@colht \textheight
\stepcounter{page}%
\let\firstmark\botmark
}
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Process Options %%%%%%%%%%%%%%%%%%%%%%%%
%
\pagestyle{headings}
\pagenumbering{arabic} 
\frenchspacing
\flushbottom
\def\sloppy{\tolerance=100000\hfuzz=\maxdimen\vfuzz=\maxdimen} 
\vbadness=12000
\hbadness=12000
\sloppy
\usepackage{rotating_rv} 
\usepackage{amsfonts} 
\usepackage{amssymb,amsmath} 
\onecolumn
\endinput
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


