# Version v07  2017-02-08
# by Mathias Zechmeister
# this script allows to scroll and zoom with keyboard keys
# (wxt does not support KP_Add and KP_Subtract, i.e Num + and Num -, and ctrl+char bindings)
# bindings can be reset with:  bind! or reset bind

DELTA_X(a)=a*(GPVAL_X_MAX-GPVAL_X_MIN)
DELTA_Y(a)=a*(GPVAL_Y_MAX-GPVAL_Y_MIN)
SCR_X=.25     # x scroll fraction
SCR_Y=.25     # y scroll fraction
ZOOM_X=0.05   # x zoom factor
ZOOM_Y=0.05   # y zoom factor
ZOOM_XX=2.    # fast x zoom factor
ZOOM_YY=2.    # fast y zoom factor
X_NEW(a) =  0.5*GPVAL_X_MIN*(1-a) + 0.5*GPVAL_X_MAX*(1+a) # = (x1+x2)/2 - a*(x2-x1)/2
Y_NEW(a) =  0.5*GPVAL_Y_MIN*(1-a) + 0.5*GPVAL_Y_MAX*(1+a)

# SCROLLING
bind "Right"       "set xrange [GPVAL_X_MIN+DELTA_X(SCR_X):GPVAL_X_MAX+DELTA_X(SCR_X)]; replot" # scroll right
bind "Left"        "set xrange [GPVAL_X_MIN-DELTA_X(SCR_X):GPVAL_X_MAX-DELTA_X(SCR_X)]; replot" # scroll left
bind "Up"          "set yrange [GPVAL_Y_MIN+DELTA_Y(SCR_Y):GPVAL_Y_MAX+DELTA_Y(SCR_Y)]; replot" # scroll up
bind "Down"        "set yrange [GPVAL_Y_MIN-DELTA_Y(SCR_Y):GPVAL_Y_MAX-DELTA_Y(SCR_Y)]; replot" # scroll down
# SCROLL PAGEWISE
bind "PageUp"      "set yrange [GPVAL_Y_MIN+DELTA_Y(1.)   :GPVAL_Y_MAX+DELTA_Y(1.)];    replot" # page up
bind "Alt-Up"      "set yrange [GPVAL_Y_MIN+DELTA_Y(1.)   :GPVAL_Y_MAX+DELTA_Y(1.)];    replot" # page up
bind "PageDown"    "set yrange [GPVAL_Y_MIN-DELTA_Y(1.)   :GPVAL_Y_MAX-DELTA_Y(1.)];    replot" # page down
bind "Alt-Down"    "set yrange [GPVAL_Y_MIN-DELTA_Y(1.)   :GPVAL_Y_MAX-DELTA_Y(1.)];    replot" # page down
bind "Home"        "set xrange [GPVAL_X_MIN-DELTA_X(1.)   :GPVAL_X_MAX-DELTA_X(1.)];    replot" # page left
bind "Alt-Left"    "set xrange [GPVAL_X_MIN-DELTA_X(1.)   :GPVAL_X_MAX-DELTA_X(1.)];    replot" # page left
bind "End"         "set xrange [GPVAL_X_MIN+DELTA_X(1.)   :GPVAL_X_MAX+DELTA_X(1.)];    replot" # page right
bind "Alt-Right"   "set xrange [GPVAL_X_MIN+DELTA_X(1.)   :GPVAL_X_MAX+DELTA_X(1.)];    replot" # page right
# SCROLL END
bind "Ctrl-Home"   "set xrange [GPVAL_DATA_X_MIN : GPVAL_DATA_X_MIN+DELTA_X(1.)];    replot" # horizontal scroll to first data point
bind "Ctrl-End"    "set xrange [GPVAL_DATA_X_MAX-DELTA_X(1.) : GPVAL_DATA_X_MAX];    replot" # horizontal scroll to last data point
bind "Ctrl-PageUp" "set yrange [GPVAL_DATA_Y_MAX-DELTA_Y(1.) : GPVAL_DATA_Y_MAX];    replot" # scroll to top
bind "Ctrl-PageDown" "set yrange [GPVAL_DATA_Y_MIN : GPVAL_DATA_Y_MIN+DELTA_Y(1.)];  replot" # scroll to bottom

# ZOOMING
bind "Ctrl-Right"  "set xrange [GPVAL_X_MIN+DELTA_X(ZOOM_X):GPVAL_X_MAX-DELTA_X(ZOOM_X)]; replot" # zoom in x
bind "Ctrl-Left"   "set xrange [GPVAL_X_MIN-DELTA_X(ZOOM_X):GPVAL_X_MAX+DELTA_X(ZOOM_X)]; replot" # zoom out x
bind "Ctrl-Up"     "set yrange [GPVAL_Y_MIN-DELTA_Y(ZOOM_Y):GPVAL_Y_MAX+DELTA_Y(ZOOM_Y)]; replot" # zoom out y
bind "Ctrl-Down"   "set yrange [GPVAL_Y_MIN+DELTA_Y(ZOOM_Y):GPVAL_Y_MAX-DELTA_Y(ZOOM_Y)]; replot" # zoom in y
bind "Ctrl-Alt-Right"  "set xrange [X_NEW(-1./ZOOM_XX):X_NEW(1./ZOOM_XX)]; replot" # fast zoom in x
bind "Ctrl-Alt-Left"   "set xrange [X_NEW(-ZOOM_XX):X_NEW(ZOOM_XX)]; replot"       # fast zoom out x
bind "Ctrl-Alt-Up"     "set yrange [Y_NEW(-ZOOM_YY):Y_NEW(ZOOM_YY)]; replot"       # fast zoom out y
bind "Ctrl-Alt-Down"   "set yrange [Y_NEW(-1./ZOOM_YY):Y_NEW(1./ZOOM_YY)]; replot" # fast zoom in y
bind "KP_Add"      "print 'a'; set xrange [X_NEW(-1./ZOOM_XX):X_NEW(1./ZOOM_XX)];\
                    set yrange [GPVAL_Y_MIN+DELTA_Y(ZOOM_Y):GPVAL_Y_MAX-DELTA_Y(ZOOM_Y)]; replot" # zoom in both
bind "+"           "set xrange [GPVAL_X_MIN+DELTA_X(ZOOM_X):GPVAL_X_MAX-DELTA_X(ZOOM_X)];\
                    set yrange [GPVAL_Y_MIN+DELTA_Y(ZOOM_Y):GPVAL_Y_MAX-DELTA_Y(ZOOM_Y)]; replot" # zoom in both
bind "KP_Subtract" "set xrange [GPVAL_X_MIN-DELTA_X(ZOOM_X):GPVAL_X_MAX+DELTA_X(ZOOM_X)];\
                    set yrange [GPVAL_Y_MIN-DELTA_Y(ZOOM_Y):GPVAL_Y_MAX+DELTA_Y(ZOOM_Y)]; replot" # zoom out both
bind "-"           "set xrange [GPVAL_X_MIN-DELTA_X(ZOOM_X):GPVAL_X_MAX+DELTA_X(ZOOM_X)];\
                    set yrange [GPVAL_Y_MIN-DELTA_Y(ZOOM_Y):GPVAL_Y_MAX+DELTA_Y(ZOOM_Y)]; replot" # zoom out both
bind "Alt-KP_Add"  "set xrange [Y_NEW(-1./ZOOM_YY):Y_NEW(1./ZOOM_YY)];\
                    set yrange [Y_NEW(-1./ZOOM_YY):Y_NEW(1./ZOOM_YY)]; replot" # fast zoom in both
bind "Alt-+"       "set xrange [X_NEW(-1./ZOOM_XX):X_NEW(1./ZOOM_XX)];\
                    set yrange [Y_NEW(-1./ZOOM_YY):Y_NEW(1./ZOOM_YY)]; replot" # fast zoom in both
bind "Alt-KP_Subtract" "set xrange [X_NEW(-ZOOM_XX):X_NEW(ZOOM_XX)];\
                    set yrange [Y_NEW(-ZOOM_YY):Y_NEW(ZOOM_YY)]; replot" # fast zoom out both
bind "Alt--"       "set xrange [X_NEW(-ZOOM_XX):X_NEW(ZOOM_XX)];\
                    set yrange [Y_NEW(-ZOOM_YY):Y_NEW(ZOOM_YY)]; replot" # fast zoom out both

bind "C"      "set cbrange [*:*]; replot"                   # reset color range
bind "Alt-C"  "set cbrange [GPVAL_CB_MIN:GPVAL_CB_MAX];"    # save  color range (resetted by 'u', use 'U')
bind "Alt-c"  "set cbrange [GPVAL_CB_MIN:GPVAL_CB_MAX];"    # save  color range (resetted by 'u', use 'U')

# resize all
if (GPVAL_VERSION>4.2)\
bind "Ctrl-u" "set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX];\
               set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]; replot";\
bind "U"      "set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX];\
               set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]; replot;"  # reset ranges to max/min data
bind "Alt-u"  "set xrange [*:*]; set yrange [*:*]; replot"  # total reset
bind "x"      "set xrange [*:*]; replot"
bind "y"      "set yrange [*:*]; replot"
bind "z"      "set yrange [0:]; replot"; # zero yaxis

if (GPVAL_VERSION>4.2)\
bind "Ctrl-x" "set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]; replot";\
bind "X"      "set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]; replot";\
bind "Ctrl-y" "set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]; replot";\
bind "Y"      "set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]; replot";

 # reset ranges to min/max y data

# pan with mouse (press a, then )
bind "a"      "pause mouse; DELTAX=MOUSE_X;        DELTAY=MOUSE_Y;\
               pause mouse; DELTAX=DELTAX-MOUSE_X; DELTAY=DELTAY-MOUSE_Y;\
               set xrange [GPVAL_X_MIN+DELTAX:GPVAL_X_MAX+DELTAX];\
               set yrange [GPVAL_Y_MIN+DELTAY:GPVAL_Y_MAX+DELTAY]; replot"

# center to mouse (if (GPVAL_VERSION<=4.2) pause mouse;?)
bind "c"      "DELTAX=0.5*(GPVAL_X_MAX-GPVAL_X_MIN);\
               DELTAY=0.5*(GPVAL_Y_MAX-GPVAL_Y_MIN);\
               set xrange [MOUSE_X-DELTAX:MOUSE_X+DELTAX];\
               set yrange [MOUSE_Y-DELTAY:MOUSE_Y+DELTAY]; replot" # center to a previous mouse click

bind "G"      'set grid front lt GR_COL=exists("GR_COL")?(GR_COL+1):1; replot'  # toggle various grid colors
bind "Alt-G"  'set grid front lt GR_COL=exists("GR_COL")?(GR_COL-1):1; replot'  # toggle previous grid colors
bind "Alt-g"  'set grid front lt GR_COL=exists("GR_COL")?(GR_COL-1):1; replot'  # toggle previous grid colors

if (GPVAL_VERSION>4.2)\
bind "P"      'PAL_MOD=exists("PAL_MOD")?(PAL_MOD+1)%5:1;\
               eval "set palette col model ".(PAL_MOD==0?"RGB":PAL_MOD==1?"HSV":PAL_MOD==2?"CMY":PAL_MOD==3?"YIQ":"XYZ"); replot'
               # toggling various palette models

bind "B"      'B_FUNC=exists("B_FUNC")?(B_FUNC+1)%36:1; set palette rgb 3,2,B_FUNC; replot'  # toggle blue
bind "Alt-b"  'B_FUNC=exists("B_FUNC")?(B_FUNC-1)%36:1; set palette rgb 3,2,B_FUNC; replot'  # toggle blue
bind "R"      'R_FUNC=exists("R_FUNC")?(R_FUNC+1)%36:1; set palette rgb R_FUNC,3,2; replot'  # toggle red
bind "Alt-r"  'R_FUNC=exists("R_FUNC")?(R_FUNC-1)%36:1; set palette rgb R_FUNC,3,2; replot'  # toggle red

# bind "7"     is built-in toggle ratio
bind "Alt-7"  "set size square; GP_Y_CEN=0.5*(GPVAL_Y_MAX+GPVAL_Y_MIN);\
               set yrange [GP_Y_CEN-DELTA_X(.5):GP_Y_CEN+DELTA_X(.5)]; replot" # set square; keep x fixed
bind "/"      "set size square; GP_X_CEN=0.5*(GPVAL_X_MAX+GPVAL_X_MIN);\
               set xrange [GP_X_CEN-DELTA_Y(.5):GP_X_CEN+DELTA_Y(.5)]; replot" # set square; keep y fixed
               # assumes "Shift-7" == "/" German keyboard

bind "M"      'GP_CBOX=exists("GP_CBOX")?(GP_CBOX+1)%2:1; if (GP_CBOX) unset colorbox; replot;\
               else set colorbox; replot'   # toggle colorbox

bind "Alt-R"  "reset; replot" # reset with "Shift-Alt-r"

# sed -r  -n "/\\\\/N; /^bind /{ s/[^d] [Q\"][^#]*[#]?/  /; s/bind //p} " ~/zoom.gnu
bind "H"      'system("sed -r -n  \"/\\\\\\\\/N; /^bind /{ s/[^d] [\\x27\\\"][^#]*[#]?/  /; s/bind //p} \"  ~/zoom.gnu")'
bind "H"      'system("sed -r -n  \"/\\\\\\\\/N; /^bind /{ s/bind \\\"/ /; s/\\\"[^#]*[#]?/\t\t/p} \"  ~/zoom.gnu")'
# append next line if \ at line end

# awk  '/^bind /{printf("%-16s%s\n",$2,$NF""substr(RT,1,1))}' FS='bind "|" [^#]*|# *' RS='[^\\\\]\n'   zoom.gnu
# awk  "/^bind /{printf(\"%-16s%s\\n\",\$2,\$NF\"\"substr(RT,1,1))}" FS='bind "|" [^#]*|# *' RS='[^\\\\]\n'   zoom.gnu
# bind "H"      'print("awk  \x27/^bind /{printf(\"%-16s%s\\n\",$2,$NF\"\"substr(RT,1,1))}\x27 FS=\x27bind \"|\" [^#]*|# *\x27 RS=\x27[^\\\\\\\\]\n\x27   ~/zoom.gnu")'

#awk  "/^bind /{printf(\"%-16s%s\",\$2,\$NF\"\"RT)}" FS='bind \"|\" [^#]*|# *' RS='[^\\\\]\n'   zoom.gnu
bind "H"      'system("awk  \"/^bind /{printf(\\\"%-16s%s\\\",\\\$2,\\\$NF\\\"\\\"RT)}\" FS=\"bind \\\"|\\\" [^#]*|# *\" RS=\"[^\\\\\\\\\\\\\\]\\n\"   ~/zoom.gnu")' # polished help
bind "H"      'system("awk  \047/^bind /{printf(\"%-16s%s%s\",\$2,\$NF,RT)}\047 FS=\047bind \"|\" [^#]*|# *\047 RS=\047[^\\\\\\\\]\n\047   ~/zoom.gnu")' # polished help

#                      awk 'NR>1{printf("%-16s%s\n",$1,$2)}' RS='\nbind "' FS='" [^#]*|\n' zoom.gnu
