/* A Bison parser, made by GNU Bison 3.0.2.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2013 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.0.2"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* Copy the first part of user declarations.  */
#line 11 "gram.y" /* yacc.c:339  */

/*
 * 
 * evaluate expressions, commands, parameter files
 * 
 */

#define GRAMMAR

#ifndef lint
static char RCSid[] = "$Id: gram.y,v 1.31 2008/04/10 18:29:08 pturner Exp $";
#endif

#include <stdio.h>
#include <math.h>
#include <ctype.h>

#include "defines.h"
#include "globals.h"

typedef struct _symtab_entry {
    char *s;
    int type;
} symtab_entry;

#ifndef M_PI
#     define M_PI  3.14159265358979323846
#endif

#ifndef TRUE
#     define TRUE 1
#endif

#ifndef FALSE
#     define FALSE 0
#endif
/* for LINUX
char *gettxt(char *t, char *s) { return *s; }
*/

int elcirc_maxlevels;
int elcirc_gridno;
int elcirc_flowno;
int elcircmarker = 0;
int curg = 0;

int maxboxes = MAXBOXES;
int maxlines = MAXLINES;
int maxstring = MAXSTR;

double result, resx, resy;	/* return value if expression */
double nonl_parms[10];

double drand48(void);
long lrand48(void);

double rnorm(double mean, double sdev), fx(double x), normp(double b, double *s);
void yyerror(char *s);

static int interr;

static double *freelist[100]; 	/* temporary vectors */
static int fcnt;		/* number allocated */

int naxis = 0;	/* current axis */
int curline, curbox, curstring, curleg, curobject;

int gotbatch, gotparams, gotread; /* these guys attempt to avoid reentrancy problems */
int readtype, readsrc;
extern char batchfile[];
char paramfile[256], readfile[256];

static char f_string[512];	/* buffer for string to parse */
static int pos = 0;
static double *aa, *bb, *cc, *dd, *xx, *yy;
static int setindex, lxy, ls;
static int setsetno;
static int whichgraph;
static int whichset;

extern int change_gno;
extern int change_type;

int checkptr(void *ptr, char *buf);
static Isolparms *setisol = NULL; /* pointer to current Isolparms struct */
static DisplayFlow *setflow = NULL; /* pointer to current Isolparms struct */
static Props *setprops = NULL; /* pointer to current Isolparms struct */
static Zoom_box *setzoombox = NULL;
static Hist_marker *sethistbox = NULL;
static DisplaySlice *setslice = NULL;
static Elevmarker *setelevmarker = NULL;
static DisplayGrid *setgrid = NULL;
static Transect *settrans = NULL;
static ADCIRC3D *setadc3d = NULL;
static Display3dFlow *setflow3d = NULL;
static DisplayParticles *setdrogs = NULL;

/* may add these later TODO
*/


#line 168 "y.tab.c" /* yacc.c:339  */

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "y.tab.h".  */
#ifndef YY_YY_Y_TAB_H_INCLUDED
# define YY_YY_Y_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    VAR = 258,
    X = 259,
    Y = 260,
    CHRSTR = 261,
    FITPARM = 262,
    NUMBER = 263,
    ABS = 264,
    ACOS = 265,
    ASIN = 266,
    ATAN = 267,
    ATAN2 = 268,
    CEIL = 269,
    COS = 270,
    DEG = 271,
    DX = 272,
    DY = 273,
    ERF = 274,
    ERFC = 275,
    EXP = 276,
    FLOOR = 277,
    HYPOT = 278,
    INDEX = 279,
    INT = 280,
    IRAND = 281,
    LGAMMA = 282,
    LN = 283,
    LOG = 284,
    LOGISTIC = 285,
    MAXP = 286,
    MINP = 287,
    MINMAX = 288,
    MOD = 289,
    NORM = 290,
    NORMP = 291,
    PI = 292,
    RAD = 293,
    RAND = 294,
    SETNO = 295,
    SIN = 296,
    SQR = 297,
    SQRT = 298,
    TAN = 299,
    INUM = 300,
    ABORT = 301,
    ABOVE = 302,
    ABSOLUTE = 303,
    ACTIVATE = 304,
    ACTIVE = 305,
    ADCIRC = 306,
    ADCIRC3DFLOW = 307,
    ALL = 308,
    ALT = 309,
    ALTERNATE = 310,
    ALTXAXIS = 311,
    ALTYAXIS = 312,
    AMP = 313,
    ANGLE = 314,
    ANNOTATE = 315,
    APPEND = 316,
    AREA = 317,
    ARROW = 318,
    ASCEND = 319,
    AT = 320,
    ATTACH = 321,
    AUTO = 322,
    AUTOSCALE = 323,
    AUTOTICKS = 324,
    AVERAGE = 325,
    AVG = 326,
    AXES = 327,
    AXIS = 328,
    BACKBUFFER = 329,
    BACKGROUND = 330,
    BAR = 331,
    BATCH = 332,
    BATH = 333,
    BATHYMETRY = 334,
    COURANT = 335,
    BELOW = 336,
    BIN = 337,
    BINARY = 338,
    BOTH = 339,
    BOTTOM = 340,
    BOUNDARY = 341,
    BOX = 342,
    CELLS = 343,
    CENTER = 344,
    CH3D = 345,
    CHAR = 346,
    CHDIR = 347,
    CIRCLE = 348,
    CLEAR = 349,
    CLICK = 350,
    CLOCK = 351,
    CLOSE = 352,
    CM = 353,
    CMAP = 354,
    COLOR = 355,
    COLORMAP = 356,
    COMMENT = 357,
    CONC = 358,
    CONCENTRATION = 359,
    CONCENTRATIONS = 360,
    COPY = 361,
    CROSS = 362,
    CYCLE = 363,
    DAYMONTH = 364,
    DAYOFWEEKL = 365,
    DAYOFWEEKS = 366,
    DAYOFYEAR = 367,
    DAYS = 368,
    DDMMYY = 369,
    DDMONTHSYYHHMMSS = 370,
    DECIMAL = 371,
    DEF = 372,
    DEFAULT = 373,
    DEGREESLAT = 374,
    DEGREESLON = 375,
    DEGREESMMLAT = 376,
    DEGREESMMLON = 377,
    DEGREESMMSSLAT = 378,
    DEGREESMMSSLON = 379,
    DELAYP = 380,
    DELETE = 381,
    DEPTH = 382,
    DEPTHS = 383,
    DESCEND = 384,
    DEVICE = 385,
    DEVXY = 386,
    DFT = 387,
    DT = 388,
    DIAMOND = 389,
    DIFFERENCE = 390,
    DISK = 391,
    DISPLAY = 392,
    DOT = 393,
    DOUBLEBUFFER = 394,
    DOWN = 395,
    DRAW2 = 396,
    DROGUE = 397,
    DROGUES = 398,
    DRY = 399,
    DXDX = 400,
    DXP = 401,
    DYDY = 402,
    DYP = 403,
    ECHO = 404,
    EDIT = 405,
    ELA = 406,
    ELCIRC = 407,
    ELEMENT = 408,
    ELEMENTS = 409,
    ELEV = 410,
    ELEVATION = 411,
    ELEVATIONS = 412,
    ELEVMARKER = 413,
    ELLIPSE = 414,
    ELLIPSES = 415,
    ELLIPSEZ = 416,
    ELSE = 417,
    END = 418,
    ERRORBAR = 419,
    EXIT = 420,
    EXPAND = 421,
    EXPONENTIAL = 422,
    FACTOR = 423,
    FALSEP = 424,
    FAST = 425,
    FEET = 426,
    FFT = 427,
    FILEP = 428,
    FILL = 429,
    FIND = 430,
    FIXEDPOINT = 431,
    FLOW = 432,
    FLUSH = 433,
    FLUX = 434,
    FOCUS = 435,
    FOLLOWS = 436,
    FONTP = 437,
    FOREGROUND = 438,
    FORMAT = 439,
    FORT14 = 440,
    FORT63 = 441,
    FORT64 = 442,
    FORWARD = 443,
    FRAMEP = 444,
    FREQ = 445,
    FRONTBUFFER = 446,
    GENERAL = 447,
    GETP = 448,
    GOTO = 449,
    GRAPH = 450,
    GRAPHNO = 451,
    GRAPHS = 452,
    GRAPHTYPE = 453,
    GRID = 454,
    HARDCOPY = 455,
    HBAR = 456,
    HELP = 457,
    HGAP = 458,
    HIDDEN = 459,
    HISTBOX = 460,
    HISTO = 461,
    HISTORY = 462,
    HMS = 463,
    HORIZONTAL = 464,
    HOURS = 465,
    HPGLL = 466,
    HPGLP = 467,
    IF = 468,
    IGNORE = 469,
    IHL = 470,
    IMAGE = 471,
    IMAGES = 472,
    IN = 473,
    INCLUDE = 474,
    INFO = 475,
    INIT = 476,
    INITGRAPHICS = 477,
    INOUT = 478,
    INTEGRATE = 479,
    INTERP = 480,
    INUNDATION = 481,
    INVDFT = 482,
    INVFFT = 483,
    ISOLINE = 484,
    ISOLINES = 485,
    JUST = 486,
    KILL = 487,
    KM = 488,
    LABEL = 489,
    LAYOUT = 490,
    LEAVE = 491,
    LEAVEGRAPHICS = 492,
    LEFT = 493,
    LEGEND = 494,
    LENGTH = 495,
    LEVEL = 496,
    LEVELS = 497,
    LIMITS = 498,
    LINE = 499,
    LINES = 500,
    LINESTYLE = 501,
    LINETO = 502,
    LINEW = 503,
    LINEWIDTH = 504,
    LINK = 505,
    LOAD = 506,
    LOC = 507,
    LOCATE = 508,
    LOCATOR = 509,
    LOCTYPE = 510,
    LOGX = 511,
    LOGXY = 512,
    LOGY = 513,
    M = 514,
    MAG = 515,
    MAGNITUDE = 516,
    MAJOR = 517,
    MAPSCALE = 518,
    MARKER = 519,
    MARKERS = 520,
    MAXLEVELS = 521,
    METHOD = 522,
    MIFL = 523,
    MIFP = 524,
    MILES = 525,
    MINOR = 526,
    MINUTES = 527,
    MISSINGP = 528,
    MM = 529,
    MMDD = 530,
    MMDDHMS = 531,
    MMDDYY = 532,
    MMDDYYHMS = 533,
    MMSSLAT = 534,
    MMSSLON = 535,
    MMYY = 536,
    MONTHDAY = 537,
    MONTHL = 538,
    MONTHS = 539,
    MOVE = 540,
    MOVE2 = 541,
    MOVETO = 542,
    NEGATE = 543,
    NO = 544,
    NODE = 545,
    NODES = 546,
    NONE = 547,
    NORMAL = 548,
    NORTH = 549,
    NXY = 550,
    OFF = 551,
    OFFSETX = 552,
    OFFSETY = 553,
    ON = 554,
    OP = 555,
    OPEN = 556,
    ORIENT = 557,
    OUT = 558,
    PAGE = 559,
    PARA = 560,
    PARALLEL = 561,
    PARAMETERS = 562,
    PARAMS = 563,
    PARMS = 564,
    PATTERN = 565,
    PER = 566,
    PERIMETER = 567,
    PERP = 568,
    PERPENDICULAR = 569,
    PHASE = 570,
    PIE = 571,
    PIPE = 572,
    PLACE = 573,
    PLAN = 574,
    PLUS = 575,
    POINT = 576,
    POLAR = 577,
    POLY = 578,
    POLYI = 579,
    POLYO = 580,
    POP = 581,
    POWER = 582,
    PREC = 583,
    PREFIX = 584,
    PREPEND = 585,
    PRINT = 586,
    PROFILE = 587,
    PROP = 588,
    PS = 589,
    PSCOLORL = 590,
    PSCOLORP = 591,
    PSMONOL = 592,
    PSMONOP = 593,
    PUSH = 594,
    PUTP = 595,
    QUIT = 596,
    READ = 597,
    READBIN = 598,
    REDRAW = 599,
    REGION = 600,
    REGIONS = 601,
    REGNUM = 602,
    REGRESS = 603,
    REMOVE = 604,
    RENDER = 605,
    REPORT = 606,
    RESET = 607,
    REVERSE = 608,
    REWIND = 609,
    RIGHT = 610,
    RISER = 611,
    ROT = 612,
    RUN = 613,
    SALINITY = 614,
    SAMPLE = 615,
    SAVE = 616,
    SCALAR = 617,
    SCALE = 618,
    SCIENTIFIC = 619,
    SECONDS = 620,
    SET = 621,
    SETS = 622,
    SHOW = 623,
    SHRINK = 624,
    SIGMA = 625,
    SIGN = 626,
    SIZE = 627,
    SKIP = 628,
    SLAB = 629,
    SLEEP = 630,
    SLICE = 631,
    SOURCE = 632,
    SPEC = 633,
    SPECIFIED = 634,
    SPECTRUM = 635,
    SPLITS = 636,
    SQUARE = 637,
    STACK = 638,
    STACKEDBAR = 639,
    STACKEDHBAR = 640,
    STACKEDLINE = 641,
    STAGGER = 642,
    STAR = 643,
    START = 644,
    STARTSTEP = 645,
    STARTTYPE = 646,
    STATION = 647,
    STATUS = 648,
    STEP = 649,
    STOP = 650,
    STREAMLINES = 651,
    STRING = 652,
    STRINGS = 653,
    SUBTITLE = 654,
    SURFACE = 655,
    SWAPBUFFER = 656,
    SYMBOL = 657,
    SYSTEM = 658,
    TEANL = 659,
    TEXT = 660,
    TICK = 661,
    TICKLABEL = 662,
    TICKMARKS = 663,
    TICKP = 664,
    TIDALCLOCK = 665,
    TIDESTATION = 666,
    TIME = 667,
    TIMEINFO = 668,
    TIMELINE = 669,
    TITLE = 670,
    TO = 671,
    TOP = 672,
    TOTAL = 673,
    TRACK = 674,
    TRANSECT = 675,
    TRIANGLE1 = 676,
    TRIANGLE2 = 677,
    TRIANGLE3 = 678,
    TRIANGLE4 = 679,
    TRUEP = 680,
    TYPE = 681,
    UNITS = 682,
    UP = 683,
    VALUE = 684,
    VECTOR = 685,
    VEL = 686,
    VELMARKER = 687,
    VELOCITY = 688,
    VERTICAL = 689,
    VGAP = 690,
    VIEW = 691,
    VSCALE = 692,
    VX1 = 693,
    VX2 = 694,
    VY1 = 695,
    VY2 = 696,
    WEEKS = 697,
    WET = 698,
    WETDRY = 699,
    WIDTH = 700,
    WIND = 701,
    WITH = 702,
    WORLD = 703,
    WRAP = 704,
    WRITE = 705,
    WSCALE = 706,
    WX1 = 707,
    WX2 = 708,
    WY1 = 709,
    WY2 = 710,
    X0 = 711,
    X1 = 712,
    X2 = 713,
    X3 = 714,
    X4 = 715,
    X5 = 716,
    XAXES = 717,
    XAXIS = 718,
    XCOR = 719,
    XMAX = 720,
    XMIN = 721,
    XY = 722,
    XYARC = 723,
    XYBOX = 724,
    XYDX = 725,
    XYDXDX = 726,
    XYDXDY = 727,
    XYDY = 728,
    XYDYDY = 729,
    XYFIXED = 730,
    XYHILO = 731,
    XYRT = 732,
    XYSEG = 733,
    XYSTRING = 734,
    XYUV = 735,
    XYX2Y2 = 736,
    XYXX = 737,
    XYYY = 738,
    XYZ = 739,
    XYZW = 740,
    Y0 = 741,
    Y1 = 742,
    Y2 = 743,
    Y3 = 744,
    Y4 = 745,
    Y5 = 746,
    YAXES = 747,
    YAXIS = 748,
    YEARS = 749,
    YES = 750,
    YMAX = 751,
    YMIN = 752,
    ZEROXAXIS = 753,
    ZEROYAXIS = 754,
    ZOOM = 755,
    ZOOMBOX = 756,
    OR = 757,
    AND = 758,
    GT = 759,
    LT = 760,
    LE = 761,
    GE = 762,
    EQ = 763,
    NE = 764,
    UMINUS = 765,
    NOT = 766
  };
#endif
/* Tokens.  */
#define VAR 258
#define X 259
#define Y 260
#define CHRSTR 261
#define FITPARM 262
#define NUMBER 263
#define ABS 264
#define ACOS 265
#define ASIN 266
#define ATAN 267
#define ATAN2 268
#define CEIL 269
#define COS 270
#define DEG 271
#define DX 272
#define DY 273
#define ERF 274
#define ERFC 275
#define EXP 276
#define FLOOR 277
#define HYPOT 278
#define INDEX 279
#define INT 280
#define IRAND 281
#define LGAMMA 282
#define LN 283
#define LOG 284
#define LOGISTIC 285
#define MAXP 286
#define MINP 287
#define MINMAX 288
#define MOD 289
#define NORM 290
#define NORMP 291
#define PI 292
#define RAD 293
#define RAND 294
#define SETNO 295
#define SIN 296
#define SQR 297
#define SQRT 298
#define TAN 299
#define INUM 300
#define ABORT 301
#define ABOVE 302
#define ABSOLUTE 303
#define ACTIVATE 304
#define ACTIVE 305
#define ADCIRC 306
#define ADCIRC3DFLOW 307
#define ALL 308
#define ALT 309
#define ALTERNATE 310
#define ALTXAXIS 311
#define ALTYAXIS 312
#define AMP 313
#define ANGLE 314
#define ANNOTATE 315
#define APPEND 316
#define AREA 317
#define ARROW 318
#define ASCEND 319
#define AT 320
#define ATTACH 321
#define AUTO 322
#define AUTOSCALE 323
#define AUTOTICKS 324
#define AVERAGE 325
#define AVG 326
#define AXES 327
#define AXIS 328
#define BACKBUFFER 329
#define BACKGROUND 330
#define BAR 331
#define BATCH 332
#define BATH 333
#define BATHYMETRY 334
#define COURANT 335
#define BELOW 336
#define BIN 337
#define BINARY 338
#define BOTH 339
#define BOTTOM 340
#define BOUNDARY 341
#define BOX 342
#define CELLS 343
#define CENTER 344
#define CH3D 345
#define CHAR 346
#define CHDIR 347
#define CIRCLE 348
#define CLEAR 349
#define CLICK 350
#define CLOCK 351
#define CLOSE 352
#define CM 353
#define CMAP 354
#define COLOR 355
#define COLORMAP 356
#define COMMENT 357
#define CONC 358
#define CONCENTRATION 359
#define CONCENTRATIONS 360
#define COPY 361
#define CROSS 362
#define CYCLE 363
#define DAYMONTH 364
#define DAYOFWEEKL 365
#define DAYOFWEEKS 366
#define DAYOFYEAR 367
#define DAYS 368
#define DDMMYY 369
#define DDMONTHSYYHHMMSS 370
#define DECIMAL 371
#define DEF 372
#define DEFAULT 373
#define DEGREESLAT 374
#define DEGREESLON 375
#define DEGREESMMLAT 376
#define DEGREESMMLON 377
#define DEGREESMMSSLAT 378
#define DEGREESMMSSLON 379
#define DELAYP 380
#define DELETE 381
#define DEPTH 382
#define DEPTHS 383
#define DESCEND 384
#define DEVICE 385
#define DEVXY 386
#define DFT 387
#define DT 388
#define DIAMOND 389
#define DIFFERENCE 390
#define DISK 391
#define DISPLAY 392
#define DOT 393
#define DOUBLEBUFFER 394
#define DOWN 395
#define DRAW2 396
#define DROGUE 397
#define DROGUES 398
#define DRY 399
#define DXDX 400
#define DXP 401
#define DYDY 402
#define DYP 403
#define ECHO 404
#define EDIT 405
#define ELA 406
#define ELCIRC 407
#define ELEMENT 408
#define ELEMENTS 409
#define ELEV 410
#define ELEVATION 411
#define ELEVATIONS 412
#define ELEVMARKER 413
#define ELLIPSE 414
#define ELLIPSES 415
#define ELLIPSEZ 416
#define ELSE 417
#define END 418
#define ERRORBAR 419
#define EXIT 420
#define EXPAND 421
#define EXPONENTIAL 422
#define FACTOR 423
#define FALSEP 424
#define FAST 425
#define FEET 426
#define FFT 427
#define FILEP 428
#define FILL 429
#define FIND 430
#define FIXEDPOINT 431
#define FLOW 432
#define FLUSH 433
#define FLUX 434
#define FOCUS 435
#define FOLLOWS 436
#define FONTP 437
#define FOREGROUND 438
#define FORMAT 439
#define FORT14 440
#define FORT63 441
#define FORT64 442
#define FORWARD 443
#define FRAMEP 444
#define FREQ 445
#define FRONTBUFFER 446
#define GENERAL 447
#define GETP 448
#define GOTO 449
#define GRAPH 450
#define GRAPHNO 451
#define GRAPHS 452
#define GRAPHTYPE 453
#define GRID 454
#define HARDCOPY 455
#define HBAR 456
#define HELP 457
#define HGAP 458
#define HIDDEN 459
#define HISTBOX 460
#define HISTO 461
#define HISTORY 462
#define HMS 463
#define HORIZONTAL 464
#define HOURS 465
#define HPGLL 466
#define HPGLP 467
#define IF 468
#define IGNORE 469
#define IHL 470
#define IMAGE 471
#define IMAGES 472
#define IN 473
#define INCLUDE 474
#define INFO 475
#define INIT 476
#define INITGRAPHICS 477
#define INOUT 478
#define INTEGRATE 479
#define INTERP 480
#define INUNDATION 481
#define INVDFT 482
#define INVFFT 483
#define ISOLINE 484
#define ISOLINES 485
#define JUST 486
#define KILL 487
#define KM 488
#define LABEL 489
#define LAYOUT 490
#define LEAVE 491
#define LEAVEGRAPHICS 492
#define LEFT 493
#define LEGEND 494
#define LENGTH 495
#define LEVEL 496
#define LEVELS 497
#define LIMITS 498
#define LINE 499
#define LINES 500
#define LINESTYLE 501
#define LINETO 502
#define LINEW 503
#define LINEWIDTH 504
#define LINK 505
#define LOAD 506
#define LOC 507
#define LOCATE 508
#define LOCATOR 509
#define LOCTYPE 510
#define LOGX 511
#define LOGXY 512
#define LOGY 513
#define M 514
#define MAG 515
#define MAGNITUDE 516
#define MAJOR 517
#define MAPSCALE 518
#define MARKER 519
#define MARKERS 520
#define MAXLEVELS 521
#define METHOD 522
#define MIFL 523
#define MIFP 524
#define MILES 525
#define MINOR 526
#define MINUTES 527
#define MISSINGP 528
#define MM 529
#define MMDD 530
#define MMDDHMS 531
#define MMDDYY 532
#define MMDDYYHMS 533
#define MMSSLAT 534
#define MMSSLON 535
#define MMYY 536
#define MONTHDAY 537
#define MONTHL 538
#define MONTHS 539
#define MOVE 540
#define MOVE2 541
#define MOVETO 542
#define NEGATE 543
#define NO 544
#define NODE 545
#define NODES 546
#define NONE 547
#define NORMAL 548
#define NORTH 549
#define NXY 550
#define OFF 551
#define OFFSETX 552
#define OFFSETY 553
#define ON 554
#define OP 555
#define OPEN 556
#define ORIENT 557
#define OUT 558
#define PAGE 559
#define PARA 560
#define PARALLEL 561
#define PARAMETERS 562
#define PARAMS 563
#define PARMS 564
#define PATTERN 565
#define PER 566
#define PERIMETER 567
#define PERP 568
#define PERPENDICULAR 569
#define PHASE 570
#define PIE 571
#define PIPE 572
#define PLACE 573
#define PLAN 574
#define PLUS 575
#define POINT 576
#define POLAR 577
#define POLY 578
#define POLYI 579
#define POLYO 580
#define POP 581
#define POWER 582
#define PREC 583
#define PREFIX 584
#define PREPEND 585
#define PRINT 586
#define PROFILE 587
#define PROP 588
#define PS 589
#define PSCOLORL 590
#define PSCOLORP 591
#define PSMONOL 592
#define PSMONOP 593
#define PUSH 594
#define PUTP 595
#define QUIT 596
#define READ 597
#define READBIN 598
#define REDRAW 599
#define REGION 600
#define REGIONS 601
#define REGNUM 602
#define REGRESS 603
#define REMOVE 604
#define RENDER 605
#define REPORT 606
#define RESET 607
#define REVERSE 608
#define REWIND 609
#define RIGHT 610
#define RISER 611
#define ROT 612
#define RUN 613
#define SALINITY 614
#define SAMPLE 615
#define SAVE 616
#define SCALAR 617
#define SCALE 618
#define SCIENTIFIC 619
#define SECONDS 620
#define SET 621
#define SETS 622
#define SHOW 623
#define SHRINK 624
#define SIGMA 625
#define SIGN 626
#define SIZE 627
#define SKIP 628
#define SLAB 629
#define SLEEP 630
#define SLICE 631
#define SOURCE 632
#define SPEC 633
#define SPECIFIED 634
#define SPECTRUM 635
#define SPLITS 636
#define SQUARE 637
#define STACK 638
#define STACKEDBAR 639
#define STACKEDHBAR 640
#define STACKEDLINE 641
#define STAGGER 642
#define STAR 643
#define START 644
#define STARTSTEP 645
#define STARTTYPE 646
#define STATION 647
#define STATUS 648
#define STEP 649
#define STOP 650
#define STREAMLINES 651
#define STRING 652
#define STRINGS 653
#define SUBTITLE 654
#define SURFACE 655
#define SWAPBUFFER 656
#define SYMBOL 657
#define SYSTEM 658
#define TEANL 659
#define TEXT 660
#define TICK 661
#define TICKLABEL 662
#define TICKMARKS 663
#define TICKP 664
#define TIDALCLOCK 665
#define TIDESTATION 666
#define TIME 667
#define TIMEINFO 668
#define TIMELINE 669
#define TITLE 670
#define TO 671
#define TOP 672
#define TOTAL 673
#define TRACK 674
#define TRANSECT 675
#define TRIANGLE1 676
#define TRIANGLE2 677
#define TRIANGLE3 678
#define TRIANGLE4 679
#define TRUEP 680
#define TYPE 681
#define UNITS 682
#define UP 683
#define VALUE 684
#define VECTOR 685
#define VEL 686
#define VELMARKER 687
#define VELOCITY 688
#define VERTICAL 689
#define VGAP 690
#define VIEW 691
#define VSCALE 692
#define VX1 693
#define VX2 694
#define VY1 695
#define VY2 696
#define WEEKS 697
#define WET 698
#define WETDRY 699
#define WIDTH 700
#define WIND 701
#define WITH 702
#define WORLD 703
#define WRAP 704
#define WRITE 705
#define WSCALE 706
#define WX1 707
#define WX2 708
#define WY1 709
#define WY2 710
#define X0 711
#define X1 712
#define X2 713
#define X3 714
#define X4 715
#define X5 716
#define XAXES 717
#define XAXIS 718
#define XCOR 719
#define XMAX 720
#define XMIN 721
#define XY 722
#define XYARC 723
#define XYBOX 724
#define XYDX 725
#define XYDXDX 726
#define XYDXDY 727
#define XYDY 728
#define XYDYDY 729
#define XYFIXED 730
#define XYHILO 731
#define XYRT 732
#define XYSEG 733
#define XYSTRING 734
#define XYUV 735
#define XYX2Y2 736
#define XYXX 737
#define XYYY 738
#define XYZ 739
#define XYZW 740
#define Y0 741
#define Y1 742
#define Y2 743
#define Y3 744
#define Y4 745
#define Y5 746
#define YAXES 747
#define YAXIS 748
#define YEARS 749
#define YES 750
#define YMAX 751
#define YMIN 752
#define ZEROXAXIS 753
#define ZEROYAXIS 754
#define ZOOM 755
#define ZOOMBOX 756
#define OR 757
#define AND 758
#define GT 759
#define LT 760
#define LE 761
#define GE 762
#define EQ 763
#define NE 764
#define UMINUS 765
#define NOT 766

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE YYSTYPE;
union YYSTYPE
{
#line 113 "gram.y" /* yacc.c:355  */

    double val;
    int ival;
    double *ptr;
    int func;
    int pset;
    char *str;

#line 1239 "y.tab.c" /* yacc.c:355  */
};
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_Y_TAB_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 1254 "y.tab.c" /* yacc.c:358  */

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

#if !defined _Noreturn \
     && (!defined __STDC_VERSION__ || __STDC_VERSION__ < 201112)
# if defined _MSC_VER && 1200 <= _MSC_VER
#  define _Noreturn __declspec (noreturn)
# else
#  define _Noreturn YY_ATTRIBUTE ((__noreturn__))
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif


#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  522
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   6977

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  526
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  72
/* YYNRULES -- Number of rules.  */
#define YYNRULES  797
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  1798

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   766

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint16 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     519,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,   515,     2,     2,
     523,   524,   513,   511,   520,   512,   525,   514,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   502,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   521,     2,   522,   516,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   110,   111,   112,   113,   114,
     115,   116,   117,   118,   119,   120,   121,   122,   123,   124,
     125,   126,   127,   128,   129,   130,   131,   132,   133,   134,
     135,   136,   137,   138,   139,   140,   141,   142,   143,   144,
     145,   146,   147,   148,   149,   150,   151,   152,   153,   154,
     155,   156,   157,   158,   159,   160,   161,   162,   163,   164,
     165,   166,   167,   168,   169,   170,   171,   172,   173,   174,
     175,   176,   177,   178,   179,   180,   181,   182,   183,   184,
     185,   186,   187,   188,   189,   190,   191,   192,   193,   194,
     195,   196,   197,   198,   199,   200,   201,   202,   203,   204,
     205,   206,   207,   208,   209,   210,   211,   212,   213,   214,
     215,   216,   217,   218,   219,   220,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
     235,   236,   237,   238,   239,   240,   241,   242,   243,   244,
     245,   246,   247,   248,   249,   250,   251,   252,   253,   254,
     255,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   347,   348,   349,   350,   351,   352,   353,   354,
     355,   356,   357,   358,   359,   360,   361,   362,   363,   364,
     365,   366,   367,   368,   369,   370,   371,   372,   373,   374,
     375,   376,   377,   378,   379,   380,   381,   382,   383,   384,
     385,   386,   387,   388,   389,   390,   391,   392,   393,   394,
     395,   396,   397,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,   412,   413,   414,
     415,   416,   417,   418,   419,   420,   421,   422,   423,   424,
     425,   426,   427,   428,   429,   430,   431,   432,   433,   434,
     435,   436,   437,   438,   439,   440,   441,   442,   443,   444,
     445,   446,   447,   448,   449,   450,   451,   452,   453,   454,
     455,   456,   457,   458,   459,   460,   461,   462,   463,   464,
     465,   466,   467,   468,   469,   470,   471,   472,   473,   474,
     475,   476,   477,   478,   479,   480,   481,   482,   483,   484,
     485,   486,   487,   488,   489,   490,   491,   492,   493,   494,
     495,   496,   497,   498,   499,   500,   501,   503,   504,   505,
     506,   507,   508,   509,   510,   517,   518
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   668,   668,   669,   670,   671,   672,   673,   674,   675,
     676,   677,   678,   679,   680,   681,   682,   683,   684,   685,
     686,   687,   688,   689,   690,   691,   692,   695,   696,   697,
     698,   699,   700,   701,   710,   711,   712,   713,   714,   715,
     716,   717,   729,   730,   731,   732,   733,   734,   741,   742,
     743,   744,   745,   746,   747,   748,   760,   761,   762,   763,
     764,   765,   770,   771,   772,   773,   774,   775,   776,   777,
     778,   779,   780,   781,   782,   802,   803,   804,   805,   806,
     807,   808,   809,   810,   811,   812,   813,   817,   818,   824,
     830,   831,   832,   833,   834,   835,   836,   837,   838,   839,
     840,   841,   842,   843,   844,   854,   863,   864,   865,   869,
     871,   873,   875,   875,   876,   876,   880,   883,   884,   885,
     889,   890,   899,   904,   905,   906,   910,   914,   923,   932,
     941,   951,   960,   969,   978,   987,   991,   992,   992,   993,
     994,   995,   996,   997,   998,   999,  1000,  1001,  1002,  1003,
    1004,  1010,  1011,  1012,  1013,  1014,  1019,  1028,  1036,  1044,
    1052,  1060,  1068,  1076,  1082,  1086,  1090,  1094,  1098,  1102,
    1106,  1113,  1120,  1128,  1136,  1136,  1137,  1137,  1138,  1138,
    1139,  1139,  1140,  1145,  1146,  1147,  1148,  1149,  1150,  1151,
    1152,  1159,  1164,  1169,  1177,  1178,  1179,  1180,  1181,  1182,
    1183,  1184,  1185,  1186,  1187,  1188,  1189,  1189,  1190,  1190,
    1191,  1191,  1192,  1192,  1193,  1193,  1194,  1195,  1196,  1197,
    1198,  1199,  1200,  1201,  1202,  1203,  1204,  1205,  1206,  1207,
    1208,  1212,  1213,  1214,  1215,  1216,  1216,  1217,  1218,  1219,
    1224,  1236,  1237,  1238,  1239,  1240,  1241,  1242,  1243,  1243,
    1244,  1244,  1245,  1246,  1247,  1248,  1249,  1249,  1250,  1255,
    1260,  1269,  1270,  1271,  1271,  1272,  1283,  1294,  1295,  1296,
    1297,  1298,  1299,  1300,  1301,  1302,  1303,  1311,  1319,  1327,
    1336,  1339,  1349,  1350,  1351,  1352,  1356,  1357,  1358,  1359,
    1359,  1360,  1360,  1361,  1365,  1372,  1373,  1373,  1374,  1379,
    1380,  1381,  1382,  1383,  1384,  1385,  1386,  1387,  1388,  1389,
    1390,  1391,  1392,  1399,  1404,  1409,  1417,  1418,  1418,  1419,
    1424,  1425,  1426,  1427,  1428,  1429,  1430,  1431,  1432,  1433,
    1434,  1441,  1446,  1451,  1459,  1460,  1461,  1462,  1463,  1464,
    1469,  1480,  1485,  1486,  1487,  1488,  1489,  1490,  1495,  1506,
    1510,  1511,  1512,  1513,  1514,  1515,  1520,  1531,  1535,  1536,
    1537,  1538,  1539,  1540,  1545,  1549,  1550,  1555,  1560,  1561,
    1562,  1563,  1564,  1565,  1566,  1567,  1571,  1572,  1573,  1574,
    1575,  1576,  1577,  1578,  1579,  1580,  1581,  1582,  1587,  1591,
    1592,  1592,  1593,  1598,  1599,  1600,  1601,  1602,  1603,  1604,
    1605,  1606,  1613,  1618,  1623,  1631,  1632,  1633,  1634,  1635,
    1636,  1637,  1638,  1639,  1640,  1641,  1642,  1643,  1644,  1645,
    1646,  1650,  1651,  1652,  1653,  1654,  1659,  1669,  1670,  1671,
    1672,  1673,  1674,  1675,  1676,  1677,  1681,  1685,  1686,  1693,
    1694,  1695,  1696,  1697,  1704,  1705,  1706,  1707,  1708,  1711,
    1714,  1715,  1718,  1722,  1725,  1728,  1729,  1733,  1737,  1738,
    1739,  1740,  1741,  1742,  1743,  1744,  1745,  1746,  1747,  1748,
    1749,  1750,  1751,  1752,  1757,  1762,  1767,  1771,  1772,  1773,
    1774,  1775,  1779,  1780,  1781,  1782,  1783,  1784,  1788,  1789,
    1790,  1794,  1795,  1796,  1797,  1798,  1799,  1800,  1805,  1806,
    1807,  1808,  1809,  1818,  1823,  1824,  1825,  1826,  1827,  1828,
    1829,  1830,  1831,  1832,  1833,  1834,  1835,  1836,  1837,  1838,
    1839,  1840,  1841,  1842,  1843,  1844,  1845,  1846,  1847,  1848,
    1849,  1850,  1851,  1852,  1856,  1857,  1861,  1862,  1863,  1864,
    1865,  1866,  1867,  1868,  1869,  1870,  1871,  1872,  1873,  1874,
    1875,  1876,  1877,  1878,  1879,  1880,  1881,  1882,  1883,  1884,
    1885,  1886,  1887,  1888,  1889,  1890,  1894,  1895,  1896,  1897,
    1898,  1899,  1903,  1904,  1905,  1906,  1907,  1911,  1912,  1913,
    1914,  1926,  1927,  1931,  1932,  1933,  1942,  1943,  1944,  1945,
    1946,  1950,  1951,  1952,  1956,  1957,  1958,  1971,  1972,  1973,
    1974,  1975,  1976,  1977,  1978,  1979,  1980,  1981,  1982,  1983,
    1984,  1985,  1986,  1987,  1988,  1989,  1990,  1991,  1992,  1993,
    1994,  1995,  1996,  1997,  2000,  2001,  2005,  2006,  2007,  2008,
    2012,  2013,  2016,  2017,  2028,  2029,  2032,  2033,  2034,  2037,
    2038,  2039,  2040,  2041,  2042,  2043,  2047,  2048,  2049,  2050,
    2064,  2078,  2086,  2097,  2106,  2115,  2124,  2133,  2142,  2151,
    2160,  2169,  2178,  2187,  2196,  2205,  2214,  2223,  2236,  2251,
    2266,  2279,  2288,  2297,  2306,  2315,  2324,  2333,  2342,  2351,
    2360,  2369,  2378,  2387,  2396,  2405,  2414,  2423,  2432,  2441,
    2450,  2459,  2468,  2477,  2486,  2495,  2504,  2513,  2522,  2531,
    2540,  2549,  2558,  2567,  2576,  2585,  2594,  2604,  2613,  2622,
    2631,  2640,  2649,  2658,  2667,  2676,  2685,  2694,  2703,  2712,
    2721,  2730,  2739,  2748,  2757,  2767,  2768,  2771,  2774,  2777,
    2780,  2783,  2792,  2795,  2798,  2801,  2804,  2807,  2810,  2813,
    2816,  2819,  2822,  2825,  2828,  2831,  2834,  2837,  2840,  2843,
    2846,  2849,  2852,  2855,  2858,  2861,  2864,  2867,  2870,  2873,
    2876,  2879,  2882,  2885,  2888,  2891,  2894,  2897,  2900,  2903,
    2906,  2909,  2912,  2916,  2919,  2922,  2925,  2928,  2932,  2935,
    2938,  2941,  2944,  2947,  2950,  2953,  2957,  2964,  2967,  2970,
    2973,  2976,  2979,  2982,  2985,  2988,  2991,  2994
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "VAR", "X", "Y", "CHRSTR", "FITPARM",
  "NUMBER", "ABS", "ACOS", "ASIN", "ATAN", "ATAN2", "CEIL", "COS", "DEG",
  "DX", "DY", "ERF", "ERFC", "EXP", "FLOOR", "HYPOT", "INDEX", "INT",
  "IRAND", "LGAMMA", "LN", "LOG", "LOGISTIC", "MAXP", "MINP", "MINMAX",
  "MOD", "NORM", "NORMP", "PI", "RAD", "RAND", "SETNO", "SIN", "SQR",
  "SQRT", "TAN", "INUM", "ABORT", "ABOVE", "ABSOLUTE", "ACTIVATE",
  "ACTIVE", "ADCIRC", "ADCIRC3DFLOW", "ALL", "ALT", "ALTERNATE",
  "ALTXAXIS", "ALTYAXIS", "AMP", "ANGLE", "ANNOTATE", "APPEND", "AREA",
  "ARROW", "ASCEND", "AT", "ATTACH", "AUTO", "AUTOSCALE", "AUTOTICKS",
  "AVERAGE", "AVG", "AXES", "AXIS", "BACKBUFFER", "BACKGROUND", "BAR",
  "BATCH", "BATH", "BATHYMETRY", "COURANT", "BELOW", "BIN", "BINARY",
  "BOTH", "BOTTOM", "BOUNDARY", "BOX", "CELLS", "CENTER", "CH3D", "CHAR",
  "CHDIR", "CIRCLE", "CLEAR", "CLICK", "CLOCK", "CLOSE", "CM", "CMAP",
  "COLOR", "COLORMAP", "COMMENT", "CONC", "CONCENTRATION",
  "CONCENTRATIONS", "COPY", "CROSS", "CYCLE", "DAYMONTH", "DAYOFWEEKL",
  "DAYOFWEEKS", "DAYOFYEAR", "DAYS", "DDMMYY", "DDMONTHSYYHHMMSS",
  "DECIMAL", "DEF", "DEFAULT", "DEGREESLAT", "DEGREESLON", "DEGREESMMLAT",
  "DEGREESMMLON", "DEGREESMMSSLAT", "DEGREESMMSSLON", "DELAYP", "DELETE",
  "DEPTH", "DEPTHS", "DESCEND", "DEVICE", "DEVXY", "DFT", "DT", "DIAMOND",
  "DIFFERENCE", "DISK", "DISPLAY", "DOT", "DOUBLEBUFFER", "DOWN", "DRAW2",
  "DROGUE", "DROGUES", "DRY", "DXDX", "DXP", "DYDY", "DYP", "ECHO", "EDIT",
  "ELA", "ELCIRC", "ELEMENT", "ELEMENTS", "ELEV", "ELEVATION",
  "ELEVATIONS", "ELEVMARKER", "ELLIPSE", "ELLIPSES", "ELLIPSEZ", "ELSE",
  "END", "ERRORBAR", "EXIT", "EXPAND", "EXPONENTIAL", "FACTOR", "FALSEP",
  "FAST", "FEET", "FFT", "FILEP", "FILL", "FIND", "FIXEDPOINT", "FLOW",
  "FLUSH", "FLUX", "FOCUS", "FOLLOWS", "FONTP", "FOREGROUND", "FORMAT",
  "FORT14", "FORT63", "FORT64", "FORWARD", "FRAMEP", "FREQ", "FRONTBUFFER",
  "GENERAL", "GETP", "GOTO", "GRAPH", "GRAPHNO", "GRAPHS", "GRAPHTYPE",
  "GRID", "HARDCOPY", "HBAR", "HELP", "HGAP", "HIDDEN", "HISTBOX", "HISTO",
  "HISTORY", "HMS", "HORIZONTAL", "HOURS", "HPGLL", "HPGLP", "IF",
  "IGNORE", "IHL", "IMAGE", "IMAGES", "IN", "INCLUDE", "INFO", "INIT",
  "INITGRAPHICS", "INOUT", "INTEGRATE", "INTERP", "INUNDATION", "INVDFT",
  "INVFFT", "ISOLINE", "ISOLINES", "JUST", "KILL", "KM", "LABEL", "LAYOUT",
  "LEAVE", "LEAVEGRAPHICS", "LEFT", "LEGEND", "LENGTH", "LEVEL", "LEVELS",
  "LIMITS", "LINE", "LINES", "LINESTYLE", "LINETO", "LINEW", "LINEWIDTH",
  "LINK", "LOAD", "LOC", "LOCATE", "LOCATOR", "LOCTYPE", "LOGX", "LOGXY",
  "LOGY", "M", "MAG", "MAGNITUDE", "MAJOR", "MAPSCALE", "MARKER",
  "MARKERS", "MAXLEVELS", "METHOD", "MIFL", "MIFP", "MILES", "MINOR",
  "MINUTES", "MISSINGP", "MM", "MMDD", "MMDDHMS", "MMDDYY", "MMDDYYHMS",
  "MMSSLAT", "MMSSLON", "MMYY", "MONTHDAY", "MONTHL", "MONTHS", "MOVE",
  "MOVE2", "MOVETO", "NEGATE", "NO", "NODE", "NODES", "NONE", "NORMAL",
  "NORTH", "NXY", "OFF", "OFFSETX", "OFFSETY", "ON", "OP", "OPEN",
  "ORIENT", "OUT", "PAGE", "PARA", "PARALLEL", "PARAMETERS", "PARAMS",
  "PARMS", "PATTERN", "PER", "PERIMETER", "PERP", "PERPENDICULAR", "PHASE",
  "PIE", "PIPE", "PLACE", "PLAN", "PLUS", "POINT", "POLAR", "POLY",
  "POLYI", "POLYO", "POP", "POWER", "PREC", "PREFIX", "PREPEND", "PRINT",
  "PROFILE", "PROP", "PS", "PSCOLORL", "PSCOLORP", "PSMONOL", "PSMONOP",
  "PUSH", "PUTP", "QUIT", "READ", "READBIN", "REDRAW", "REGION", "REGIONS",
  "REGNUM", "REGRESS", "REMOVE", "RENDER", "REPORT", "RESET", "REVERSE",
  "REWIND", "RIGHT", "RISER", "ROT", "RUN", "SALINITY", "SAMPLE", "SAVE",
  "SCALAR", "SCALE", "SCIENTIFIC", "SECONDS", "SET", "SETS", "SHOW",
  "SHRINK", "SIGMA", "SIGN", "SIZE", "SKIP", "SLAB", "SLEEP", "SLICE",
  "SOURCE", "SPEC", "SPECIFIED", "SPECTRUM", "SPLITS", "SQUARE", "STACK",
  "STACKEDBAR", "STACKEDHBAR", "STACKEDLINE", "STAGGER", "STAR", "START",
  "STARTSTEP", "STARTTYPE", "STATION", "STATUS", "STEP", "STOP",
  "STREAMLINES", "STRING", "STRINGS", "SUBTITLE", "SURFACE", "SWAPBUFFER",
  "SYMBOL", "SYSTEM", "TEANL", "TEXT", "TICK", "TICKLABEL", "TICKMARKS",
  "TICKP", "TIDALCLOCK", "TIDESTATION", "TIME", "TIMEINFO", "TIMELINE",
  "TITLE", "TO", "TOP", "TOTAL", "TRACK", "TRANSECT", "TRIANGLE1",
  "TRIANGLE2", "TRIANGLE3", "TRIANGLE4", "TRUEP", "TYPE", "UNITS", "UP",
  "VALUE", "VECTOR", "VEL", "VELMARKER", "VELOCITY", "VERTICAL", "VGAP",
  "VIEW", "VSCALE", "VX1", "VX2", "VY1", "VY2", "WEEKS", "WET", "WETDRY",
  "WIDTH", "WIND", "WITH", "WORLD", "WRAP", "WRITE", "WSCALE", "WX1",
  "WX2", "WY1", "WY2", "X0", "X1", "X2", "X3", "X4", "X5", "XAXES",
  "XAXIS", "XCOR", "XMAX", "XMIN", "XY", "XYARC", "XYBOX", "XYDX",
  "XYDXDX", "XYDXDY", "XYDY", "XYDYDY", "XYFIXED", "XYHILO", "XYRT",
  "XYSEG", "XYSTRING", "XYUV", "XYX2Y2", "XYXX", "XYYY", "XYZ", "XYZW",
  "Y0", "Y1", "Y2", "Y3", "Y4", "Y5", "YAXES", "YAXIS", "YEARS", "YES",
  "YMAX", "YMIN", "ZEROXAXIS", "ZEROYAXIS", "ZOOM", "ZOOMBOX", "'='", "OR",
  "AND", "GT", "LT", "LE", "GE", "EQ", "NE", "'+'", "'-'", "'*'", "'/'",
  "'%'", "'^'", "UMINUS", "NOT", "'\\n'", "','", "'['", "']'", "'('",
  "')'", "'.'", "$accept", "list", "annotation", "animation", "misc",
  "$@1", "$@2", "models", "$@3", "$@4", "$@5", "$@6", "$@7", "flowprops",
  "$@8", "$@9", "$@10", "$@11", "$@12", "elevmarker", "$@13", "grid",
  "$@14", "$@15", "$@16", "drogues", "$@17", "isolines", "$@18", "$@19",
  "histboxes", "$@20", "zoomboxes", "$@21", "wscale", "vscale", "mapscale",
  "tidalclock", "timeinfo", "timeline", "slice", "$@22", "props", "graph",
  "setaxis", "axis", "allaxes", "axesprops", "axisfeature", "tickattr",
  "ticklabeldesc", "ticklabelattr", "axislabeldesc", "axisbardesc",
  "sourcetype", "justchoice", "graphtype", "inoutchoice", "signchoice",
  "formatchoice", "horv", "flowonoff", "onoff", "worldview", "torf",
  "filltype", "opchoice", "units", "asgn", "vasgn", "vexpr", "expr", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   347,   348,   349,   350,   351,   352,   353,   354,
     355,   356,   357,   358,   359,   360,   361,   362,   363,   364,
     365,   366,   367,   368,   369,   370,   371,   372,   373,   374,
     375,   376,   377,   378,   379,   380,   381,   382,   383,   384,
     385,   386,   387,   388,   389,   390,   391,   392,   393,   394,
     395,   396,   397,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,   412,   413,   414,
     415,   416,   417,   418,   419,   420,   421,   422,   423,   424,
     425,   426,   427,   428,   429,   430,   431,   432,   433,   434,
     435,   436,   437,   438,   439,   440,   441,   442,   443,   444,
     445,   446,   447,   448,   449,   450,   451,   452,   453,   454,
     455,   456,   457,   458,   459,   460,   461,   462,   463,   464,
     465,   466,   467,   468,   469,   470,   471,   472,   473,   474,
     475,   476,   477,   478,   479,   480,   481,   482,   483,   484,
     485,   486,   487,   488,   489,   490,   491,   492,   493,   494,
     495,   496,   497,   498,   499,   500,   501,   502,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,   517,   518,   519,   520,   521,   522,   523,   524,
     525,   526,   527,   528,   529,   530,   531,   532,   533,   534,
     535,   536,   537,   538,   539,   540,   541,   542,   543,   544,
     545,   546,   547,   548,   549,   550,   551,   552,   553,   554,
     555,   556,   557,   558,   559,   560,   561,   562,   563,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     575,   576,   577,   578,   579,   580,   581,   582,   583,   584,
     585,   586,   587,   588,   589,   590,   591,   592,   593,   594,
     595,   596,   597,   598,   599,   600,   601,   602,   603,   604,
     605,   606,   607,   608,   609,   610,   611,   612,   613,   614,
     615,   616,   617,   618,   619,   620,   621,   622,   623,   624,
     625,   626,   627,   628,   629,   630,   631,   632,   633,   634,
     635,   636,   637,   638,   639,   640,   641,   642,   643,   644,
     645,   646,   647,   648,   649,   650,   651,   652,   653,   654,
     655,   656,   657,   658,   659,   660,   661,   662,   663,   664,
     665,   666,   667,   668,   669,   670,   671,   672,   673,   674,
     675,   676,   677,   678,   679,   680,   681,   682,   683,   684,
     685,   686,   687,   688,   689,   690,   691,   692,   693,   694,
     695,   696,   697,   698,   699,   700,   701,   702,   703,   704,
     705,   706,   707,   708,   709,   710,   711,   712,   713,   714,
     715,   716,   717,   718,   719,   720,   721,   722,   723,   724,
     725,   726,   727,   728,   729,   730,   731,   732,   733,   734,
     735,   736,   737,   738,   739,   740,   741,   742,   743,   744,
     745,   746,   747,   748,   749,   750,   751,   752,   753,   754,
     755,   756,    61,   757,   758,   759,   760,   761,   762,   763,
     764,    43,    45,    42,    47,    37,    94,   765,   766,    10,
      44,    91,    93,    40,    41,    46
};
# endif

#define YYPACT_NINF -954

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-954)))

#define YYTABLE_NINF -798

#define yytable_value_is_error(Yytable_value) \
  (!!((Yytable_value) == (-798)))

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
    1797,  -452,  -423,  -954,  -954,  -448,  -443,  -431,  -427,  -415,
    -407,  -401,  -402,  -379,  -375,  -372,  -368,  -319,  -312,  -310,
    -371,  -308,  -304,  -287,  -261,    14,  -254,  -241,  -235,  -226,
    -219,  -190,  -225,  -208,  -191,  -179,  -181,  -180,  -167,  -165,
     687,  -954,  -954,  -954,   676,  -269,  2451,    94,   -34,   318,
     354,   360,   364,  -954,   239,   372,   728,    60,  -954,   285,
    -954,   177,   378,   -29,   -51,   741,   -43,  -125,   187,  3898,
      18,   114,  2904,   111,  -195,  -954,    32,   504,  -954,  -954,
    -954,  -954,    -4,  -954,  -954,  -954,  -954,    53,  -954,   266,
    -954,   420,   958,  -109,    -5,  -954,  -954,  2314,    11,   425,
     687,   165,  2357,   996,    34,  2068,   403,  -954,  -954,  -954,
    -954,   558,  3586,  -195,   485,  -954,  -954,  -954,  -954,   676,
    -954,   676,  -954,  -954,  -954,  4192,   655,  4234,  4234,  4234,
     433,   -85,   -80,   -74,   -73,   -70,   -68,   -67,   -65,   -64,
     -63,   -56,   -52,   -44,   -42,   -41,   -39,   -28,   -25,   -22,
     -40,  -954,   -18,   -14,  6430,  6413,  -954,  4234,  4192,  4234,
    4234,  4234,  4234,  4234,  4234,  4234,  4234,  4234,  4234,  4234,
    4234,  4234,   464,  4234,  4234,  -954,  4234,  4234,  4234,  4234,
    4234,  4234,  4234,  4234,  4234,  4234,  4234,  -954,   466,   119,
     571,   456,   487,  -161,   491,  -954,   323,  -954,   181,   547,
     548,  -195,   555,   570,  -954,  -954,  -954,  -954,   573,   572,
      65,    59,    64,    82,    83,    85,    87,   102,  -954,  -954,
    -954,   103,   104,   106,   107,   109,  -954,   116,   117,   118,
     123,   125,   126,   127,   130,   131,   134,   136,  -954,  -954,
    -954,  -954,   138,   153,   167,   168,   629,  -954,    54,   169,
     644,   684,  -358,  4192,  4192,  4192,  -954,  4955,  -954,  -954,
    -954,  -954,  -954,  -954,   173,   175,   176,  -195,   365,  -954,
     689,  1007,   701,   716,  -954,  4192,  -195,  -195,  4192,  -358,
     703,  -195,   365,  -954,  -954,  -954,  -954,  -954,   612,   706,
    -195,   707,   710,   711,  -954,  -954,  -134,   463,  -151,  -195,
    -137,   736,  -954,   -40,  -954,  -954,   586,   591,   595,   596,
    -195,   597,  -195,   598,   365,   730,   731,   537,   646,   739,
    4192,  -358,   743,   744,   753,  4192,  4192,  4192,   365,  4192,
     758,   -46,   169,  4973,  -954,  4192,   339,  3015,  -143,   342,
    4192,   764,   365,  -954,  -954,    -1,   767,  -954,   169,   769,
     770,  -358,  -954,  4991,  -195,  -954,    84,   771,   772,  4192,
    -358,   -32,   774,   101,  -954,  -954,  -954,   775,  -954,  -954,
    -954,  6447,   630,   782,    30,   785,  -954,   270,   789,    84,
     376,  -954,   796,   797,  -173,   708,   798,  4192,  -358,   799,
    4192,  4192,  4192,   365,  -954,  -954,  -954,  -954,  -954,  4192,
     437,   807,   810,   812,   169,   817,   823,  -358,   826,    24,
    -954,  5009,  -954,   827,   829,   830,   831,  -954,  -954,   832,
     746,  4192,  -358,   429,  -954,  -954,   472,   842,   843,   846,
     847,   848,  -358,   849,   854,  -954,  5027,   853,   763,   858,
    4192,  -358,   861,   864,   865,   866,   868,   872,  -954,  -954,
    -954,   879,   880,   882,   883,   884,   889,   892,   893,  5045,
     894,   895,  4192,  -358,   898,   101,  -954,  -954,   900,   901,
     905,     2,   908,   909,  -954,   911,   912,   913,   915,   916,
     920,   921,  4192,  4192,  4192,  4192,  5063,  -954,   923,   928,
    4192,  -358,   929,   101,  -954,  -954,  -954,  -954,  5081,   930,
     931,  -158,   840,   933,  4192,  -358,   938,   939,  4192,  4192,
    4192,   940,   365,  -954,    65,   426,  -954,   435,  -954,   446,
    2091,  1085,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,   308,   650,   835,  1061,  -954,  -954,  -954,  -954,
    4234,  4234,  4234,  4234,  4234,  4234,  4234,  4234,  4234,  4234,
    4234,  4234,  4234,  -954,  4192,  4192,  4192,  4192,  4192,  4192,
    4192,  4192,  4234,  4234,  4234,  4234,  4192,  4234,  -954,  2503,
    6461,  4915,  2311,  1251,  3373,  1403,  4387,  1715,  4409,  1762,
    1273,  5099,  4431,  1912,  4453,  2130,  4475,  2172,  4497,  2214,
    4519,  2335,  4541,  2446,  1425,  5117,  4563,  2472,   428,  4585,
    2615,  4607,  2665,  4629,  2727,  1462,  5135,  6359,  5153,  6377,
    5171,  6395,  5189,  4651,  2778,  4673,  2814,  4695,  2864,  4717,
    2886,  4739,  2968,  4761,  3243,    84,  -954,   786,  -954,     6,
    -245,  -195,  -954,  -954,  -954,  -954,  -954,    84,  -954,   957,
    -954,  -954,    84,   -10,   959,  -954,  -954,  -954,  -954,  -954,
    -954,   449,  4192,  4192,  4192,  4192,  4192,  4192,  4192,  4192,
    4192,  4192,  4192,  4192,  4192,  4192,   962,  4192,  4192,  4192,
    4192,  4192,  4192,  4192,  4192,  4192,  4192,  4192,  4192,  4192,
    -954,   963,  -954,   968,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  3280,  4192,  4192,  4192,  4192,  4192,  4192,   969,
     972,   973,  -954,  -954,   960,   976,   977,  -150,   982,   890,
     983,   985,  4192,  -358,   989,   990,   991,   993,   995,  4192,
    4192,  4192,   365,  -954,   471,   193,   997,   129,  1000,  1002,
    1004,  1005,   850,   998,  1012,  1013,  1017,   828,    84,  5207,
    -954,  -954,  5225,  -954,  -954,  -954,  -954,  1018,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,   -33,  -954,  2507,  1020,  1021,
    4192,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -195,
      84,  -233,  -195,   365,  -195,  -954,  -195,  -954,  -195,  -954,
    -954,  -954,  1022,   -84,  -195,  -954,  1023,  -954,  5243,  -954,
     495,  1027,   517,  5261,  5279,  5297,  -954,  3302,  -954,  1031,
    1032,  1033,  4192,  5315,  1034,   -35,  1038,  -195,  -176,  -358,
    1040,   365,  -954,  5333,  1043,  1042,  -954,  -954,  1048,  1060,
    -954,  -954,  -954,  1049,  1052,  -954,  -954,  -954,  -954,  4192,
    -954,  -954,  -954,  -954,  5351,  -954,     1,   690,  1055,   332,
    1056,  2507,  1057,  1059,  1063,     5,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  1064,   679,     8,  1066,   657,  -234,   658,
    1071,  1019,  1076,  1075,  -954,  -954,  1077,  -954,  -954,  -195,
    -954,  1079,  -954,  5369,  -954,   569,  5387,  5405,  5423,  -954,
    5441,  1086,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    1089,  1090,   -36,  1091,  4192,  -954,  -954,  -954,  -954,  -954,
    1092,  5459,  -954,  1093,  1094,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  4192,  -954,  1096,  -954,  5477,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  4192,  -954,  -954,  5495,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  1098,  1100,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  6447,  6447,  6447,  6447,  4192,
    -954,  -954,  5513,  -954,  -954,  -954,  4192,  -954,  -954,  -195,
    -954,  1101,  -954,  5531,  -954,   590,  -954,  5549,  5567,  5585,
    -954,  -954,  -954,  -954,  1103,  1105,  1106,  -954,  -954,  -954,
     745,  1108,  1110,   -47,  -243,  1111,     4,  -954,   600,  1113,
    1117,   752,  1119,  1120,  1347,  1121,   -47,  -128,  1122,   -36,
    1125,  1132,   -20,  1138,  1139,  3633,  3675,   -26,  1141,   835,
    -954,  -954,   620,  -195,  4192,  4192,  -195,  -954,  1143,  1144,
    -954,  1145,  3067,  3532,  4192,  4192,   -36,  -954,  1146,  1148,
     -24,  -954,  -954,  -954,  1153,  6461,  1526,  1556,  1556,  1556,
    1556,  1556,  1556,  -213,   -45,  -213,   -45,   641,  -199,   641,
    -199,   641,  -199,  2936,  1127,  2265,  2265,  2265,  2265,  2265,
    2265,  -213,   -45,  -213,   -45,   641,  -114,   641,  -106,   642,
     641,  -102,   662,  -954,   645,  -954,  -954,  -954,  -954,  -954,
    -954,  4234,  4192,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  4234,  4234,  -954,  -954,   651,
    -954,  -954,  -954,  -954,  -954,  -954,  4192,  4192,  4234,  4192,
    4234,  4192,  4234,  4192,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  -954,  -145,  -195,  -195,
    -195,  -195,  -195,  -195,  -954,  -195,  -195,  -954,    84,  -954,
    -954,    84,  -954,  1151,  1159,  1160,    -2,  -216,  4192,  -954,
    -954,  1161,  4935,  3329,  3351,  3440,  3485,  5603,  3507,  3628,
    3654,  3776,  3798,  3851,  5621,  3887,   647,  3928,  3950,  3972,
    5639,  5657,  5675,  5693,  3994,  4016,  4038,  4060,  4082,  4104,
    -954,  -954,    36,    36,   642,   642,   642,  5711,   653,   660,
     663,  -954,  -954,  -954,  -195,  -954,  -954,  1174,  -954,  -954,
    5729,  -954,  -954,    84,   665,  -954,    84,  -954,    84,  -954,
    5747,  5765,  5783,  -954,  4192,  1184,  -195,  -195,  -954,  -954,
    -954,  -954,  -954,   673,  -954,   674,  1189,  -954,  -954,   801,
    -954,  1191,  -954,  4192,  4192,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  2507,   680,  -954,  5801,  -954,
    -954,  -195,  1072,  -954,  -954,  -954,  -954,  -954,   -81,  1196,
    -954,  -954,  -954,  4192,  1201,  -954,  1203,  4192,  4192,  4192,
    4192,  -954,  -954,  -954,  6447,  4192,  -954,  1204,  -954,   694,
    -954,  -954,  -954,  -954,  -954,   695,  -954,  4192,  -954,  1044,
    -954,  4192,  -954,  -954,  5819,  4192,  -954,  1208,  1210,  1211,
    -954,  1212,  1213,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  1215,  1026,  1218,   838,  1222,  1223,  1170,  1224,  1225,
    1173,  1229,  -954,  -954,  -954,  -954,  -954,  -954,  4192,  1230,
    4192,  4192,  4192,  4192,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  6447,  -954,  4192,  -954,
    -954,  6447,  -954,  4192,  5837,  4192,  -954,  -954,  5855,  4192,
    5873,  -954,  -954,  4192,  1232,  4192,  4192,  4192,  -954,  -954,
    -954,  1235,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,   717,  -954,  -954,  1239,  -954,  -954,  1238,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  -954,   -23,  6447,   -17,
    6447,  -954,  -954,  -954,  -954,  4192,  -954,  6447,  6447,  -954,
    -954,  -954,  -954,  1241,  -195,  1242,  1244,  1245,  -954,  6447,
    1247,  -195,  1250,  1252,  1257,  -954,  6447,  6447,  6447,  -954,
    -954,  -954,  -954,  -954,  4192,  4783,  4145,  4805,  4187,  4827,
    4213,  5891,  5909,  4849,  4255,  4871,  4277,  4893,  4299,  1258,
    1259,  4192,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  1263,  -954,  -954,  5927,
     874,  -954,  -954,  4192,  -954,  4192,  1264,  1265,  1266,  -954,
    -954,  4192,  -954,  1267,  -954,  -954,  4192,  4192,  4192,  5945,
    -954,  -954,  -954,  4192,  1268,  -954,  1269,  -954,  6447,  6447,
    -954,  1270,  4192,  1147,  1274,  1275,  -954,  -954,  6447,  -954,
    -954,  6447,  5963,  6447,    86,  6447,  -954,  1276,  1277,  6447,
    1281,  6447,  4192,  6447,  -954,  -954,  -954,  -954,  -954,  -954,
    1282,   876,  1283,  -954,   904,  -954,   922,   761,  -954,   925,
    6447,  -954,  6447,  5981,  6447,  5999,  6447,  6447,  4192,  6447,
    4192,  6447,  4192,  6447,  -954,  6447,  6017,  6447,  -954,  1286,
    -954,  -954,  -954,  -954,  -954,  -954,  6447,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  -954,  6447,  -954,  -954,
    -954,  -954,  -954,   800,  4192,  4192,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  6035,  -954,  4192,  1304,  4321,  6053,
     802,   806,   808,  6447,  -954,  6447,  6071,  6447,  4192,  6089,
    -954,   917,  -954,  6447,  1308,  -954,  -954,  4192,  4192,  -954,
    -954,  -954,  6107,  1315,  1319,   935,  1323,  1325,  1326,  1331,
    4192,  4192,  6125,  6143,  6161,  4192,  -954,  4343,  4365,  4192,
    6447,  -954,  -954,  4192,  1332,  1334,  1335,  4192,  6447,  4192,
    1336,  -954,  6179,  6447,  4192,  -954,   974,  1338,   953,   955,
     964,   956,  6197,  6215,  4192,  4192,  4192,  6233,  -954,  -954,
    6447,  6447,  -954,  -954,  -954,  6251,  6269,  -954,  4192,  6447,
    1344,   981,  1348,  1354,  1355,  1357,  4192,  4192,  6447,  6447,
    6447,  4192,  4192,  4192,  6447,  -954,  1359,   999,  1001,   979,
    1003,  6447,   793,  6447,  6447,  6287,   -69,  1360,  1362,  1363,
    1367,  4192,  4192,  1369,  1371,  1228,  1327,  1011,    -9,  6305,
    6447,  -954,  1328,  -954,  -954,  1378,  -954,  1379,  4192,  -954,
    1339,  -954,  6323,  -954,  4192,  6341,  4192,  6447
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint16 yydefact[] =
{
       0,     0,   653,   726,   725,     0,     0,     0,     0,     0,
       0,     0,   683,   684,   685,     0,     0,     0,     0,     0,
     694,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   707,   708,   709,   695,     0,     0,     0,     0,
       0,   484,   485,   259,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   434,   263,     0,     0,   235,    95,     0,
      78,     0,     0,     0,     0,   248,   296,     0,     0,     0,
     289,     0,     0,     0,     0,   114,     0,     0,   433,   107,
     432,    93,     0,    91,    90,    76,    77,    81,   112,     0,
      96,     0,   390,     0,     0,    75,    86,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   757,   758,   759,
     760,     0,     0,     0,     0,   761,   762,   763,   764,     0,
     482,     0,   483,   486,   487,     0,   317,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   478,     0,     0,     0,     0,    26,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   104,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   210,     0,     0,
       0,   206,     0,     0,     0,   212,     0,   120,     0,     0,
       0,     0,     0,     0,   631,   630,   488,   491,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   741,   742,
     743,     0,     0,     0,     0,     0,   765,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   778,   779,
     780,   766,     0,     0,     0,     0,     0,    41,     0,    32,
       0,     0,     0,     0,     0,     0,    31,     0,    89,    28,
      44,   437,    56,   105,     0,     0,     0,     0,     0,   108,
       0,   180,     0,   137,   123,     0,     0,     0,     0,     0,
       0,     0,     0,   430,   429,   426,   428,   427,     0,     0,
       0,     0,     0,     0,   458,    80,   258,     0,     0,     0,
       0,     0,   465,   479,   481,   241,   256,     0,   250,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   725,     0,     0,   272,     0,     0,   291,     0,     0,
       0,     0,     0,   423,   424,     0,     0,    55,    46,     0,
       0,     0,    45,     0,     0,   425,     0,     0,     0,     0,
       0,     0,     0,     0,   357,   350,   100,     0,    97,    98,
      99,   101,     0,     0,     0,     0,   118,     0,     0,     0,
       0,    92,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   396,   581,   582,   431,   435,     0,
       0,     0,     0,     0,    60,     0,     0,     0,     0,     0,
      59,     0,   453,     0,     0,     0,     0,    87,   117,     0,
       0,     0,     0,     0,   364,   358,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   365,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   388,   376,
     448,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   349,   342,     0,    29,
       0,     0,     0,     0,   421,     0,     0,    42,     0,    57,
       0,     0,     0,     0,     0,     0,     0,    79,     0,     0,
       0,     0,     0,     0,   341,   334,   489,   490,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   323,   653,     0,   724,   654,   722,   654,
       0,     0,     1,     8,     7,    25,     9,    14,    10,    13,
      18,    15,    16,    24,    23,    22,    19,    21,    20,    17,
      11,    12,     0,     0,     0,     0,   477,   502,     3,     4,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     6,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     5,   651,
     652,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   205,     0,   629,     0,
       0,     0,   628,   626,   627,   194,   208,     0,   218,     0,
     214,   217,     0,     0,     0,   492,   495,   497,   494,   493,
      85,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      37,   637,   638,   636,    38,    35,    36,   633,   632,    34,
     797,   795,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   262,   264,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   178,     0,   176,   174,     0,     0,
       0,     0,     0,   124,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     232,   234,     0,   237,   238,   233,   236,     0,   462,   463,
     460,   461,   459,   635,   634,     0,   469,     0,     0,     0,
       0,   472,   470,   466,   587,   589,   588,   586,   590,   471,
     749,   750,   751,   752,   753,   754,   755,   756,   480,     0,
       0,     0,     0,     0,     0,   243,     0,   255,     0,   249,
     299,   306,     0,     0,     0,   303,     0,   307,     0,   301,
       0,     0,     0,     0,     0,     0,   297,     0,   106,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   267,     0,     0,     0,   277,   276,     0,     0,
     273,   290,    52,     0,     0,    51,    50,    49,    48,     0,
     103,   115,   352,   351,     0,   354,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   353,   647,   649,   648,
     646,   356,   102,     0,   265,     0,     0,     0,     0,     0,
       0,   157,     0,     0,    88,   113,     0,   393,   398,     0,
     397,     0,   399,     0,   394,     0,     0,     0,     0,   391,
       0,     0,    64,    74,    66,    67,    63,    62,    65,    68,
       0,     0,     0,     0,     0,   456,   454,   457,   455,   359,
       0,     0,   362,     0,     0,   370,   372,   375,   373,   369,
     368,   371,   366,     0,   384,     0,   377,     0,   386,   382,
     379,   381,   380,   383,   378,   451,   449,   452,   450,   445,
     444,   447,   446,     0,   345,   343,     0,   346,   344,   348,
     119,    30,   261,   122,     0,     0,   231,   422,   242,   295,
      43,   389,    58,   116,   316,   440,   439,   442,   441,     0,
     337,   335,     0,   338,   336,   340,     0,   320,   327,     0,
     324,     0,   328,     0,   321,     0,   326,     0,     0,     0,
     325,   318,   723,   796,     0,     0,     0,   501,   577,   566,
       0,     0,     0,     0,     0,     0,     0,   500,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   499,
     534,   536,     0,     0,     0,     0,     0,   593,     0,     0,
     591,     0,     0,     0,     0,     0,     0,   592,     0,     0,
       0,   498,   514,   503,   721,   654,   720,   714,   715,   716,
     717,   718,   719,   656,   654,   660,   654,   664,   654,   668,
     654,   674,   654,   794,   793,   787,   788,   789,   790,   791,
     792,   657,   654,   661,   654,   665,   654,   669,   654,   732,
     672,   654,   727,   676,   675,   677,   735,   678,   736,   679,
     737,     0,     0,   681,   739,   682,   740,   686,   744,   687,
     745,   688,   746,   689,   747,     0,     0,   696,   767,   697,
     698,   769,   699,   770,   700,   771,     0,     0,     0,     0,
       0,     0,     0,     0,   705,   776,   706,   777,   710,   781,
     711,   782,   712,   783,   713,   784,   211,     0,     0,     0,
       0,     0,     0,     0,   195,     0,     0,   204,     0,   207,
     216,     0,   213,     0,     0,     0,     0,     0,     0,   219,
     496,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      39,    40,   728,   729,   730,   731,   733,     0,     0,     0,
       0,   125,   183,   187,     0,   186,   169,     0,   164,   188,
       0,   184,   165,     0,     0,   168,     0,   167,     0,   166,
       0,     0,     0,   181,     0,     0,     0,     0,   147,   143,
     139,   141,   144,     0,   149,     0,     0,   142,   151,     0,
     152,     0,   138,     0,     0,   464,   467,   468,   606,   611,
     610,   612,   601,   607,   597,   620,   616,   621,   617,   622,
     618,   598,   600,   613,   604,   614,   602,   615,   623,   619,
     603,   605,   609,   608,   599,     0,     0,   476,     0,   244,
     257,     0,     0,   247,   251,   254,   253,   252,     0,     0,
     310,   302,   308,     0,     0,   309,     0,     0,     0,     0,
       0,   286,   288,   287,   283,     0,   275,     0,   270,     0,
     269,   624,   625,   268,   284,     0,   292,     0,   278,   280,
     274,     0,    53,    54,     0,     0,   418,     0,     0,     0,
     405,   637,   636,   414,   415,   409,   408,   407,   406,   410,
     412,     0,     0,     0,     0,     0,     0,   159,     0,     0,
     161,     0,   158,   260,    82,   136,   395,   400,     0,     0,
       0,     0,     0,     0,    73,    72,    71,   639,   640,   645,
     644,   641,   642,   643,    69,    70,    61,   360,     0,   361,
     374,   367,   385,     0,     0,     0,   163,   135,     0,     0,
       0,   322,   329,     0,     0,     0,     0,     0,   578,   579,
     580,     0,   575,   574,   585,   584,   583,   572,   568,   567,
     576,     0,   569,   570,     0,   547,   542,     0,   563,   562,
     541,   540,   560,   548,   544,   546,   545,   564,   551,   539,
     543,   595,   596,   594,   552,   549,   550,     0,   553,     0,
     554,   537,   538,   559,   535,     0,   515,   512,   511,   510,
     519,   513,   520,     0,     0,     0,     0,     0,   504,   506,
       0,     0,     0,     0,     0,   505,   507,   508,   509,   529,
     516,   532,   530,   531,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   198,   199,   196,   201,   200,   197,   202,   203,
     209,   215,   221,   226,   225,   220,     0,   224,   223,     0,
      83,   727,   734,     0,   768,     0,     0,     0,     0,   185,
     189,     0,   179,     0,   177,   175,     0,     0,     0,     0,
     145,   146,   148,     0,     0,   155,     0,   140,   240,   239,
     473,     0,     0,     0,     0,     0,   304,   311,   314,   298,
     300,   313,     0,   315,   785,   282,   271,     0,     0,   285,
       0,   281,     0,   355,   420,   419,   411,   416,   417,   413,
       0,     0,     0,   126,     0,   160,     0,     0,   162,     0,
     403,   392,   402,     0,   404,     0,   363,   387,     0,   347,
       0,   339,     0,   332,   319,   331,     0,   333,   573,     0,
     565,   561,   556,   555,   558,   557,   533,   521,   527,   525,
     523,   517,   522,   528,   526,   524,   518,   650,   680,   738,
     690,   692,   691,   693,     0,     0,   702,   773,   703,   774,
     704,   775,   229,   228,     0,   222,     0,     0,     0,     0,
       0,     0,     0,   192,   182,   191,     0,   193,     0,     0,
     154,     0,   474,   475,     0,   245,   305,     0,     0,   294,
     293,   279,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   571,     0,     0,     0,
     227,    84,   748,     0,     0,     0,     0,     0,   153,     0,
       0,   246,     0,   786,     0,   121,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   701,   772,
     230,    33,   109,   111,   110,     0,     0,   150,     0,    47,
       0,     0,     0,     0,     0,     0,     0,     0,   443,   438,
      94,     0,     0,     0,   312,   266,     0,     0,     0,     0,
       0,   401,     0,   330,   190,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   133,   170,     0,   131,     0,
     156,   129,   127,   134,   171,     0,   132,     0,     0,   128,
     172,   130,     0,   173,     0,     0,     0,   436
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,    17,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  -954,  -354,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,   -19,  -954,  -954,  1329,  1333,   -12,  1099,  -954,
    -954,   353,  -954,  -954,  -954,   369,  -954,  -954,  -954,  -837,
    -954,  -954,   697,    -3,  -297,   540,  -953,  -410,  -954,  -954,
      61,     0
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,   130,   131,   132,   133,   379,   356,   134,   748,  1238,
    1236,  1233,   732,   197,   647,  1168,   635,   652,  1171,   135,
     282,   136,   314,   793,   790,   137,   268,   138,   342,   831,
     139,   328,   140,   512,   141,   142,   143,   144,   145,   146,
     147,   393,   364,   148,   149,   150,   151,   206,   546,  1061,
    1039,  1040,  1017,  1007,   397,  1427,   779,  1062,  1454,  1295,
    1333,   645,   207,   699,   766,   694,  1394,   871,   152,   153,
     154,  1065
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
     155,   772,   851,   398,  1525,    41,    42,   842,  1509,  1346,
     963,  1387,  1431,  1360,  1364,  1165,  1309,   412,   763,  1565,
     175,    44,  1173,   315,  1356,   885,   334,   395,  1451,   649,
     834,   856,   909,  1331,  1266,   763,   542,  1158,   875,   296,
     450,  1461,  1424,  1492,  1622,  1388,   257,   372,  1389,  1390,
    1624,   335,  1786,   259,   819,   959,  1368,   316,  1773,   857,
     208,   377,  1428,   204,  1159,  1327,   205,   156,   858,   333,
    1429,  1432,   353,   274,  1527,   159,  1448,   371,   697,   157,
     160,  1444,   424,   985,   448,   763,   835,   466,   763,   209,
     698,   889,   161,   275,   317,   494,   162,   411,   158,   650,
     258,   204,   436,  1489,   205,   459,   989,   496,   163,   497,
     276,   413,   486,  -741,  1224,   876,   164,   418,  1787,   774,
     775,   776,   165,   204,   910,   498,   205,   517,   519,   521,
     378,   318,   357,  1160,   451,  1250,  -742,  1251,   204,   373,
    -743,   205,   859,  1174,  -765,  1510,   204,   297,   374,   205,
     860,   166,   861,   836,   691,   167,   837,   580,   581,   583,
     585,   587,   589,   591,   593,   595,   597,   599,   601,   603,
     605,   607,  1774,   610,   612,   298,   614,   616,   618,   620,
     622,   624,   626,   628,   630,   632,   634,  1441,   516,   518,
     520,  1425,   336,   414,   543,   375,  1301,   277,   911,   867,
     820,  1166,  1391,   821,   168,   299,   319,   877,   396,   320,
     260,   169,   321,   170,   862,   171,   452,   863,   579,   172,
     582,   584,   586,   588,   590,   592,   594,   596,   598,   600,
     602,   604,   606,  1369,   609,   611,   173,   613,   615,   617,
     619,   621,   623,   625,   627,   629,   631,   633,  1678,   713,
    1445,  1528,   288,   700,   701,   702,   204,   337,  1332,   205,
     415,   204,   174,   756,   205,   419,   964,   204,  1452,   177,
     205,  1161,   358,  1453,   764,   749,   753,   289,   752,   912,
    1175,  1156,   178,   453,   359,   322,   204,   360,   179,   205,
    -778,   764,   765,  1169,   878,   799,   864,   180,  1172,   323,
     560,   561,   204,   562,   181,   205,  1446,  -779,  1426,   816,
     343,   344,   278,    69,    70,   279,   576,   577,   809,  1392,
     808,  1162,  1511,   841,  -780,   813,   814,   815,   204,   817,
     777,   205,  1176,   182,   868,   823,  -766,   833,   778,   420,
     839,   764,   183,   184,   764,  1267,   692,   263,   848,   261,
     280,   290,  1462,  1365,  1493,  1623,   185,   855,   186,   854,
     869,  1625,   264,   262,   693,   361,   324,   544,   265,   545,
     865,   843,   266,  1347,   899,   870,   267,  1361,   269,   338,
     283,  1393,  1433,   416,   339,   894,   295,   893,  1245,   879,
     896,   897,   898,   325,  1262,   362,   913,   300,   329,   900,
     376,  -730,   577,   330,   907,   326,   454,   340,  1004,  -731,
     577,   119,   120,  -733,   577,   354,  1177,   421,   380,   922,
     422,   921,   965,   291,   327,   844,   292,  1348,   381,   930,
     880,   417,  1351,   522,   523,  1163,  1300,  1246,   938,   524,
     937,   121,   122,   399,   341,   525,   526,   123,   124,   527,
     881,   528,   529,  1247,   530,   531,   532,  1178,  1560,   363,
     957,   204,   956,   533,   205,  1526,   284,   534,   574,   575,
     576,   577,   608,   204,   636,   535,   205,   536,   537,   637,
     538,   285,   975,   976,   977,   978,   281,   646,   983,   204,
     982,   539,   205,  1001,   540,   648,   301,   541,   361,   651,
     653,   548,   994,   460,   993,   549,  1310,   210,   997,   998,
     999,     3,     4,   211,   212,   213,   214,   215,   216,   217,
     218,   219,   220,   221,   222,   223,   224,   225,   226,   227,
     228,   229,   230,   231,   232,   233,   234,   176,   235,   236,
     237,   238,   239,   240,   241,   242,   243,   244,   245,   705,
     706,   576,   707,   654,  1005,   655,   656,  1006,  1074,  1076,
    1078,  1080,  1082,   658,  1083,  1084,  1085,  1086,  1087,  1088,
    1089,  1090,  1092,  1094,  1096,  1098,  1099,  1101,   659,   660,
     661,   204,   663,   423,   205,   488,   662,   664,   802,   564,
     565,   566,   567,   568,   569,   570,   571,   703,   704,   705,
     706,   576,   707,   293,   204,   665,   666,   205,   667,   468,
     668,  1064,  1066,  1067,  1068,  1069,  1070,  1071,  1072,  1073,
    1075,  1077,  1079,  1081,   692,   669,   670,   671,   204,   672,
     673,   205,   674,  1091,  1093,  1095,  1097,   690,  1100,   675,
     676,   677,  1352,   461,   366,   469,   678,   767,   679,   680,
     681,   286,   695,   682,   683,   462,  1009,   684,   463,   685,
     638,   686,  1182,  1183,  1184,  1185,  1186,  1187,  1188,  1189,
    1190,  1191,  1192,  1193,  1194,  1195,   687,  1197,  1198,  1199,
    1200,  1201,  1202,  1203,  1204,  1205,  1206,  1207,  1208,  1209,
     688,   689,   696,   709,   301,   710,   711,   714,   361,   204,
     332,   470,   205,  1212,  1213,  1214,  1215,  1216,  1217,   733,
     471,   754,   757,  1243,   758,   760,   472,    67,   761,   762,
    1231,   499,  1230,   789,   734,   489,   639,   367,   791,  1240,
    1241,  1242,   792,   794,   796,   798,   361,   490,   800,   801,
     491,  1010,   368,   256,   803,   187,   806,   807,   640,   305,
    1011,   810,   811,   473,   474,   500,   287,   475,   294,   204,
     302,   812,   205,   476,   818,   824,   464,   198,   838,   352,
    1298,   355,   840,   365,  1304,   845,   199,   846,   847,   852,
     853,   204,   866,   872,   205,   873,   187,   188,   874,   394,
     883,   768,   501,   882,   410,   884,   886,   641,   425,   435,
     449,   804,   477,   467,   887,   888,   892,   895,   891,   901,
     487,   495,  1336,   189,  1520,   902,   903,  1521,   361,   306,
     904,   307,  1324,   513,   190,   905,  1334,   308,   188,   502,
     465,   906,  1012,   204,   908,   915,   205,   916,   917,   918,
     919,   923,   191,  1018,   924,   192,   920,   547,   492,  1344,
     925,   926,   927,   735,   189,   928,   929,   931,   200,   369,
     932,   934,   642,   935,   193,   190,   936,   643,   309,   939,
     644,   736,   940,   941,   942,   201,   943,   194,   310,  1542,
     944,  1013,  1544,   191,  1545,  1014,   192,   945,   946,   769,
     947,   948,   949,   737,  1019,   311,  1020,   950,   657,  1015,
     951,   952,   954,   955,   503,   193,   958,   504,   960,   961,
     505,   738,   493,   962,  1396,   312,   966,   967,   194,   968,
     969,   970,   202,   971,   972,   203,  1021,   270,   973,   974,
     770,   980,   370,  1401,   478,  1022,   981,   984,   987,   988,
     991,   992,   107,   108,   109,   110,   995,   996,  1000,   176,
    -797,   204,  1129,  1404,   205,   479,   115,   116,   117,   118,
     739,  -795,   480,  1157,   712,  1170,  1221,  1180,  1016,  1181,
    1196,  1210,   204,   750,   751,   205,  1211,  1218,   755,  1408,
    1219,  1220,   740,   506,  1222,  1223,  1410,   759,  1226,  1228,
    1227,  1244,   271,  1229,   771,  1232,   773,  1235,  1234,  1237,
     547,  1239,   195,  1249,  1257,  1261,   741,   795,  1252,   797,
    1253,  1566,  1254,  1255,   805,  1314,   253,  1023,   507,  1024,
    1258,  1259,   254,  1256,   382,  1260,  1265,   255,  1296,  1297,
    1308,  1312,   313,  1315,   832,  1458,  1460,  1316,  1025,  1321,
    1322,  1323,  1326,   195,  1467,  1468,  1329,   196,  1335,  1338,
    1339,   850,  1479,  1486,  1487,  1488,  1340,  1342,   383,   481,
    1343,   742,  1349,  1350,  1355,  1357,  1026,  1358,  1363,  1042,
    1027,  1359,  1362,   715,  1366,   743,   744,  1367,  1370,  1371,
    1372,   890,  1373,  1374,  1028,  1375,   272,  1377,   196,  1379,
    1043,   508,  1044,  1045,  1384,   384,   437,  1385,  1386,  1395,
    1397,  1399,  1400,   509,  1402,   745,  1406,   716,  1407,  1412,
    1414,  1418,  1496,  1419,  1420,  1046,  1422,  1421,  1423,  1430,
    1434,  1435,   510,  1436,  1437,  1498,  1500,  1438,  1439,  1442,
    1447,   204,   385,  1449,   205,  1029,  1501,  1502,  1450,  1504,
    1465,  1506,   746,  1508,   717,  1047,  1455,  1456,   273,  1463,
     747,  1470,  1471,  1472,  1490,   511,  1491,   562,   707,  1522,
    -734,  1048,   718,  1030,  1494,  1031,  -768,  1523,  1524,  1530,
     438,  1534,  1495,  1536,   780,   781,   782,   783,  1529,  1049,
    1537,   719,  1540,  1538,   720,  1543,  1497,  1499,   784,   785,
     786,   787,  1550,  1553,  1554,  1555,  1556,  1557,   990,  1503,
    1561,  1505,  1771,  1507,  1567,  1564,  1032,   386,  1033,  1569,
     387,  1570,  1576,   388,  1577,  1578,  1584,  1580,  1585,  1586,
    1587,  1588,  1034,  1589,  1035,  1590,  1591,  1592,  1593,  1594,
    1036,  1595,  1596,  1597,  1598,  1599,   439,  1619,  1601,  1008,
    1614,  1041,  1063,  1618,  1549,  1620,  1621,  1657,   440,  1627,
    1629,   441,  1630,  1631,   204,  1632,   721,   205,  1634,   722,
    1635,  1037,   723,  1558,  1559,  1636,  1652,  1653,   724,  1655,
    1038,  1684,  1660,  1661,  1662,  1664,  1670,  1671,  1672,  1050,
    1674,  1688,  1675,  1676,  1679,  1680,   389,  1681,  1683,  1783,
    1710,  1685,   204,  1686,  1696,   205,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
    1051,  1687,  1701,  1568,  1689,  -748,  1711,  1571,  1572,  1573,
    1574,  1715,  1704,  1052,   442,  1575,  1705,  1716,  1706,   361,
    1717,  1718,  1053,  1719,  1720,   725,  1164,  1579,  1167,  1721,
    1732,  1581,  1733,  1734,  1737,  1583,  1741,  1740,  1742,  1743,
    1179,  1745,  1755,  1744,  1756,  1440,  1757,   204,  1054,  1055,
     205,  1056,  1758,  1759,  1057,  1760,   726,  1766,  1775,   727,
    1776,  1777,  1767,  1769,  1768,  1778,  1770,  1781,  1600,  1782,
    1602,  1603,  1604,  1605,  1785,   443,  1790,  1791,  1784,  1789,
     444,   445,  1464,   303,   390,  1443,     0,   304,  1606,  1354,
    1793,     0,   788,  1607,     0,  1609,   391,     0,     0,  1611,
       0,     0,     0,  1613,  1225,  1615,  1616,  1617,     0,     0,
       0,     0,     0,   446,     0,   392,     0,     0,     0,     0,
       0,     0,  1248,  1058,     0,     0,     0,   728,     0,  1059,
       0,   447,     0,   729,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1341,   730,  1268,  1269,  1270,  1271,
       0,  1272,  1273,  1274,     0,  1626,  1275,  1276,  1277,  1278,
    1279,  1280,     0,     0,   731,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1299,  1060,  1302,  1303,
       0,  1305,     0,  1306,  1637,  1307,     0,     0,     0,     0,
       0,  1311,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1654,     0,     0,  1281,     0,     0,     0,     0,     0,
       0,     0,  1328,     0,  1330,     0,     0,     0,     0,     0,
       0,     0,     0,  1658,     0,  1659,     0,     0,     0,  1282,
       0,  1663,     0,     0,     0,     0,  1665,  1666,  1667,     0,
       0,     0,     0,  1669,     0,  1283,  1353,     0,     0,     0,
       0,     0,  1673,   564,   565,   566,   567,   568,   569,   570,
     571,   703,   704,   705,   706,   576,   707,     0,     0,     0,
       0,     0,  1682,     0,     0,     0,  1376,     0,   564,   565,
     566,   567,   568,   569,   570,   571,   572,   573,   574,   575,
     576,   577,     0,     0,     0,     0,     0,     0,  1692,  1003,
    1693,     0,  1694,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1284,  1285,  1286,  1287,  1288,  1289,  1290,  1291,
    1292,  1293,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,  1697,  1698,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1700,   551,   552,   553,
     554,   555,   556,   557,   558,   559,   560,   561,  1708,   562,
       0,     0,     0,     0,  1294,     0,     0,  1712,  1713,     0,
       0,     0,     0,     0,     0,     0,  1411,     0,     0,     0,
    1722,  1723,     0,     0,     0,  1727,     0,     0,     0,  1730,
       0,     0,     0,  1731,     0,     0,     0,  1735,     0,  1736,
       0,     0,     0,     0,  1739,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1748,  1749,  1750,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1041,     0,  1754,     0,
    1466,     0,     0,  1469,     0,     0,  1761,  1762,     0,  1478,
    1485,  1763,  1764,  1765,   564,   565,   566,   567,   568,   569,
     570,   571,   572,   573,   574,   575,   576,   577,     0,     0,
       0,  1779,  1780,     0,     0,  1104,   550,   551,   552,   553,
     554,   555,   556,   557,   558,   559,   560,   561,  1792,   562,
       0,     0,     0,  1111,  1795,     0,  1797,    -2,     1,     0,
       2,     0,     0,     0,     3,     4,     5,     6,     7,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
       0,    29,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    39,     0,     0,     0,     0,     0,     0,    40,     0,
       0,     0,     0,    41,    42,  1512,  1513,  1514,  1515,  1516,
    1517,     0,  1518,  1519,     0,    43,     0,     0,     0,    44,
       0,     0,     0,     0,    45,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    46,     0,     0,     0,     0,    47,
       0,    48,     0,     0,    49,     0,    50,    51,    52,     0,
       0,     0,     0,     0,     0,    53,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   575,   576,   577,
       0,  1539,     0,     0,     0,     0,     0,  1106,   550,   551,
     552,   553,   554,   555,   556,   557,   558,   559,   560,   561,
      54,   562,     0,  1551,  1552,  1125,    55,     0,     0,    56,
       0,     0,     0,     0,     0,    57,     0,     0,     0,     0,
       0,     0,     0,    58,     0,   550,   551,   552,   553,   554,
     555,   556,   557,   558,   559,   560,   561,    59,   562,     0,
       0,     0,  1136,     0,     0,    60,    61,     0,     0,     0,
       0,    62,     0,    63,    64,     0,    65,     0,  1563,     0,
       0,     0,    66,     0,     0,     0,     0,     0,     0,     0,
      67,     0,     0,     0,     0,     0,    68,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    69,    70,     0,    71,
       0,   552,   553,   554,   555,   556,   557,   558,   559,   560,
     561,    72,   562,     0,     0,     0,     0,    73,     0,     0,
       0,    74,     0,     0,     0,     0,     0,     0,    75,     0,
      76,  -798,  -798,  -798,  -798,  -798,  -798,   558,   559,   560,
     561,   210,   562,     0,     0,     3,     4,   211,   212,   213,
     214,   215,   216,   217,   218,   219,   220,   221,   222,   223,
     224,   225,   226,   227,   228,   229,   230,   231,   232,   233,
     234,    77,   235,   236,   237,   238,   239,   240,   241,   242,
     243,   244,   245,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    78,     0,     0,     0,     0,    79,     0,
       0,     0,     0,     0,     0,     0,    80,     0,    81,    82,
       0,    83,     0,     0,     0,     0,     0,     0,     0,    84,
      85,    86,     0,     0,     0,    87,     0,     0,     0,    88,
       0,     0,     0,    89,     0,     0,    90,     0,     0,     0,
       0,  1628,    91,    92,    93,     0,     0,     0,  1633,     0,
      94,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    95,    96,     0,    97,     0,    98,     0,     0,     0,
      99,   100,     0,     0,     0,     0,     0,   101,     0,     0,
     102,   103,   104,     0,     0,     0,     0,     0,   564,   565,
     566,   567,   568,   569,   570,   571,   572,   573,   574,   575,
     576,   577,     0,   105,   106,   107,   108,   109,   110,  1108,
       0,     0,     0,     0,   111,   112,   113,     0,   114,   115,
     116,   117,   118,     0,     0,     0,     0,     0,     0,   119,
     120,     0,     0,     0,   332,   564,   565,   566,   567,   568,
     569,   570,   571,   572,   573,   574,   575,   576,   577,     0,
       0,    67,     0,     0,     0,     0,  1110,     0,     0,   121,
     122,     0,     0,     0,     0,   123,   124,   125,   126,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   127,
       0,     0,     0,     0,     0,   128,   -27,   210,     0,     0,
     129,     3,     4,   211,   212,   213,   214,   215,   216,   217,
     218,   219,   220,   221,   222,   223,   224,   225,   226,   227,
     228,   229,   230,   231,   232,   233,   234,     0,   235,   236,
     237,   238,   239,   240,   241,   242,   243,   244,   245,     0,
     210,     0,     0,     0,     3,     4,   211,   212,   213,   214,
     215,   216,   217,   218,   219,   220,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
       0,   235,   236,   237,   238,   239,   240,   241,   242,   243,
     244,   245,     0,     0,     0,   400,     0,     0,     0,     0,
       0,     0,     0,     0,   401,   564,   565,   566,   567,   568,
     569,   570,   571,   572,   573,   574,   575,   576,   577,     0,
       0,   402,     0,     0,     0,     0,  1114,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   426,     0,
       0,     0,     0,     0,   210,     0,     0,   427,     3,     4,
     211,   212,   213,   214,   215,   216,   217,   218,   219,   220,
     221,   222,   223,   224,   225,   226,   227,   228,   229,   230,
     231,   232,   233,   234,     0,   235,   236,   237,   238,   239,
     240,   241,   242,   243,   244,   245,   403,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   107,   108,   109,   110,
     404,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     115,   116,   117,   118,     0,     0,     0,    67,     0,     0,
       0,     0,     0,   455,   456,     0,     0,     0,     0,   428,
       0,   429,     0,     0,     0,   405,     0,     0,     0,     0,
       0,   246,     0,   332,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   406,   457,   458,     0,     0,   247,   407,
      67,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     253,     0,     0,     0,     0,     0,   254,     0,   430,     0,
       0,   255,     0,     0,   550,   551,   552,   553,   554,   555,
     556,   557,   558,   559,   560,   561,   431,   562,     0,     0,
     204,     0,   432,   205,     0,  1002,  1268,  1269,  1270,  1271,
       0,  1272,  1273,  1274,     0,   248,  1275,  1276,  1277,  1278,
    1279,  1280,     0,   564,   565,   566,   567,   568,   569,   570,
     571,   572,   573,   574,   575,   576,   577,   249,     0,     0,
       0,     0,     0,   204,  1116,     0,   205,     0,     0,     0,
       0,     0,     0,     0,    67,     0,     0,     0,     0,     0,
       0,   408,     0,     0,  1281,   564,   565,   566,   567,   568,
     569,   570,   571,   572,   573,   574,   575,   576,   577,     0,
       0,     0,     0,     0,     0,     0,  1118,   250,     0,  1282,
     251,     0,     0,     0,     0,     0,   252,     0,     0,     0,
       0,     0,     0,     0,   433,  1283,   409,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   575,   576,
     577,     0,     0,     0,     0,     0,     0,     0,  1120,     0,
       0,     0,     0,     0,     0,     0,   434,   204,     0,     0,
     205,     0,   107,   108,   109,   110,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   115,   116,   117,   118,
    -798,  -798,  -798,  -798,  -798,  -798,   703,   704,   705,   706,
     576,   707,  1284,  1285,  1286,  1287,  1288,  1289,  1290,  1291,
    1292,  1293,     0,     0,     0,   107,   108,   109,   110,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   115,
     116,   117,   118,     0,   550,   551,   552,   553,   554,   555,
     556,   557,   558,   559,   560,   561,   253,   562,     0,     0,
       0,     0,   254,     0,  1294,  1103,     0,   255,   564,   565,
     566,   567,   568,   569,   570,   571,   572,   573,   574,   575,
     576,   577,     0,     0,     0,     0,     0,     0,     0,  1122,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   253,
       0,     0,     0,     0,     0,   254,     0,     0,     0,     0,
     255,     0,     0,     0,     0,     0,     0,     0,     0,   107,
     108,   109,   110,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   115,   116,   117,   118,   210,     0,     0,
       0,     3,     4,   211,   212,   213,   214,   215,   216,   217,
     218,   219,   220,   221,   222,   223,   224,   225,   226,   227,
     228,   229,   230,   231,   232,   233,   234,     0,   235,   236,
     237,   238,   239,   240,   241,   242,   243,   244,   245,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     575,   576,   577,   253,     0,     0,     0,   345,     0,   254,
    1124,     0,     0,     0,   255,   564,   565,   566,   567,   568,
     569,   570,   571,   572,   573,   574,   575,   576,   577,     0,
       0,     0,     0,     0,     0,     0,  1128,     0,     0,     0,
       0,     0,     0,     0,   346,     0,   550,   551,   552,   553,
     554,   555,   556,   557,   558,   559,   560,   561,   210,   562,
       0,   347,     3,     4,   211,   212,   213,   214,   215,   216,
     217,   218,   219,   220,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,     0,   235,
     236,   237,   238,   239,   240,   241,   242,   243,   244,   245,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     210,     0,     0,     0,     3,     4,   211,   212,   213,   214,
     215,   216,   217,   218,   219,   220,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
     348,   235,   236,   237,   238,   239,   240,   241,   242,   243,
     244,   245,     0,     0,     0,     0,     0,    67,   564,   565,
     566,   567,   568,   569,   570,   571,   572,   573,   574,   575,
     576,   577,     0,     0,     0,     0,     0,     0,     0,  1131,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     349,     0,     0,   350,     0,     0,     0,     0,     0,   351,
       0,     0,     0,     0,     0,     0,     0,  1473,   564,   565,
     566,   567,   568,   569,   570,   571,   572,   573,   574,   575,
     576,   577,     0,     0,     0,     0,     0,     0,     0,  1133,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     204,     0,     0,   205,   825,     0,     0,     0,     0,     0,
       0,   332,     0,     0,     0,     0,     0,     0,   826,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    67,     0,
     564,   565,   566,   567,   568,   569,   570,   571,   572,   573,
     574,   575,   576,   577,     0,     0,     0,     0,     0,   827,
     828,  1135,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   332,     0,     0,  1474,     0,     0,     0,
     829,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      67,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   575,   576,   577,     0,     0,     0,     0,     0,
       0,     0,  1145,     0,     0,     0,     0,     0,     0,     0,
       0,   204,     0,  1475,   205,     0,  1476,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   575,   576,
     577,     0,     0,     0,     0,     0,     0,     0,  1147,     0,
       0,     0,   107,   108,   109,   110,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   115,   116,   117,   118,
       0,     0,     0,   204,     0,     0,   205,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   575,   576,
     577,     0,     0,     0,     0,     0,     0,   830,  1149,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     575,   576,   577,     0,     0,     0,     0,     0,     0,     0,
    1151,     0,     0,     0,     0,     0,   253,     0,     0,     0,
       0,     0,   254,     0,     0,     0,     0,   255,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1477,
     565,   566,   567,   568,   569,   570,   571,   703,   704,   705,
     706,   576,   707,   107,   108,   109,   110,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   115,   116,   117,
     118,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   575,   576,   577,     0,     0,     0,     0,     0,
       0,     0,  1153,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   107,   108,   109,   110,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   115,
     116,   117,   118,     0,     0,     0,     0,   253,     0,     0,
       0,     0,     0,   254,     0,   210,     0,     0,   255,     3,
       4,   211,   212,   213,   214,   215,   216,   217,   218,   219,
     220,   221,   222,   223,   224,   225,   226,   227,   228,   229,
     230,   231,   232,   233,   234,     0,   235,   236,   237,   238,
     239,   240,   241,   242,   243,   244,   245,     0,     0,   253,
       0,     0,     0,     0,     0,   254,     0,     0,     0,   210,
     255,     0,     0,     3,     4,   211,   212,   213,   214,   215,
     216,   217,   218,   219,   220,   221,   222,   223,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   233,   234,     0,
     235,   236,   237,   238,   239,   240,   241,   242,   243,   244,
     245,     0,  1480,     0,     0,     0,   210,     0,     0,     0,
       3,     4,   211,   212,   213,   214,   215,   216,   217,   218,
     219,   220,   221,   222,   223,   224,   225,   226,   227,   228,
     229,   230,   231,   232,   233,   234,     0,   235,   236,   237,
     238,   239,   240,   241,   242,   243,   244,   245,   210,     0,
       0,     0,     3,     4,   211,   212,   213,   214,   215,   216,
     217,   218,   219,   220,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,     0,   235,
     236,   237,   238,   239,   240,   241,   242,   243,   244,   245,
       0,     0,     0,     0,     0,     0,     0,     0,   332,     0,
       0,  1481,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    67,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   575,   576,   577,
       0,     0,     0,     0,     0,     0,     0,  1155,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1482,     0,
       0,  1483,   332,   564,   565,   566,   567,   568,   569,   570,
     571,   703,   704,   705,   706,   576,   707,     0,     0,    67,
       0,     0,     0,     0,  1003,   564,   565,   566,   567,   568,
     569,   570,   571,   703,   704,   705,   706,   576,   707,     0,
       0,     0,     0,     0,     0,     0,  1320,     0,   204,   332,
       0,   205,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,    67,     0,     0,     0,
       0,     0,     0,  1532,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,   332,     0,     0,     0,  1106,   550,   551,   552,   553,
     554,   555,   556,   557,   558,   559,   560,   561,    67,   562,
       0,     0,     0,     0,     0,     0,     0,  1105,     0,     0,
       0,   210,     0,     0,  1484,     3,   331,   211,   212,   213,
     214,   215,   216,   217,   218,   219,   220,   221,   222,   223,
     224,   225,   226,   227,   228,   229,   230,   231,   232,   233,
     234,     0,   235,   236,   237,   238,   239,   240,   241,   242,
     243,   244,   245,   564,   565,   566,   567,   568,   569,   570,
     571,   703,   704,   705,   706,   576,   707,     0,     0,     0,
       0,     0,     0,     0,  1108,     0,     0,     0,     0,     0,
     107,   108,   109,   110,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   115,   116,   117,   118,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,     0,     0,     0,     0,  1110,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,   107,   108,   109,   110,     0,     0,
       0,  1114,     0,     0,     0,     0,     0,     0,   115,   116,
     117,   118,     0,     0,   253,     0,     0,     0,     0,     0,
     254,   482,   483,     0,     0,   255,     0,     0,     0,  1457,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   107,   108,   109,   110,     0,     0,     0,     0,     0,
       0,     0,   484,   485,     0,   115,   116,   117,   118,     0,
       0,     0,     0,     0,   332,     0,     0,     0,   253,     0,
       0,  1459,     0,     0,   254,     0,     0,     0,     0,   255,
       0,    67,     0,   107,   108,   109,   110,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   115,   116,   117,
     118,   564,   565,   566,   567,   568,   569,   570,   571,   703,
     704,   705,   706,   576,   707,   253,     0,     0,     0,     0,
       0,   254,  1116,     0,     0,     0,   255,   564,   565,   566,
     567,   568,   569,   570,   571,   703,   704,   705,   706,   576,
     707,     0,     0,     0,     0,     0,     0,     0,  1118,     0,
       0,     0,     0,     0,     0,     0,     0,   253,     0,     0,
       0,     0,     0,   254,     0,   210,     0,     0,   255,     3,
       4,   211,   212,   213,   214,   215,   216,   217,   218,   219,
     220,   221,   222,   223,   224,   225,   226,   227,   228,   229,
     230,   231,   232,   233,   234,     0,   235,   236,   237,   238,
     239,   240,   241,   242,   243,   244,   245,   514,     0,     0,
       0,     3,     4,     5,     6,     7,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,   515,    26,    27,    28,     0,    29,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    39,   564,
     565,   566,   567,   568,   569,   570,   571,   703,   704,   705,
     706,   576,   707,     0,     0,     0,     0,     0,     0,     0,
    1120,   564,   565,   566,   567,   568,   569,   570,   571,   703,
     704,   705,   706,   576,   707,     0,     0,     0,     0,     0,
       0,     0,  1122,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   107,   108,   109,   110,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     115,   116,   117,   118,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,     0,     0,     0,     0,  1124,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   332,     0,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,    67,     0,     0,     0,     0,
     253,  1128,     0,     0,     0,     0,   254,     0,     0,     0,
       0,   255,     0,     0,     0,     0,     0,     0,     0,     0,
     332,   564,   565,   566,   567,   568,   569,   570,   571,   703,
     704,   705,   706,   576,   707,     0,     0,    67,     0,     0,
       0,     0,  1131,   564,   565,   566,   567,   568,   569,   570,
     571,   703,   704,   705,   706,   576,   707,     0,     0,     0,
       0,     0,     0,     0,  1133,   564,   565,   566,   567,   568,
     569,   570,   571,   703,   704,   705,   706,   576,   707,     0,
       0,     0,     0,     0,     0,     0,  1135,   564,   565,   566,
     567,   568,   569,   570,   571,   703,   704,   705,   706,   576,
     707,     0,     0,     0,     0,     0,     0,     0,  1145,   564,
     565,   566,   567,   568,   569,   570,   571,   703,   704,   705,
     706,   576,   707,     0,     0,     0,     0,     0,     0,     0,
    1147,   564,   565,   566,   567,   568,   569,   570,   571,   703,
     704,   705,   706,   576,   707,     0,     0,     0,     0,     0,
       0,     0,  1149,   564,   565,   566,   567,   568,   569,   570,
     571,   703,   704,   705,   706,   576,   707,     0,     0,     0,
       0,     0,     0,     0,  1151,   564,   565,   566,   567,   568,
     569,   570,   571,   703,   704,   705,   706,   576,   707,     0,
       0,     0,     0,     0,     0,     0,  1153,   564,   565,   566,
     567,   568,   569,   570,   571,   703,   704,   705,   706,   576,
     707,     0,     0,     0,     0,     0,     0,     0,  1155,     0,
     107,   108,   109,   110,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   115,   116,   117,   118,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,     0,     0,     0,     0,  1639,
       0,     0,   107,   108,   109,   110,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   115,   116,   117,   118,
     564,   565,   566,   567,   568,   569,   570,   571,   572,   573,
     574,   575,   576,   577,   253,     0,     0,     0,     0,     0,
     254,  1641,     0,     0,     0,   255,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   575,   576,   577,
       0,     0,     0,     0,     0,     0,     0,  1643,     0,     0,
       0,     0,     0,     0,     0,     0,   127,     0,     0,     0,
       0,     0,   128,     0,     0,     0,     0,   129,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,     0,     0,     0,     0,  1647,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,     0,     0,     0,
       0,  1649,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,     0,
       0,     0,     0,  1651,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,     0,     0,     0,     0,  1702,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,     0,     0,     0,     0,  1728,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,     0,     0,     0,     0,  1729,
     550,   551,   552,   553,   554,   555,   556,   557,   558,   559,
     560,   561,     0,   562,     0,     0,     0,     0,     0,     0,
       0,  1107,   550,   551,   552,   553,   554,   555,   556,   557,
     558,   559,   560,   561,     0,   562,     0,     0,     0,     0,
       0,     0,     0,  1109,   550,   551,   552,   553,   554,   555,
     556,   557,   558,   559,   560,   561,     0,   562,     0,     0,
       0,     0,     0,     0,     0,  1113,   550,   551,   552,   553,
     554,   555,   556,   557,   558,   559,   560,   561,     0,   562,
       0,     0,     0,     0,     0,     0,     0,  1115,   550,   551,
     552,   553,   554,   555,   556,   557,   558,   559,   560,   561,
       0,   562,     0,     0,     0,     0,     0,     0,     0,  1117,
     550,   551,   552,   553,   554,   555,   556,   557,   558,   559,
     560,   561,     0,   562,     0,     0,     0,     0,     0,     0,
       0,  1119,   550,   551,   552,   553,   554,   555,   556,   557,
     558,   559,   560,   561,     0,   562,     0,     0,     0,     0,
       0,     0,     0,  1121,   550,   551,   552,   553,   554,   555,
     556,   557,   558,   559,   560,   561,     0,   562,     0,     0,
       0,     0,     0,     0,     0,  1123,   550,   551,   552,   553,
     554,   555,   556,   557,   558,   559,   560,   561,     0,   562,
       0,     0,     0,     0,     0,     0,     0,  1127,   550,   551,
     552,   553,   554,   555,   556,   557,   558,   559,   560,   561,
       0,   562,     0,     0,     0,     0,     0,     0,     0,  1130,
     550,   551,   552,   553,   554,   555,   556,   557,   558,   559,
     560,   561,     0,   562,     0,     0,     0,     0,     0,     0,
       0,  1132,   550,   551,   552,   553,   554,   555,   556,   557,
     558,   559,   560,   561,     0,   562,     0,     0,     0,     0,
       0,     0,     0,  1134,   550,   551,   552,   553,   554,   555,
     556,   557,   558,   559,   560,   561,     0,   562,     0,     0,
       0,     0,     0,     0,     0,  1144,   550,   551,   552,   553,
     554,   555,   556,   557,   558,   559,   560,   561,     0,   562,
       0,     0,     0,     0,     0,     0,     0,  1146,   550,   551,
     552,   553,   554,   555,   556,   557,   558,   559,   560,   561,
       0,   562,     0,     0,     0,     0,     0,     0,     0,  1148,
     550,   551,   552,   553,   554,   555,   556,   557,   558,   559,
     560,   561,     0,   562,     0,     0,     0,     0,     0,     0,
       0,  1150,   550,   551,   552,   553,   554,   555,   556,   557,
     558,   559,   560,   561,     0,   562,     0,     0,     0,     0,
       0,     0,     0,  1152,   550,   551,   552,   553,   554,   555,
     556,   557,   558,   559,   560,   561,     0,   562,     0,     0,
       0,     0,     0,     0,     0,  1154,   550,   551,   552,   553,
     554,   555,   556,   557,   558,   559,   560,   561,     0,   562,
       0,     0,     0,     0,     0,     0,     0,  1638,   550,   551,
     552,   553,   554,   555,   556,   557,   558,   559,   560,   561,
       0,   562,     0,     0,     0,     0,     0,     0,     0,  1640,
     550,   551,   552,   553,   554,   555,   556,   557,   558,   559,
     560,   561,     0,   562,     0,     0,     0,     0,     0,     0,
       0,  1642,   550,   551,   552,   553,   554,   555,   556,   557,
     558,   559,   560,   561,     0,   562,     0,     0,     0,     0,
       0,     0,     0,  1646,   550,   551,   552,   553,   554,   555,
     556,   557,   558,   559,   560,   561,     0,   562,     0,     0,
       0,     0,     0,     0,     0,  1648,   550,   551,   552,   553,
     554,   555,   556,   557,   558,   559,   560,   561,     0,   562,
       0,     0,     0,     0,     0,     0,     0,  1650,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,     0,     0,  1102,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,     0,     0,  1531,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,   708,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,   822,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,   849,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,   914,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,   933,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,   953,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,   979,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,   986,   564,   565,   566,   567,   568,   569,   570,   571,
     572,   573,   574,   575,   576,   577,     0,     0,     0,  1112,
     564,   565,   566,   567,   568,   569,   570,   571,   572,   573,
     574,   575,   576,   577,     0,     0,     0,  1126,   564,   565,
     566,   567,   568,   569,   570,   571,   572,   573,   574,   575,
     576,   577,     0,     0,     0,  1137,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   575,   576,   577,
       0,     0,     0,  1139,   564,   565,   566,   567,   568,   569,
     570,   571,   572,   573,   574,   575,   576,   577,     0,     0,
       0,  1141,   564,   565,   566,   567,   568,   569,   570,   571,
     572,   573,   574,   575,   576,   577,     0,     0,     0,  1143,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1263,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1264,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1313,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1317,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,  1318,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1319,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1325,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1337,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1345,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,  1378,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1380,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1381,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1382,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1383,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,  1398,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1403,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1405,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1409,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1413,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,  1415,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1416,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1417,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1112,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1533,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,  1137,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1139,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1141,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1143,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1535,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,  1541,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1546,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1547,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1548,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1562,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,  1582,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1608,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1610,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1612,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1644,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,  1645,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1656,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1668,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1677,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1690,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,  1691,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1695,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1699,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1703,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1707,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,  1709,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1714,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1724,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1725,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1726,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,  1738,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1746,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1747,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1751,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1752,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,  1753,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1772,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1788,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1794,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1796,   550,   551,   552,   553,   554,   555,   556,   557,
     558,   559,   560,   561,     0,   562,     0,     0,     0,  1138,
     550,   551,   552,   553,   554,   555,   556,   557,   558,   559,
     560,   561,     0,   562,     0,     0,     0,  1140,   550,   551,
     552,   553,   554,   555,   556,   557,   558,   559,   560,   561,
       0,   562,     0,     0,     0,  1142,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   575,   576,   577,
       0,     0,   578,   550,   551,   552,   553,   554,   555,   556,
     557,   558,   559,   560,   561,     0,   562,     0,     0,   563,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,   564,   565,   566,   567,   568,   569,
     570,   571,   572,   573,   574,   575,   576,   577
};

static const yytype_int16 yycheck[] =
{
       0,   298,   356,     8,     6,    56,    57,     8,   153,     8,
       8,    47,     8,     8,     6,   260,   100,     6,   169,   100,
       6,    72,    32,    66,   861,   379,     8,   136,    48,   190,
     173,    63,     8,   209,    67,   169,    76,    31,     8,    68,
       6,    67,    89,    67,    67,    81,    46,    51,    84,    85,
      67,    33,    61,    87,   100,   465,   290,   100,   127,    91,
     329,     8,   305,   296,    58,   100,   299,   519,   100,    69,
     313,    67,    72,    56,   290,   523,  1029,    77,   436,   502,
     523,   209,   101,   493,   103,   169,   229,   106,   169,   358,
     448,   264,   523,    33,   137,   114,   523,    97,   521,   260,
       6,   296,   102,  1056,   299,   105,   264,   119,   523,   121,
      50,   100,   112,   515,   264,    85,   523,   100,   127,   256,
     257,   258,   523,   296,   100,   125,   299,   127,   128,   129,
      77,   174,   100,   127,   100,     6,   515,     8,   296,   143,
     515,   299,   174,   153,   515,   290,   296,   176,   152,   299,
     182,   523,   184,   296,   100,   523,   299,   157,   158,   159,
     160,   161,   162,   163,   164,   165,   166,   167,   168,   169,
     170,   171,   241,   173,   174,   204,   176,   177,   178,   179,
     180,   181,   182,   183,   184,   185,   186,  1024,   127,   128,
     129,   238,   174,   182,   234,   199,   429,   137,   174,    98,
     246,   446,   238,   249,   523,   234,   249,   177,   317,   252,
     244,   523,   255,   523,   246,   523,   182,   249,   157,   523,
     159,   160,   161,   162,   163,   164,   165,   166,   167,   168,
     169,   170,   171,   467,   173,   174,   523,   176,   177,   178,
     179,   180,   181,   182,   183,   184,   185,   186,   162,   268,
     378,   467,    75,   253,   254,   255,   296,   239,   434,   299,
     249,   296,   523,   282,   299,   100,   264,   296,   288,   523,
     299,   265,   240,   293,   425,   275,   279,   100,   278,   255,
     290,   635,   523,   249,   252,   328,   296,   255,   523,   299,
     515,   425,   426,   647,   264,   314,   328,   523,   652,   342,
     513,   514,   296,   516,   523,   299,   434,   515,   355,   328,
     196,   197,   252,   229,   230,   255,   515,   516,   321,   355,
     320,   315,   467,   342,   515,   325,   326,   327,   296,   329,
     467,   299,   342,   523,   233,   335,   515,   337,   475,   174,
     340,   425,   523,   523,   425,   378,   292,    29,   351,   383,
     290,   174,   378,   345,   378,   378,   523,   360,   523,   359,
     259,   378,     8,   397,   310,   333,   409,   407,     8,   409,
     402,   372,     8,   372,   393,   274,   137,   372,     6,   361,
      95,   417,   378,   372,   366,   388,     8,   387,   195,   359,
     390,   391,   392,   436,   748,   363,   372,   426,   523,   399,
     404,   515,   516,   216,   407,   448,   372,   389,   100,   515,
     516,   462,   463,   515,   516,   304,   426,   252,   152,   422,
     255,   421,   420,   246,   467,   426,   249,   426,     8,   432,
     400,     6,   100,     0,   519,   429,   790,   244,   441,   519,
     440,   492,   493,   448,   426,   519,   519,   498,   499,   519,
     420,   519,   519,   260,   519,   519,   519,   467,  1295,   427,
     463,   296,   462,   519,   299,   467,   181,   519,   513,   514,
     515,   516,     8,   296,     8,   519,   299,   519,   519,   360,
     519,   196,   482,   483,   484,   485,   426,    31,   491,   296,
     490,   519,   299,   512,   519,     8,   525,   519,   333,     8,
     177,   519,   505,   100,   504,   519,   803,     3,   508,   509,
     510,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,   523,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,   513,
     514,   515,   516,   372,   246,     8,     8,   249,   558,   559,
     560,   561,   562,     8,   564,   565,   566,   567,   568,   569,
     570,   571,   572,   573,   574,   575,   576,   577,     8,     6,
       8,   296,   523,   418,   299,   100,   521,   523,    51,   503,
     504,   505,   506,   507,   508,   509,   510,   511,   512,   513,
     514,   515,   516,   426,   296,   523,   523,   299,   523,    51,
     523,   550,   551,   552,   553,   554,   555,   556,   557,   558,
     559,   560,   561,   562,   292,   523,   523,   523,   296,   523,
     523,   299,   523,   572,   573,   574,   575,     8,   577,   523,
     523,   523,   310,   240,   140,    87,   523,   184,   523,   523,
     523,   366,     8,   523,   523,   252,     6,   523,   255,   523,
      89,   523,   662,   663,   664,   665,   666,   667,   668,   669,
     670,   671,   672,   673,   674,   675,   523,   677,   678,   679,
     680,   681,   682,   683,   684,   685,   686,   687,   688,   689,
     523,   523,     8,   520,   525,   520,   520,     8,   333,   296,
     196,   143,   299,   703,   704,   705,   706,   707,   708,     8,
     152,     8,   100,   732,     8,     8,   158,   213,     8,     8,
     723,    66,   722,   137,     8,   240,   155,   223,   137,   729,
     730,   731,   137,   137,   137,   137,   333,   252,     8,     8,
     255,    91,   238,    46,   207,    58,   100,     8,   177,     8,
     100,     8,     8,   195,   196,   100,    59,   199,    61,   296,
      63,     8,   299,   205,     6,   426,   363,    91,   426,    72,
     770,    74,     8,    76,   793,     8,   100,     8,     8,     8,
       8,   296,     8,     8,   299,   155,    58,   100,     6,    92,
     520,   328,   137,     8,    97,     6,   420,   226,   101,   102,
     103,   264,   244,   106,     8,     8,     8,     8,   100,   372,
     113,   114,   831,   126,  1168,     8,     6,  1171,   333,    78,
       8,    80,   822,   126,   137,     8,   829,    86,   100,   174,
     427,     8,   182,   296,     8,     8,   299,     8,     8,     8,
       8,   412,   155,     8,   372,   158,   100,   150,   363,   849,
       8,     8,     6,   137,   126,     8,     8,     8,   182,   355,
       6,     8,   291,   100,   177,   137,     8,   296,   127,     8,
     299,   155,     8,     8,     8,   199,     8,   190,   137,  1233,
       8,   231,  1236,   155,  1238,   235,   158,     8,     8,   426,
       8,     8,     8,   177,    59,   154,    61,     8,   201,   249,
       8,     8,     8,     8,   249,   177,     8,   252,     8,     8,
     255,   195,   427,     8,   914,   174,     8,     8,   190,     8,
       8,     8,   246,     8,     8,   249,    91,   199,     8,     8,
     467,     8,   428,   933,   376,   100,     8,     8,     8,     8,
     100,     8,   438,   439,   440,   441,     8,     8,     8,   523,
     515,   296,   524,   953,   299,   397,   452,   453,   454,   455,
     244,   515,   404,   177,   267,     8,     6,     8,   318,   520,
       8,     8,   296,   276,   277,   299,     8,     8,   281,   979,
       8,     8,   266,   328,     8,     8,   986,   290,     6,     6,
     100,   520,   264,     8,   297,     6,   299,     6,     8,     6,
     303,     6,   315,     6,     6,   177,   290,   310,     8,   312,
       8,  1308,     8,     8,   317,   520,   512,   182,   363,   184,
       8,     8,   518,   173,    66,     8,     8,   523,     8,     8,
       8,     8,   291,     6,   337,  1035,  1036,   520,   203,     8,
       8,     8,     8,   315,  1044,  1045,     8,   360,     8,     6,
       8,   354,  1052,  1053,  1054,  1055,     8,     8,   100,   501,
       8,   345,   372,     8,     8,     8,   231,     8,   389,     8,
     235,     8,     8,    66,     8,   359,   360,   420,   420,     8,
      61,   384,     6,     8,   249,     8,   358,     8,   360,   520,
      29,   436,    31,    32,     8,   137,   100,     8,     8,     8,
       8,     8,     8,   448,     8,   389,     8,   100,     8,     8,
     520,     8,  1112,     8,     8,    54,     8,   372,     8,     8,
     520,     8,   467,     6,   372,  1125,  1126,     8,     8,     8,
       8,   296,   174,     8,   299,   300,  1136,  1137,     6,  1139,
     520,  1141,   426,  1143,   137,    84,     8,     8,   420,     8,
     434,     8,     8,     8,     8,   500,     8,   516,   516,     8,
     515,   100,   155,   328,   502,   330,   515,     8,     8,     8,
     174,   524,  1111,   520,   438,   439,   440,   441,  1178,   118,
     520,   174,     8,   520,   177,   520,  1125,  1126,   452,   453,
     454,   455,     8,   520,   520,     6,   395,     6,   501,  1138,
     520,  1140,   409,  1142,     8,   133,   371,   249,   373,     8,
     252,     8,     8,   255,   520,   520,     8,   173,     8,     8,
       8,     8,   387,     8,   389,   199,     8,   389,     6,     6,
     395,    61,     8,     8,    61,     6,   240,   520,     8,   542,
       8,   544,   545,     8,  1244,     6,     8,   373,   252,     8,
       8,   255,     8,     8,   296,     8,   249,   299,     8,   252,
       8,   426,   255,  1263,  1264,     8,     8,     8,   261,     6,
     435,   395,     8,     8,     8,     8,     8,     8,     8,   218,
     133,   520,     8,     8,     8,     8,   328,     6,     6,    61,
     373,     8,   296,   389,     8,   299,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
     249,   389,     8,  1313,   389,   515,     8,  1317,  1318,  1319,
    1320,     6,   520,   262,   328,  1325,   520,     8,   520,   333,
     395,     8,   271,     8,     8,   328,   639,  1337,   641,     8,
       8,  1341,     8,     8,     8,  1345,     8,   373,   395,   394,
     653,   395,     8,   389,   373,     8,     8,   296,   297,   298,
     299,   300,     8,     8,   303,     8,   359,     8,     8,   362,
       8,     8,   373,   394,   373,     8,   373,     8,  1378,     8,
    1380,  1381,  1382,  1383,   373,   389,     8,     8,    61,    61,
     394,   395,  1039,    64,   436,  1026,    -1,    64,  1398,   859,
      61,    -1,   303,  1403,    -1,  1405,   448,    -1,    -1,  1409,
      -1,    -1,    -1,  1413,   717,  1415,  1416,  1417,    -1,    -1,
      -1,    -1,    -1,   427,    -1,   467,    -1,    -1,    -1,    -1,
      -1,    -1,   735,   372,    -1,    -1,    -1,   430,    -1,   378,
      -1,   445,    -1,   436,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   394,   448,   109,   110,   111,   112,
      -1,   114,   115,   116,    -1,  1465,   119,   120,   121,   122,
     123,   124,    -1,    -1,   467,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   789,   426,   791,   792,
      -1,   794,    -1,   796,  1494,   798,    -1,    -1,    -1,    -1,
      -1,   804,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1511,    -1,    -1,   167,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   825,    -1,   827,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1533,    -1,  1535,    -1,    -1,    -1,   192,
      -1,  1541,    -1,    -1,    -1,    -1,  1546,  1547,  1548,    -1,
      -1,    -1,    -1,  1553,    -1,   208,   859,    -1,    -1,    -1,
      -1,    -1,  1562,   503,   504,   505,   506,   507,   508,   509,
     510,   511,   512,   513,   514,   515,   516,    -1,    -1,    -1,
      -1,    -1,  1582,    -1,    -1,    -1,   889,    -1,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,    -1,    -1,    -1,  1608,   524,
    1610,    -1,  1612,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   275,   276,   277,   278,   279,   280,   281,   282,
     283,   284,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,  1644,  1645,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1656,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,  1668,   516,
      -1,    -1,    -1,    -1,   327,    -1,    -1,  1677,  1678,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   989,    -1,    -1,    -1,
    1690,  1691,    -1,    -1,    -1,  1695,    -1,    -1,    -1,  1699,
      -1,    -1,    -1,  1703,    -1,    -1,    -1,  1707,    -1,  1709,
      -1,    -1,    -1,    -1,  1714,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1724,  1725,  1726,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1039,    -1,  1738,    -1,
    1043,    -1,    -1,  1046,    -1,    -1,  1746,  1747,    -1,  1052,
    1053,  1751,  1752,  1753,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,  1771,  1772,    -1,    -1,   524,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,  1788,   516,
      -1,    -1,    -1,   520,  1794,    -1,  1796,     0,     1,    -1,
       3,    -1,    -1,    -1,     7,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      -1,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    44,    -1,    -1,    -1,    -1,    -1,    -1,    51,    -1,
      -1,    -1,    -1,    56,    57,  1158,  1159,  1160,  1161,  1162,
    1163,    -1,  1165,  1166,    -1,    68,    -1,    -1,    -1,    72,
      -1,    -1,    -1,    -1,    77,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    87,    -1,    -1,    -1,    -1,    92,
      -1,    94,    -1,    -1,    97,    -1,    99,   100,   101,    -1,
      -1,    -1,    -1,    -1,    -1,   108,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,  1224,    -1,    -1,    -1,    -1,    -1,   524,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     143,   516,    -1,  1246,  1247,   520,   149,    -1,    -1,   152,
      -1,    -1,    -1,    -1,    -1,   158,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   166,    -1,   503,   504,   505,   506,   507,
     508,   509,   510,   511,   512,   513,   514,   180,   516,    -1,
      -1,    -1,   520,    -1,    -1,   188,   189,    -1,    -1,    -1,
      -1,   194,    -1,   196,   197,    -1,   199,    -1,  1301,    -1,
      -1,    -1,   205,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     213,    -1,    -1,    -1,    -1,    -1,   219,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   229,   230,    -1,   232,
      -1,   505,   506,   507,   508,   509,   510,   511,   512,   513,
     514,   244,   516,    -1,    -1,    -1,    -1,   250,    -1,    -1,
      -1,   254,    -1,    -1,    -1,    -1,    -1,    -1,   261,    -1,
     263,   505,   506,   507,   508,   509,   510,   511,   512,   513,
     514,     3,   516,    -1,    -1,     7,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    27,    28,    29,    30,    31,
      32,   304,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    44,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   326,    -1,    -1,    -1,    -1,   331,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   339,    -1,   341,   342,
      -1,   344,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   352,
     353,   354,    -1,    -1,    -1,   358,    -1,    -1,    -1,   362,
      -1,    -1,    -1,   366,    -1,    -1,   369,    -1,    -1,    -1,
      -1,  1474,   375,   376,   377,    -1,    -1,    -1,  1481,    -1,
     383,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   394,   395,    -1,   397,    -1,   399,    -1,    -1,    -1,
     403,   404,    -1,    -1,    -1,    -1,    -1,   410,    -1,    -1,
     413,   414,   415,    -1,    -1,    -1,    -1,    -1,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,   436,   437,   438,   439,   440,   441,   524,
      -1,    -1,    -1,    -1,   447,   448,   449,    -1,   451,   452,
     453,   454,   455,    -1,    -1,    -1,    -1,    -1,    -1,   462,
     463,    -1,    -1,    -1,   196,   503,   504,   505,   506,   507,
     508,   509,   510,   511,   512,   513,   514,   515,   516,    -1,
      -1,   213,    -1,    -1,    -1,    -1,   524,    -1,    -1,   492,
     493,    -1,    -1,    -1,    -1,   498,   499,   500,   501,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   512,
      -1,    -1,    -1,    -1,    -1,   518,   519,     3,    -1,    -1,
     523,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    -1,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,    -1,
       3,    -1,    -1,    -1,     7,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      -1,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    44,    -1,    -1,    -1,    91,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   100,   503,   504,   505,   506,   507,
     508,   509,   510,   511,   512,   513,   514,   515,   516,    -1,
      -1,   117,    -1,    -1,    -1,    -1,   524,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    91,    -1,
      -1,    -1,    -1,    -1,     3,    -1,    -1,   100,     7,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    -1,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,   182,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   438,   439,   440,   441,
     196,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     452,   453,   454,   455,    -1,    -1,    -1,   213,    -1,    -1,
      -1,    -1,    -1,   465,   466,    -1,    -1,    -1,    -1,   182,
      -1,   184,    -1,    -1,    -1,   231,    -1,    -1,    -1,    -1,
      -1,   100,    -1,   196,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   249,   496,   497,    -1,    -1,   117,   255,
     213,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     512,    -1,    -1,    -1,    -1,    -1,   518,    -1,   231,    -1,
      -1,   523,    -1,    -1,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   249,   516,    -1,    -1,
     296,    -1,   255,   299,    -1,   524,   109,   110,   111,   112,
      -1,   114,   115,   116,    -1,   174,   119,   120,   121,   122,
     123,   124,    -1,   503,   504,   505,   506,   507,   508,   509,
     510,   511,   512,   513,   514,   515,   516,   196,    -1,    -1,
      -1,    -1,    -1,   296,   524,    -1,   299,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   213,    -1,    -1,    -1,    -1,    -1,
      -1,   357,    -1,    -1,   167,   503,   504,   505,   506,   507,
     508,   509,   510,   511,   512,   513,   514,   515,   516,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   524,   246,    -1,   192,
     249,    -1,    -1,    -1,    -1,    -1,   255,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   357,   208,   402,   503,   504,   505,
     506,   507,   508,   509,   510,   511,   512,   513,   514,   515,
     516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   389,   296,    -1,    -1,
     299,    -1,   438,   439,   440,   441,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   452,   453,   454,   455,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,   275,   276,   277,   278,   279,   280,   281,   282,
     283,   284,    -1,    -1,    -1,   438,   439,   440,   441,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   452,
     453,   454,   455,    -1,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   512,   516,    -1,    -1,
      -1,    -1,   518,    -1,   327,   524,    -1,   523,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   512,
      -1,    -1,    -1,    -1,    -1,   518,    -1,    -1,    -1,    -1,
     523,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   438,
     439,   440,   441,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   452,   453,   454,   455,     3,    -1,    -1,
      -1,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    -1,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,   503,
     504,   505,   506,   507,   508,   509,   510,   511,   512,   513,
     514,   515,   516,   512,    -1,    -1,    -1,    63,    -1,   518,
     524,    -1,    -1,    -1,   523,   503,   504,   505,   506,   507,
     508,   509,   510,   511,   512,   513,   514,   515,   516,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   524,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   100,    -1,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,     3,   516,
      -1,   117,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    -1,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
       3,    -1,    -1,    -1,     7,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
     196,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    44,    -1,    -1,    -1,    -1,    -1,   213,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     246,    -1,    -1,   249,    -1,    -1,    -1,    -1,    -1,   255,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   100,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     296,    -1,    -1,   299,   189,    -1,    -1,    -1,    -1,    -1,
      -1,   196,    -1,    -1,    -1,    -1,    -1,    -1,   203,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   213,    -1,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,    -1,    -1,   234,
     235,   524,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   196,    -1,    -1,   199,    -1,    -1,    -1,
     255,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     213,   503,   504,   505,   506,   507,   508,   509,   510,   511,
     512,   513,   514,   515,   516,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   524,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   296,    -1,   246,   299,    -1,   249,   503,   504,   505,
     506,   507,   508,   509,   510,   511,   512,   513,   514,   515,
     516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,    -1,
      -1,    -1,   438,   439,   440,   441,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   452,   453,   454,   455,
      -1,    -1,    -1,   296,    -1,    -1,   299,   503,   504,   505,
     506,   507,   508,   509,   510,   511,   512,   513,   514,   515,
     516,    -1,    -1,    -1,    -1,    -1,    -1,   372,   524,   503,
     504,   505,   506,   507,   508,   509,   510,   511,   512,   513,
     514,   515,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     524,    -1,    -1,    -1,    -1,    -1,   512,    -1,    -1,    -1,
      -1,    -1,   518,    -1,    -1,    -1,    -1,   523,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   372,
     504,   505,   506,   507,   508,   509,   510,   511,   512,   513,
     514,   515,   516,   438,   439,   440,   441,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   452,   453,   454,
     455,   503,   504,   505,   506,   507,   508,   509,   510,   511,
     512,   513,   514,   515,   516,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   524,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   438,   439,   440,   441,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   452,
     453,   454,   455,    -1,    -1,    -1,    -1,   512,    -1,    -1,
      -1,    -1,    -1,   518,    -1,     3,    -1,    -1,   523,     7,
       8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    30,    31,    32,    -1,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    -1,    -1,   512,
      -1,    -1,    -1,    -1,    -1,   518,    -1,    -1,    -1,     3,
     523,    -1,    -1,     7,     8,     9,    10,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    -1,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    -1,   100,    -1,    -1,    -1,     3,    -1,    -1,    -1,
       7,     8,     9,    10,    11,    12,    13,    14,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    -1,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    43,    44,     3,    -1,
      -1,    -1,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    -1,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   196,    -1,
      -1,   199,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   213,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   246,    -1,
      -1,   249,   196,   503,   504,   505,   506,   507,   508,   509,
     510,   511,   512,   513,   514,   515,   516,    -1,    -1,   213,
      -1,    -1,    -1,    -1,   524,   503,   504,   505,   506,   507,
     508,   509,   510,   511,   512,   513,   514,   515,   516,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   524,    -1,   296,   196,
      -1,   299,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,   213,    -1,    -1,    -1,
      -1,    -1,    -1,   524,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   196,    -1,    -1,    -1,   524,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   213,   516,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,    -1,    -1,
      -1,     3,    -1,    -1,   372,     7,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    27,    28,    29,    30,    31,
      32,    -1,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    44,   503,   504,   505,   506,   507,   508,   509,
     510,   511,   512,   513,   514,   515,   516,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   524,    -1,    -1,    -1,    -1,    -1,
     438,   439,   440,   441,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   452,   453,   454,   455,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,   438,   439,   440,   441,    -1,    -1,
      -1,   524,    -1,    -1,    -1,    -1,    -1,    -1,   452,   453,
     454,   455,    -1,    -1,   512,    -1,    -1,    -1,    -1,    -1,
     518,   465,   466,    -1,    -1,   523,    -1,    -1,    -1,   426,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   438,   439,   440,   441,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   496,   497,    -1,   452,   453,   454,   455,    -1,
      -1,    -1,    -1,    -1,   196,    -1,    -1,    -1,   512,    -1,
      -1,   426,    -1,    -1,   518,    -1,    -1,    -1,    -1,   523,
      -1,   213,    -1,   438,   439,   440,   441,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   452,   453,   454,
     455,   503,   504,   505,   506,   507,   508,   509,   510,   511,
     512,   513,   514,   515,   516,   512,    -1,    -1,    -1,    -1,
      -1,   518,   524,    -1,    -1,    -1,   523,   503,   504,   505,
     506,   507,   508,   509,   510,   511,   512,   513,   514,   515,
     516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   512,    -1,    -1,
      -1,    -1,    -1,   518,    -1,     3,    -1,    -1,   523,     7,
       8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    30,    31,    32,    -1,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,     3,    -1,    -1,
      -1,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    -1,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,   503,
     504,   505,   506,   507,   508,   509,   510,   511,   512,   513,
     514,   515,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     524,   503,   504,   505,   506,   507,   508,   509,   510,   511,
     512,   513,   514,   515,   516,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   524,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   438,   439,   440,   441,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     452,   453,   454,   455,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   524,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   196,    -1,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,   213,    -1,    -1,    -1,    -1,
     512,   524,    -1,    -1,    -1,    -1,   518,    -1,    -1,    -1,
      -1,   523,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     196,   503,   504,   505,   506,   507,   508,   509,   510,   511,
     512,   513,   514,   515,   516,    -1,    -1,   213,    -1,    -1,
      -1,    -1,   524,   503,   504,   505,   506,   507,   508,   509,
     510,   511,   512,   513,   514,   515,   516,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   524,   503,   504,   505,   506,   507,
     508,   509,   510,   511,   512,   513,   514,   515,   516,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   524,   503,   504,   505,
     506,   507,   508,   509,   510,   511,   512,   513,   514,   515,
     516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,   503,
     504,   505,   506,   507,   508,   509,   510,   511,   512,   513,
     514,   515,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     524,   503,   504,   505,   506,   507,   508,   509,   510,   511,
     512,   513,   514,   515,   516,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   524,   503,   504,   505,   506,   507,   508,   509,
     510,   511,   512,   513,   514,   515,   516,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   524,   503,   504,   505,   506,   507,
     508,   509,   510,   511,   512,   513,   514,   515,   516,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   524,   503,   504,   505,
     506,   507,   508,   509,   510,   511,   512,   513,   514,   515,
     516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,    -1,
     438,   439,   440,   441,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   452,   453,   454,   455,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,
      -1,    -1,   438,   439,   440,   441,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   452,   453,   454,   455,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,   512,    -1,    -1,    -1,    -1,    -1,
     518,   524,    -1,    -1,    -1,   523,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   512,    -1,    -1,    -1,
      -1,    -1,   518,    -1,    -1,    -1,    -1,   523,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   524,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   524,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   524,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,    -1,   516,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   524,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,    -1,   516,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   524,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,    -1,   516,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   524,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,    -1,   516,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
      -1,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,    -1,   516,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   524,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,    -1,   516,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   524,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,    -1,   516,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   524,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,    -1,   516,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
      -1,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,    -1,   516,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   524,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,    -1,   516,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   524,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,    -1,   516,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   524,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,    -1,   516,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
      -1,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,    -1,   516,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   524,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,    -1,   516,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   524,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,    -1,   516,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   524,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,    -1,   516,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
      -1,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,    -1,   516,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   524,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,    -1,   516,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   524,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,    -1,   516,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   524,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,    -1,   516,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,    -1,    -1,   522,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,    -1,    -1,   522,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,    -1,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,    -1,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
      -1,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,   519,   503,   504,   505,   506,   507,   508,   509,
     510,   511,   512,   513,   514,    -1,   516,    -1,    -1,   519,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint16 yystos[] =
{
       0,     1,     3,     7,     8,     9,    10,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      51,    56,    57,    68,    72,    77,    87,    92,    94,    97,
      99,   100,   101,   108,   143,   149,   152,   158,   166,   180,
     188,   189,   194,   196,   197,   199,   205,   213,   219,   229,
     230,   232,   244,   250,   254,   261,   263,   304,   326,   331,
     339,   341,   342,   344,   352,   353,   354,   358,   362,   366,
     369,   375,   376,   377,   383,   394,   395,   397,   399,   403,
     404,   410,   413,   414,   415,   436,   437,   438,   439,   440,
     441,   447,   448,   449,   451,   452,   453,   454,   455,   462,
     463,   492,   493,   498,   499,   500,   501,   512,   518,   523,
     527,   528,   529,   530,   533,   545,   547,   551,   553,   556,
     558,   560,   561,   562,   563,   564,   565,   566,   569,   570,
     571,   572,   594,   595,   596,   597,   519,   502,   521,   523,
     523,   523,   523,   523,   523,   523,   523,   523,   523,   523,
     523,   523,   523,   523,   523,     6,   523,   523,   523,   523,
     523,   523,   523,   523,   523,   523,   523,    58,   100,   126,
     137,   155,   158,   177,   190,   315,   360,   539,    91,   100,
     182,   199,   246,   249,   296,   299,   573,   588,   329,   358,
       3,     9,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    30,    31,    32,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,   100,   117,   174,   196,
     246,   249,   255,   512,   518,   523,   588,   597,     6,    87,
     244,   383,   397,    29,     8,     8,     8,   137,   552,     6,
     199,   264,   358,   420,   539,    33,    50,   137,   252,   255,
     290,   426,   546,    95,   181,   196,   366,   588,    75,   100,
     174,   246,   249,   426,   588,     8,    68,   176,   204,   234,
     426,   525,   588,   571,   572,     8,    78,    80,    86,   127,
     137,   154,   174,   291,   548,    66,   100,   137,   174,   249,
     252,   255,   328,   342,   409,   436,   448,   467,   557,   523,
     216,     8,   196,   597,     8,    33,   174,   239,   361,   366,
     389,   426,   554,   196,   197,    63,   100,   117,   196,   246,
     249,   255,   588,   597,   304,   588,   532,   100,   240,   252,
     255,   333,   363,   427,   568,   588,   140,   223,   238,   355,
     428,   597,    51,   143,   152,   199,   404,     8,    77,   531,
     152,     8,    66,   100,   137,   174,   249,   252,   255,   328,
     436,   448,   467,   567,   588,   136,   317,   580,     8,   448,
      91,   100,   117,   182,   196,   231,   249,   255,   357,   402,
     588,   597,     6,   100,   182,   249,   372,     6,   539,   100,
     174,   252,   255,   418,   568,   588,    91,   100,   182,   184,
     231,   249,   255,   357,   389,   588,   597,   100,   174,   240,
     252,   255,   328,   389,   394,   395,   427,   445,   568,   588,
       6,   100,   182,   249,   372,   465,   466,   496,   497,   597,
     100,   240,   252,   255,   363,   427,   568,   588,    51,    87,
     143,   152,   158,   195,   196,   199,   205,   244,   376,   397,
     404,   501,   465,   466,   496,   497,   597,   588,   100,   240,
     252,   255,   363,   427,   568,   588,   573,   573,   597,    66,
     100,   137,   174,   249,   252,   255,   328,   363,   436,   448,
     467,   500,   559,   588,     3,    29,   596,   597,   596,   597,
     596,   597,     0,   519,   519,   519,   519,   519,   519,   519,
     519,   519,   519,   519,   519,   519,   519,   519,   519,   519,
     519,   519,    76,   234,   407,   409,   574,   588,   519,   519,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   516,   519,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,   519,   596,
     597,   597,   596,   597,   596,   597,   596,   597,   596,   597,
     596,   597,   596,   597,   596,   597,   596,   597,   596,   597,
     596,   597,   596,   597,   596,   597,   596,   597,     8,   596,
     597,   596,   597,   596,   597,   596,   597,   596,   597,   596,
     597,   596,   597,   596,   597,   596,   597,   596,   597,   596,
     597,   596,   597,   596,   597,   542,     8,   360,    89,   155,
     177,   226,   291,   296,   299,   587,    31,   540,     8,   190,
     260,     8,   543,   177,   372,     8,     8,   588,     8,     8,
       6,     8,   521,   523,   523,   523,   523,   523,   523,   523,
     523,   523,   523,   523,   523,   523,   523,   523,   523,   523,
     523,   523,   523,   523,   523,   523,   523,   523,   523,   523,
       8,   100,   292,   310,   591,     8,     8,   436,   448,   589,
     597,   597,   597,   511,   512,   513,   514,   516,   520,   520,
     520,   520,   588,   568,     8,    66,   100,   137,   155,   174,
     177,   249,   252,   255,   261,   328,   359,   362,   430,   436,
     448,   467,   538,     8,     8,   137,   155,   177,   195,   244,
     266,   290,   345,   359,   360,   389,   426,   434,   534,   597,
     588,   588,   597,   589,     8,   588,   568,   100,     8,   588,
       8,     8,     8,   169,   425,   426,   590,   184,   328,   426,
     467,   588,   590,   588,   256,   257,   258,   467,   475,   582,
     438,   439,   440,   441,   452,   453,   454,   455,   574,   137,
     550,   137,   137,   549,   137,   588,   137,   588,   137,   568,
       8,     8,    51,   207,   264,   588,   100,     8,   597,   589,
       8,     8,     8,   597,   597,   597,   568,   597,     6,   100,
     246,   249,   520,   597,   426,   189,   203,   234,   235,   255,
     372,   555,   588,   597,   173,   229,   296,   299,   426,   597,
       8,   568,     8,   372,   426,     8,     8,     8,   589,   520,
     588,   553,     8,     8,   597,   589,    63,    91,   100,   174,
     182,   184,   246,   249,   328,   402,     8,    98,   233,   259,
     274,   593,     8,   155,     6,     8,    85,   177,   264,   359,
     400,   420,     8,   520,     6,   553,   420,     8,     8,   264,
     588,   100,     8,   597,   589,     8,   597,   597,   597,   568,
     597,   372,     8,     6,     8,     8,     8,   589,     8,     8,
     100,   174,   255,   372,   520,     8,     8,     8,     8,     8,
     100,   597,   589,   412,   372,     8,     8,     6,     8,     8,
     589,     8,     6,   520,     8,   100,     8,   597,   589,     8,
       8,     8,     8,     8,     8,     8,     8,     8,     8,     8,
       8,     8,     8,   520,     8,     8,   597,   589,     8,   593,
       8,     8,     8,     8,   264,   420,     8,     8,     8,     8,
       8,     8,     8,     8,     8,   597,   597,   597,   597,   520,
       8,     8,   597,   589,     8,   593,   520,     8,     8,   264,
     588,   100,     8,   597,   589,     8,     8,   597,   597,   597,
       8,   568,   524,   524,   100,   246,   249,   579,   588,     6,
      91,   100,   182,   231,   235,   249,   318,   578,     8,    59,
      61,    91,   100,   182,   184,   203,   231,   235,   249,   300,
     328,   330,   371,   373,   387,   389,   395,   426,   435,   576,
     577,   588,     8,    29,    31,    32,    54,    84,   100,   118,
     218,   249,   262,   271,   297,   298,   300,   303,   372,   378,
     426,   575,   583,   588,   596,   597,   596,   596,   596,   596,
     596,   596,   596,   596,   597,   596,   597,   596,   597,   596,
     597,   596,   597,   597,   597,   597,   597,   597,   597,   597,
     597,   596,   597,   596,   597,   596,   597,   596,   597,   597,
     596,   597,   522,   524,   524,   524,   524,   524,   524,   524,
     524,   520,   520,   524,   524,   524,   524,   524,   524,   524,
     524,   524,   524,   524,   524,   520,   520,   524,   524,   524,
     524,   524,   524,   524,   524,   524,   520,   520,   520,   520,
     520,   520,   520,   520,   524,   524,   524,   524,   524,   524,
     524,   524,   524,   524,   524,   524,   553,   177,    31,    58,
     127,   265,   315,   429,   588,   260,   446,   588,   541,   553,
       8,   544,   553,    32,   153,   290,   342,   426,   467,   588,
       8,   520,   597,   597,   597,   597,   597,   597,   597,   597,
     597,   597,   597,   597,   597,   597,     8,   597,   597,   597,
     597,   597,   597,   597,   597,   597,   597,   597,   597,   597,
       8,     8,   597,   597,   597,   597,   597,   597,     8,     8,
       8,     6,     8,     8,   264,   588,     6,   100,     6,     8,
     597,   589,     6,   537,     8,     6,   536,     6,   535,     6,
     597,   597,   597,   568,   520,   195,   244,   260,   588,     6,
       6,     8,     8,     8,     8,     8,   173,     6,     8,     8,
       8,   177,   553,   520,   520,     8,    67,   378,   109,   110,
     111,   112,   114,   115,   116,   119,   120,   121,   122,   123,
     124,   167,   192,   208,   275,   276,   277,   278,   279,   280,
     281,   282,   283,   284,   327,   585,     8,     8,   597,   588,
     553,   429,   588,   588,   568,   588,   588,   588,     8,   100,
     590,   588,     8,   520,   520,     6,   520,   520,   520,   520,
     524,     8,     8,     8,   597,   520,     8,   100,   588,     8,
     588,   209,   434,   586,   589,     8,   568,   520,     6,     8,
       8,   394,     8,     8,   597,   520,     8,   372,   426,   372,
       8,   100,   310,   588,   591,     8,   585,     8,     8,     8,
       8,   372,     8,   389,     6,   345,     8,   420,   290,   467,
     420,     8,    61,     6,     8,     8,   588,     8,   520,   520,
     520,   520,   520,   520,     8,     8,     8,    47,    81,    84,
      85,   238,   355,   417,   592,     8,   597,     8,   520,     8,
       8,   597,     8,   520,   597,   520,     8,     8,   597,   520,
     597,   588,     8,   520,   520,   520,   520,   520,     8,     8,
       8,   372,     8,     8,    89,   238,   355,   581,   305,   313,
       8,     8,    67,   378,   520,     8,     6,   372,     8,     8,
       8,   585,     8,   581,   209,   378,   434,     8,   592,     8,
       6,    48,   288,   293,   584,     8,     8,   426,   597,   426,
     597,    67,   378,     8,   577,   520,   588,   597,   597,   588,
       8,     8,     8,   100,   199,   246,   249,   372,   588,   597,
     100,   199,   246,   249,   372,   588,   597,   597,   597,   592,
       8,     8,    67,   378,   502,   596,   597,   596,   597,   596,
     597,   597,   597,   596,   597,   596,   597,   596,   597,   153,
     290,   467,   588,   588,   588,   588,   588,   588,   588,   588,
     553,   553,     8,     8,     8,     6,   467,   290,   467,   597,
       8,   522,   524,   520,   524,   520,   520,   520,   520,   588,
       8,   520,   553,   520,   553,   553,   520,   520,   520,   597,
       8,   588,   588,   520,   520,     6,   395,     6,   597,   597,
     585,   520,   520,   588,   133,   100,   590,     8,   597,     8,
       8,   597,   597,   597,   597,   597,     8,   520,   520,   597,
     173,   597,   520,   597,     8,     8,     8,     8,     8,     8,
     199,     8,   389,     6,     6,    61,     8,     8,    61,     6,
     597,     8,   597,   597,   597,   597,   597,   597,   520,   597,
     520,   597,   520,   597,     8,   597,   597,   597,     8,   520,
       6,     8,    67,   378,    67,   378,   597,     8,   588,     8,
       8,     8,     8,   588,     8,     8,     8,   597,   524,   524,
     524,   524,   524,   524,   520,   520,   524,   524,   524,   524,
     524,   524,     8,     8,   597,     6,   520,   373,   597,   597,
       8,     8,     8,   597,     8,   597,   597,   597,   520,   597,
       8,     8,     8,   597,   133,     8,     8,   520,   162,     8,
       8,     6,   597,     6,   395,     8,   389,   389,   520,   389,
     520,   520,   597,   597,   597,   520,     8,   597,   597,   520,
     597,     8,   524,   520,   520,   520,   520,   520,   597,   520,
     373,     8,   597,   597,   520,     6,     8,   395,     8,     8,
       8,     8,   597,   597,   520,   520,   520,   597,   524,   524,
     597,   597,     8,     8,     8,   597,   597,     8,   520,   597,
     373,     8,   395,   394,   389,   395,   520,   520,   597,   597,
     597,   520,   520,   520,   597,     8,   373,     8,     8,     8,
       8,   597,   597,   597,   597,   597,     8,   373,   373,   394,
     373,   409,   520,   127,   241,     8,     8,     8,     8,   597,
     597,     8,     8,    61,    61,   373,    61,   127,   520,    61,
       8,     8,   597,    61,   520,   597,   520,   597
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint16 yyr1[] =
{
       0,   526,   527,   527,   527,   527,   527,   527,   527,   527,
     527,   527,   527,   527,   527,   527,   527,   527,   527,   527,
     527,   527,   527,   527,   527,   527,   527,   528,   528,   528,
     528,   528,   528,   528,   528,   528,   528,   528,   528,   528,
     528,   528,   528,   528,   528,   528,   528,   528,   528,   528,
     528,   528,   528,   528,   528,   528,   528,   528,   528,   528,
     528,   528,   528,   528,   528,   528,   528,   528,   528,   528,
     528,   528,   528,   528,   528,   529,   529,   529,   529,   529,
     529,   529,   529,   529,   529,   529,   529,   530,   530,   530,
     530,   530,   530,   530,   530,   530,   530,   530,   530,   530,
     530,   530,   530,   530,   530,   530,   530,   530,   530,   530,
     530,   530,   531,   530,   532,   530,   533,   533,   533,   533,
     533,   533,   533,   533,   533,   533,   533,   533,   533,   533,
     533,   533,   533,   533,   533,   533,   533,   534,   533,   533,
     533,   533,   533,   533,   533,   533,   533,   533,   533,   533,
     533,   533,   533,   533,   533,   533,   533,   533,   533,   533,
     533,   533,   533,   533,   533,   533,   533,   533,   533,   533,
     533,   533,   533,   533,   535,   533,   536,   533,   537,   533,
     538,   533,   533,   533,   533,   533,   533,   533,   533,   533,
     533,   533,   533,   533,   539,   539,   539,   539,   539,   539,
     539,   539,   539,   539,   539,   539,   540,   539,   541,   539,
     542,   539,   543,   539,   544,   539,   539,   539,   539,   539,
     539,   539,   539,   539,   539,   539,   539,   539,   539,   539,
     539,   545,   545,   545,   545,   546,   545,   545,   545,   545,
     545,   547,   547,   547,   547,   547,   547,   547,   548,   547,
     549,   547,   547,   547,   547,   547,   550,   547,   547,   547,
     547,   551,   551,   552,   551,   551,   551,   553,   553,   553,
     553,   553,   553,   553,   553,   553,   553,   553,   553,   553,
     553,   553,   553,   553,   553,   553,   553,   553,   553,   554,
     553,   555,   553,   553,   553,   556,   557,   556,   556,   556,
     556,   556,   556,   556,   556,   556,   556,   556,   556,   556,
     556,   556,   556,   556,   556,   556,   558,   559,   558,   558,
     558,   558,   558,   558,   558,   558,   558,   558,   558,   558,
     558,   558,   558,   558,   560,   560,   560,   560,   560,   560,
     560,   560,   561,   561,   561,   561,   561,   561,   561,   561,
     562,   562,   562,   562,   562,   562,   562,   562,   563,   563,
     563,   563,   563,   563,   563,   564,   564,   564,   564,   564,
     564,   564,   564,   564,   564,   564,   565,   565,   565,   565,
     565,   565,   565,   565,   565,   565,   565,   565,   565,   566,
     567,   566,   566,   566,   566,   566,   566,   566,   566,   566,
     566,   566,   566,   566,   566,   568,   568,   568,   568,   568,
     568,   568,   568,   568,   568,   568,   568,   568,   568,   568,
     568,   569,   569,   569,   569,   569,   569,   569,   569,   569,
     569,   569,   569,   569,   569,   569,   569,   569,   569,   569,
     569,   569,   569,   569,   569,   569,   569,   569,   569,   569,
     569,   569,   569,   569,   569,   569,   569,   569,   569,   569,
     569,   569,   569,   569,   569,   569,   569,   569,   569,   569,
     569,   569,   569,   569,   569,   569,   569,   570,   570,   570,
     570,   570,   571,   571,   571,   571,   571,   571,   572,   572,
     572,   573,   573,   573,   573,   573,   573,   573,   574,   574,
     574,   574,   574,   575,   575,   575,   575,   575,   575,   575,
     575,   575,   575,   575,   575,   575,   575,   575,   575,   575,
     575,   575,   575,   575,   575,   575,   575,   575,   575,   575,
     575,   575,   575,   575,   576,   576,   577,   577,   577,   577,
     577,   577,   577,   577,   577,   577,   577,   577,   577,   577,
     577,   577,   577,   577,   577,   577,   577,   577,   577,   577,
     577,   577,   577,   577,   577,   577,   578,   578,   578,   578,
     578,   578,   578,   578,   578,   578,   578,   579,   579,   579,
     579,   580,   580,   581,   581,   581,   582,   582,   582,   582,
     582,   583,   583,   583,   584,   584,   584,   585,   585,   585,
     585,   585,   585,   585,   585,   585,   585,   585,   585,   585,
     585,   585,   585,   585,   585,   585,   585,   585,   585,   585,
     585,   585,   585,   585,   586,   586,   587,   587,   587,   587,
     588,   588,   589,   589,   590,   590,   591,   591,   591,   592,
     592,   592,   592,   592,   592,   592,   593,   593,   593,   593,
     594,   595,   595,   596,   596,   596,   596,   596,   596,   596,
     596,   596,   596,   596,   596,   596,   596,   596,   596,   596,
     596,   596,   596,   596,   596,   596,   596,   596,   596,   596,
     596,   596,   596,   596,   596,   596,   596,   596,   596,   596,
     596,   596,   596,   596,   596,   596,   596,   596,   596,   596,
     596,   596,   596,   596,   596,   596,   596,   596,   596,   596,
     596,   596,   596,   596,   596,   596,   596,   596,   596,   596,
     596,   596,   596,   596,   596,   597,   597,   597,   597,   597,
     597,   597,   597,   597,   597,   597,   597,   597,   597,   597,
     597,   597,   597,   597,   597,   597,   597,   597,   597,   597,
     597,   597,   597,   597,   597,   597,   597,   597,   597,   597,
     597,   597,   597,   597,   597,   597,   597,   597,   597,   597,
     597,   597,   597,   597,   597,   597,   597,   597,   597,   597,
     597,   597,   597,   597,   597,   597,   597,   597,   597,   597,
     597,   597,   597,   597,   597,   597,   597,   597
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     0,     2,     2,
       3,     2,     2,     8,     3,     3,     3,     3,     3,     4,
       4,     2,     2,     3,     2,     2,     2,     8,     3,     3,
       3,     3,     3,     4,     4,     2,     2,     2,     3,     2,
       2,     4,     3,     3,     3,     3,     3,     3,     3,     4,
       4,     4,     4,     4,     3,     1,     1,     1,     1,     2,
       2,     1,     4,     5,     7,     3,     1,     2,     3,     2,
       1,     1,     2,     1,     8,     1,     1,     2,     2,     2,
       2,     2,     3,     3,     2,     2,     3,     1,     2,     8,
       8,     8,     0,     3,     0,     3,     3,     2,     2,     3,
       2,     7,     3,     2,     3,     4,     5,    12,    13,    12,
      13,    11,    12,    11,    12,     4,     4,     0,     4,     4,
       5,     4,     4,     4,     4,     5,     5,     4,     5,     4,
       8,     4,     4,     7,     6,     5,    12,     3,     4,     4,
       5,     4,     5,     4,     4,     4,     4,     4,     4,     4,
      11,    12,    13,    14,     0,     5,     0,     5,     0,     5,
       0,     4,     6,     4,     4,     5,     4,     4,     4,     5,
      10,     6,     6,     6,     2,     3,     4,     4,     4,     4,
       4,     4,     4,     4,     3,     2,     0,     3,     0,     4,
       0,     3,     0,     3,     0,     4,     3,     2,     2,     3,
       4,     4,     5,     4,     4,     4,     4,     6,     5,     5,
       7,     3,     3,     3,     3,     0,     3,     3,     3,     5,
       5,     2,     3,     3,     4,     6,     7,     4,     0,     3,
       0,     4,     4,     4,     4,     3,     0,     4,     2,     1,
       4,     3,     3,     0,     3,     3,     9,     3,     4,     4,
       4,     5,     2,     3,     4,     4,     3,     3,     4,     6,
       4,     5,     5,     4,     4,     5,     4,     4,     4,     0,
       3,     0,     4,     6,     6,     3,     0,     3,     5,     3,
       5,     3,     4,     3,     5,     6,     3,     3,     4,     4,
       4,     5,     9,     5,     5,     5,     3,     0,     3,     5,
       3,     3,     4,     2,     3,     3,     3,     3,     3,     4,
       9,     5,     5,     5,     2,     3,     3,     3,     3,     5,
       3,     2,     2,     3,     3,     3,     3,     5,     3,     2,
       2,     3,     3,     3,     3,     5,     3,     2,     2,     3,
       4,     4,     3,     5,     2,     2,     3,     4,     3,     3,
       3,     3,     3,     3,     4,     3,     2,     3,     3,     3,
       3,     3,     3,     3,     3,     4,     3,     5,     2,     3,
       0,     3,     5,     3,     3,     4,     2,     3,     3,     3,
       4,     9,     5,     5,     5,     3,     3,     3,     3,     3,
       3,     4,     3,     4,     3,     3,     4,     4,     3,     4,
       4,     2,     3,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     1,     1,     1,     2,    17,     2,     8,     3,
       3,     3,     3,     8,     3,     3,     3,     3,     2,     3,
       3,     3,     3,     2,     3,     3,     3,     3,     2,     3,
       3,     3,     3,     3,     4,     2,     3,     4,     4,     3,
       3,     3,     3,     5,     6,     6,     4,     2,     1,     2,
       3,     2,     1,     1,     1,     1,     1,     1,     2,     2,
       2,     1,     2,     2,     2,     2,     3,     2,     2,     2,
       2,     2,     1,     1,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     1,     2,     2,     3,     3,     2,
       2,     3,     3,     3,     3,     3,     3,     3,     3,     2,
       2,     2,     2,     3,     1,     2,     1,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     3,     3,     3,     3,     2,
       2,     3,     2,     2,     2,     3,     1,     2,     2,     2,
       2,     4,     2,     3,     2,     2,     2,     1,     2,     2,
       2,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       6,     3,     3,     1,     1,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     4,     4,     4,     4,     4,
       6,     4,     4,     1,     1,     1,     4,     4,     4,     4,
       6,     6,     6,     6,     1,     1,     4,     4,     4,     4,
       4,     8,     6,     6,     6,     4,     4,     1,     1,     1,
       4,     4,     4,     4,     3,     3,     3,     3,     3,     3,
       3,     3,     2,     3,     2,     1,     1,     4,     3,     3,
       3,     3,     3,     3,     4,     4,     4,     4,     6,     4,
       4,     1,     1,     1,     4,     4,     4,     4,     6,     3,
       3,     3,     3,     3,     3,     3,     3,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     4,     4,     4,
       4,     4,     8,     6,     6,     6,     4,     4,     1,     1,
       1,     4,     4,     4,     4,     5,     7,     3,     3,     3,
       3,     3,     3,     3,     3,     2,     3,     2
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;                                                  \
    }                                                           \
while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256



/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)

/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule)
{
  unsigned long int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &(yyvsp[(yyi + 1) - (yynrhs)])
                                              );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
yystrlen (const char *yystr)
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
yystpcpy (char *yydest, const char *yysrc)
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
        switch (*++yyp)
          {
          case '\'':
          case ',':
            goto do_not_strip_quotes;

          case '\\':
            if (*++yyp != '\\')
              goto do_not_strip_quotes;
            /* Fall through.  */
          default:
            if (yyres)
              yyres[yyn] = *yyp;
            yyn++;
            break;

          case '"':
            if (yyres)
              yyres[yyn] = '\0';
            return yyn;
          }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
{
  YYUSE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *yyvs1 = yyvs;
        yytype_int16 *yyss1 = yyss;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * sizeof (*yyssp),
                    &yyvs1, yysize * sizeof (*yyvsp),
                    &yystacksize);

        yyss = yyss1;
        yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yytype_int16 *yyss1 = yyss;
        union yyalloc *yyptr =
          (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex ();
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 3:
#line 669 "gram.y" /* yacc.c:1646  */
    {}
#line 4684 "y.tab.c" /* yacc.c:1646  */
    break;

  case 4:
#line 670 "gram.y" /* yacc.c:1646  */
    {}
#line 4690 "y.tab.c" /* yacc.c:1646  */
    break;

  case 5:
#line 671 "gram.y" /* yacc.c:1646  */
    { result = (yyvsp[-1].val); }
#line 4696 "y.tab.c" /* yacc.c:1646  */
    break;

  case 6:
#line 672 "gram.y" /* yacc.c:1646  */
    { result = *(yyvsp[-1].ptr); }
#line 4702 "y.tab.c" /* yacc.c:1646  */
    break;

  case 18:
#line 684 "gram.y" /* yacc.c:1646  */
    { }
#line 4708 "y.tab.c" /* yacc.c:1646  */
    break;

  case 26:
#line 692 "gram.y" /* yacc.c:1646  */
    { return 1; }
#line 4714 "y.tab.c" /* yacc.c:1646  */
    break;

  case 28:
#line 696 "gram.y" /* yacc.c:1646  */
    { do_clear_boxes(); }
#line 4720 "y.tab.c" /* yacc.c:1646  */
    break;

  case 29:
#line 697 "gram.y" /* yacc.c:1646  */
    { curbox = next_box(); }
#line 4726 "y.tab.c" /* yacc.c:1646  */
    break;

  case 30:
#line 698 "gram.y" /* yacc.c:1646  */
    { curbox = (int) (yyvsp[0].val); }
#line 4732 "y.tab.c" /* yacc.c:1646  */
    break;

  case 31:
#line 699 "gram.y" /* yacc.c:1646  */
    { boxes[curbox].active = (yyvsp[0].pset); }
#line 4738 "y.tab.c" /* yacc.c:1646  */
    break;

  case 32:
#line 700 "gram.y" /* yacc.c:1646  */
    { boxes[curbox].gno = (yyvsp[0].pset); }
#line 4744 "y.tab.c" /* yacc.c:1646  */
    break;

  case 33:
#line 702 "gram.y" /* yacc.c:1646  */
    {
	    if (curbox >= 0 && curbox < maxboxes) {
		boxes[curbox].x1 = (yyvsp[-6].val);
		boxes[curbox].y1 = (yyvsp[-4].val);
		boxes[curbox].x2 = (yyvsp[-2].val);
		boxes[curbox].y2 = (yyvsp[0].val);
	    }
	}
#line 4757 "y.tab.c" /* yacc.c:1646  */
    break;

  case 34:
#line 710 "gram.y" /* yacc.c:1646  */
    { sysbox.loctype = (yyvsp[0].pset); }
#line 4763 "y.tab.c" /* yacc.c:1646  */
    break;

  case 35:
#line 711 "gram.y" /* yacc.c:1646  */
    { sysbox.lines = (int) (yyvsp[0].val); }
#line 4769 "y.tab.c" /* yacc.c:1646  */
    break;

  case 36:
#line 712 "gram.y" /* yacc.c:1646  */
    { sysbox.linew = (int) (yyvsp[0].val); }
#line 4775 "y.tab.c" /* yacc.c:1646  */
    break;

  case 37:
#line 713 "gram.y" /* yacc.c:1646  */
    { sysbox.color = (int) (yyvsp[0].val); }
#line 4781 "y.tab.c" /* yacc.c:1646  */
    break;

  case 38:
#line 714 "gram.y" /* yacc.c:1646  */
    { sysbox.fill = (yyvsp[0].pset); }
#line 4787 "y.tab.c" /* yacc.c:1646  */
    break;

  case 39:
#line 715 "gram.y" /* yacc.c:1646  */
    { sysbox.fillcolor = (int) (yyvsp[0].val); }
#line 4793 "y.tab.c" /* yacc.c:1646  */
    break;

  case 40:
#line 716 "gram.y" /* yacc.c:1646  */
    { sysbox.fillpattern = (int) (yyvsp[0].val); }
#line 4799 "y.tab.c" /* yacc.c:1646  */
    break;

  case 41:
#line 718 "gram.y" /* yacc.c:1646  */
    {
	    if (curbox >= 0 && curbox < maxboxes) {
                    boxes[curbox].loctype = sysbox.loctype;
                    boxes[curbox].color = sysbox.color;
                    boxes[curbox].linew = sysbox.linew;
                    boxes[curbox].lines = sysbox.lines;
                    boxes[curbox].fill = sysbox.fill;
                    boxes[curbox].fillcolor = sysbox.fillcolor;
                    boxes[curbox].fillpattern = sysbox.fillpattern;
	    }
	}
#line 4815 "y.tab.c" /* yacc.c:1646  */
    break;

  case 42:
#line 729 "gram.y" /* yacc.c:1646  */
    { curline = next_line(); }
#line 4821 "y.tab.c" /* yacc.c:1646  */
    break;

  case 43:
#line 730 "gram.y" /* yacc.c:1646  */
    { curline = (int) (yyvsp[0].val); }
#line 4827 "y.tab.c" /* yacc.c:1646  */
    break;

  case 44:
#line 731 "gram.y" /* yacc.c:1646  */
    { do_clear_lines(); }
#line 4833 "y.tab.c" /* yacc.c:1646  */
    break;

  case 45:
#line 732 "gram.y" /* yacc.c:1646  */
    { lines[curline].active = (yyvsp[0].pset); }
#line 4839 "y.tab.c" /* yacc.c:1646  */
    break;

  case 46:
#line 733 "gram.y" /* yacc.c:1646  */
    { lines[curline].gno = (yyvsp[0].pset); }
#line 4845 "y.tab.c" /* yacc.c:1646  */
    break;

  case 47:
#line 735 "gram.y" /* yacc.c:1646  */
    {
	    lines[curline].x1 = (yyvsp[-6].val);
	    lines[curline].y1 = (yyvsp[-4].val);
	    lines[curline].x2 = (yyvsp[-2].val);
	    lines[curline].y2 = (yyvsp[0].val);
	}
#line 4856 "y.tab.c" /* yacc.c:1646  */
    break;

  case 48:
#line 741 "gram.y" /* yacc.c:1646  */
    { sysline.loctype = (yyvsp[0].pset); }
#line 4862 "y.tab.c" /* yacc.c:1646  */
    break;

  case 49:
#line 742 "gram.y" /* yacc.c:1646  */
    { sysline.linew = (int) (yyvsp[0].val); }
#line 4868 "y.tab.c" /* yacc.c:1646  */
    break;

  case 50:
#line 743 "gram.y" /* yacc.c:1646  */
    { sysline.lines = (int) (yyvsp[0].val); }
#line 4874 "y.tab.c" /* yacc.c:1646  */
    break;

  case 51:
#line 744 "gram.y" /* yacc.c:1646  */
    { sysline.color = (int) (yyvsp[0].val); }
#line 4880 "y.tab.c" /* yacc.c:1646  */
    break;

  case 52:
#line 745 "gram.y" /* yacc.c:1646  */
    { sysline.arrow = (int) (yyvsp[0].val); }
#line 4886 "y.tab.c" /* yacc.c:1646  */
    break;

  case 53:
#line 746 "gram.y" /* yacc.c:1646  */
    { sysline.asize = (yyvsp[0].val); }
#line 4892 "y.tab.c" /* yacc.c:1646  */
    break;

  case 54:
#line 747 "gram.y" /* yacc.c:1646  */
    { sysline.atype = (int) (yyvsp[0].val); }
#line 4898 "y.tab.c" /* yacc.c:1646  */
    break;

  case 55:
#line 749 "gram.y" /* yacc.c:1646  */
    {
	    if (curline >= 0 && curline < maxlines) {
		lines[curline].lines = sysline.lines;
		lines[curline].linew = sysline.linew;
		lines[curline].color = sysline.color;
		lines[curline].arrow = sysline.arrow;
		lines[curline].asize = sysline.asize;
		lines[curline].atype = sysline.atype;
		lines[curline].loctype = sysline.loctype;
	    }
	}
#line 4914 "y.tab.c" /* yacc.c:1646  */
    break;

  case 56:
#line 760 "gram.y" /* yacc.c:1646  */
    { do_clear_text(); }
#line 4920 "y.tab.c" /* yacc.c:1646  */
    break;

  case 57:
#line 761 "gram.y" /* yacc.c:1646  */
    { curstring = next_string(); }
#line 4926 "y.tab.c" /* yacc.c:1646  */
    break;

  case 58:
#line 762 "gram.y" /* yacc.c:1646  */
    { curstring = (int) (yyvsp[0].val); }
#line 4932 "y.tab.c" /* yacc.c:1646  */
    break;

  case 59:
#line 763 "gram.y" /* yacc.c:1646  */
    { pstr[curstring].active = (yyvsp[0].pset); }
#line 4938 "y.tab.c" /* yacc.c:1646  */
    break;

  case 60:
#line 764 "gram.y" /* yacc.c:1646  */
    { pstr[curstring].gno = (yyvsp[0].pset); }
#line 4944 "y.tab.c" /* yacc.c:1646  */
    break;

  case 61:
#line 766 "gram.y" /* yacc.c:1646  */
    {
	    pstr[curstring].x = (yyvsp[-2].val);
	    pstr[curstring].y = (yyvsp[0].val);
	}
#line 4953 "y.tab.c" /* yacc.c:1646  */
    break;

  case 62:
#line 770 "gram.y" /* yacc.c:1646  */
    { sysstr.loctype = (yyvsp[0].pset); }
#line 4959 "y.tab.c" /* yacc.c:1646  */
    break;

  case 63:
#line 771 "gram.y" /* yacc.c:1646  */
    { sysstr.linew = (int) (yyvsp[0].val); }
#line 4965 "y.tab.c" /* yacc.c:1646  */
    break;

  case 64:
#line 772 "gram.y" /* yacc.c:1646  */
    { sysstr.color = (int) (yyvsp[0].val); }
#line 4971 "y.tab.c" /* yacc.c:1646  */
    break;

  case 65:
#line 773 "gram.y" /* yacc.c:1646  */
    { sysstr.rot = (int) (yyvsp[0].val); }
#line 4977 "y.tab.c" /* yacc.c:1646  */
    break;

  case 66:
#line 774 "gram.y" /* yacc.c:1646  */
    { sysstr.font = (int) (yyvsp[0].val); }
#line 4983 "y.tab.c" /* yacc.c:1646  */
    break;

  case 67:
#line 775 "gram.y" /* yacc.c:1646  */
    { sysstr.just = (int) (yyvsp[0].val); }
#line 4989 "y.tab.c" /* yacc.c:1646  */
    break;

  case 68:
#line 776 "gram.y" /* yacc.c:1646  */
    { sysstr.sym = (int) (yyvsp[0].val); }
#line 4995 "y.tab.c" /* yacc.c:1646  */
    break;

  case 69:
#line 777 "gram.y" /* yacc.c:1646  */
    { sysstr.symloc = (int) (yyvsp[0].pset); }
#line 5001 "y.tab.c" /* yacc.c:1646  */
    break;

  case 70:
#line 778 "gram.y" /* yacc.c:1646  */
    { sysstr.symsize = (double) (yyvsp[0].val); }
#line 5007 "y.tab.c" /* yacc.c:1646  */
    break;

  case 71:
#line 779 "gram.y" /* yacc.c:1646  */
    { sysstr.symfill = (int) (yyvsp[0].val); }
#line 5013 "y.tab.c" /* yacc.c:1646  */
    break;

  case 72:
#line 780 "gram.y" /* yacc.c:1646  */
    { sysstr.symcolor = (int) (yyvsp[0].val); }
#line 5019 "y.tab.c" /* yacc.c:1646  */
    break;

  case 73:
#line 781 "gram.y" /* yacc.c:1646  */
    { sysstr.charsize = (double) (yyvsp[0].val); }
#line 5025 "y.tab.c" /* yacc.c:1646  */
    break;

  case 74:
#line 783 "gram.y" /* yacc.c:1646  */
    {
	    strcpy(pstr[curstring].s, (char *) (yyvsp[0].str));
	    pstr[curstring].linew = sysstr.linew;
	    pstr[curstring].color = sysstr.color;
	    pstr[curstring].font = sysstr.font;
	    pstr[curstring].just = sysstr.just;
	    pstr[curstring].sym = sysstr.sym;
	    pstr[curstring].symloc = sysstr.symloc;
	    pstr[curstring].symfill = sysstr.symfill;
	    pstr[curstring].symcolor = sysstr.symcolor;
	    pstr[curstring].symsize = sysstr.symsize;
	    pstr[curstring].loctype = sysstr.loctype;
	    pstr[curstring].rot = sysstr.rot;
	    pstr[curstring].charsize = sysstr.charsize;
	    /*print_plotstr(pstr[curstring]);*/
	}
#line 5046 "y.tab.c" /* yacc.c:1646  */
    break;

  case 75:
#line 802 "gram.y" /* yacc.c:1646  */
    { setistep(); }
#line 5052 "y.tab.c" /* yacc.c:1646  */
    break;

  case 76:
#line 803 "gram.y" /* yacc.c:1646  */
    { setreverse(); }
#line 5058 "y.tab.c" /* yacc.c:1646  */
    break;

  case 77:
#line 804 "gram.y" /* yacc.c:1646  */
    { setrewind(); }
#line 5064 "y.tab.c" /* yacc.c:1646  */
    break;

  case 78:
#line 805 "gram.y" /* yacc.c:1646  */
    { setforward(); }
#line 5070 "y.tab.c" /* yacc.c:1646  */
    break;

  case 79:
#line 806 "gram.y" /* yacc.c:1646  */
    { set_wrap((yyvsp[0].pset) == ON); }
#line 5076 "y.tab.c" /* yacc.c:1646  */
    break;

  case 80:
#line 807 "gram.y" /* yacc.c:1646  */
    { goto_step((int) (yyvsp[0].val) - 1); }
#line 5082 "y.tab.c" /* yacc.c:1646  */
    break;

  case 81:
#line 808 "gram.y" /* yacc.c:1646  */
    { setirun(); }
#line 5088 "y.tab.c" /* yacc.c:1646  */
    break;

  case 82:
#line 809 "gram.y" /* yacc.c:1646  */
    { runsteps((int) (yyvsp[-2].val), (int) (yyvsp[0].val)); }
#line 5094 "y.tab.c" /* yacc.c:1646  */
    break;

  case 83:
#line 810 "gram.y" /* yacc.c:1646  */
    { batchrunsteps((int) (yyvsp[-2].val), (int) (yyvsp[0].val), 1); }
#line 5100 "y.tab.c" /* yacc.c:1646  */
    break;

  case 84:
#line 811 "gram.y" /* yacc.c:1646  */
    { batchrunsteps((int) (yyvsp[-4].val), (int) (yyvsp[-2].val), (int) (yyvsp[0].val)); }
#line 5106 "y.tab.c" /* yacc.c:1646  */
    break;

  case 85:
#line 812 "gram.y" /* yacc.c:1646  */
    { strcpy(batchprefix, (char *) (yyvsp[0].str)); }
#line 5112 "y.tab.c" /* yacc.c:1646  */
    break;

  case 86:
#line 813 "gram.y" /* yacc.c:1646  */
    { setistop(); }
#line 5118 "y.tab.c" /* yacc.c:1646  */
    break;

  case 87:
#line 817 "gram.y" /* yacc.c:1646  */
    { system((yyvsp[0].str)); }
#line 5124 "y.tab.c" /* yacc.c:1646  */
    break;

  case 88:
#line 819 "gram.y" /* yacc.c:1646  */
    {
	    gotbatch = 1;
	    batchfile[0] = 0;
	    strcpy(batchfile, (yyvsp[0].str));
	}
#line 5134 "y.tab.c" /* yacc.c:1646  */
    break;

  case 89:
#line 824 "gram.y" /* yacc.c:1646  */
    {
		if (chdir((char *) (yyvsp[0].str)) < 0) {
			sprintf(buf, "chdir() to %s failed", (char *) (yyvsp[0].str));
			errwin(buf);
		}
	}
#line 5145 "y.tab.c" /* yacc.c:1646  */
    break;

  case 90:
#line 830 "gram.y" /* yacc.c:1646  */
    { setreset_world(); }
#line 5151 "y.tab.c" /* yacc.c:1646  */
    break;

  case 91:
#line 831 "gram.y" /* yacc.c:1646  */
    { setredraw_world(); }
#line 5157 "y.tab.c" /* yacc.c:1646  */
    break;

  case 92:
#line 832 "gram.y" /* yacc.c:1646  */
    { sleep((int) (yyvsp[0].val)); }
#line 5163 "y.tab.c" /* yacc.c:1646  */
    break;

  case 93:
#line 833 "gram.y" /* yacc.c:1646  */
    { exit(0); }
#line 5169 "y.tab.c" /* yacc.c:1646  */
    break;

  case 94:
#line 834 "gram.y" /* yacc.c:1646  */
    { my_blowup((yyvsp[-6].val), (yyvsp[-4].val), (yyvsp[-2].val), (yyvsp[0].val)); }
#line 5175 "y.tab.c" /* yacc.c:1646  */
    break;

  case 95:
#line 835 "gram.y" /* yacc.c:1646  */
    { page(page_per, 4); }
#line 5181 "y.tab.c" /* yacc.c:1646  */
    break;

  case 96:
#line 836 "gram.y" /* yacc.c:1646  */
    { page(page_per, 5); }
#line 5187 "y.tab.c" /* yacc.c:1646  */
    break;

  case 97:
#line 837 "gram.y" /* yacc.c:1646  */
    { page(page_per, 0); }
#line 5193 "y.tab.c" /* yacc.c:1646  */
    break;

  case 98:
#line 838 "gram.y" /* yacc.c:1646  */
    { page(page_per, 1); }
#line 5199 "y.tab.c" /* yacc.c:1646  */
    break;

  case 99:
#line 839 "gram.y" /* yacc.c:1646  */
    { page(page_per, 2); }
#line 5205 "y.tab.c" /* yacc.c:1646  */
    break;

  case 100:
#line 840 "gram.y" /* yacc.c:1646  */
    { page(page_per, 3); }
#line 5211 "y.tab.c" /* yacc.c:1646  */
    break;

  case 101:
#line 841 "gram.y" /* yacc.c:1646  */
    { page_per = (yyvsp[0].val); }
#line 5217 "y.tab.c" /* yacc.c:1646  */
    break;

  case 102:
#line 842 "gram.y" /* yacc.c:1646  */
    { scrollinout_proc((int) (yyvsp[0].val)); }
#line 5223 "y.tab.c" /* yacc.c:1646  */
    break;

  case 103:
#line 843 "gram.y" /* yacc.c:1646  */
    { scrolling_islinked = (yyvsp[0].pset) == ON; }
#line 5229 "y.tab.c" /* yacc.c:1646  */
    break;

  case 104:
#line 845 "gram.y" /* yacc.c:1646  */
    {
	    if ((logfp = fopen((yyvsp[0].str), "w")) != NULL) {
		logfile = 1;
		printf("Opened logfile %s\n", (yyvsp[0].str));
	    } else {
		logfile = 0;
		printf("Failed to open logfile %s\n", (yyvsp[0].str));
	    }
	}
#line 5243 "y.tab.c" /* yacc.c:1646  */
    break;

  case 105:
#line 855 "gram.y" /* yacc.c:1646  */
    {
	    if (logfp != NULL) {
		printf("Closing logfile\n");
		logfile = 0;
		fclose(logfp);
		logfp = NULL;
	    }
	}
#line 5256 "y.tab.c" /* yacc.c:1646  */
    break;

  case 106:
#line 863 "gram.y" /* yacc.c:1646  */
    { }
#line 5262 "y.tab.c" /* yacc.c:1646  */
    break;

  case 107:
#line 864 "gram.y" /* yacc.c:1646  */
    { batchrunstep(0); }
#line 5268 "y.tab.c" /* yacc.c:1646  */
    break;

  case 108:
#line 865 "gram.y" /* yacc.c:1646  */
    {
	    if (inwin) { set_left_footer((yyvsp[0].str)); }
	    else { printf("%s\n", (yyvsp[0].str)); }
	}
#line 5277 "y.tab.c" /* yacc.c:1646  */
    break;

  case 109:
#line 870 "gram.y" /* yacc.c:1646  */
    { set_colormapdata((int) (yyvsp[-6].val), (int) (yyvsp[-4].val), (int) (yyvsp[-2].val), (int) (yyvsp[0].val)); }
#line 5283 "y.tab.c" /* yacc.c:1646  */
    break;

  case 110:
#line 872 "gram.y" /* yacc.c:1646  */
    { set_colormapdata((int) (yyvsp[-6].val), (int) (yyvsp[-4].val), (int) (yyvsp[-2].val), (int) (yyvsp[0].val)); }
#line 5289 "y.tab.c" /* yacc.c:1646  */
    break;

  case 111:
#line 874 "gram.y" /* yacc.c:1646  */
    { set_colormapdata((int) (yyvsp[-6].val), (int) (yyvsp[-4].val), (int) (yyvsp[-2].val), (int) (yyvsp[0].val)); }
#line 5295 "y.tab.c" /* yacc.c:1646  */
    break;

  case 112:
#line 875 "gram.y" /* yacc.c:1646  */
    { setisol = &(g[curg].salip); }
#line 5301 "y.tab.c" /* yacc.c:1646  */
    break;

  case 114:
#line 876 "gram.y" /* yacc.c:1646  */
    { setisol = &(g[curg].velmagip); }
#line 5307 "y.tab.c" /* yacc.c:1646  */
    break;

  case 116:
#line 880 "gram.y" /* yacc.c:1646  */
    { 
		setflow = &g[curg].flowf[(int) (yyvsp[0].val)]; 
	}
#line 5315 "y.tab.c" /* yacc.c:1646  */
    break;

  case 117:
#line 883 "gram.y" /* yacc.c:1646  */
    { }
#line 5321 "y.tab.c" /* yacc.c:1646  */
    break;

  case 118:
#line 884 "gram.y" /* yacc.c:1646  */
    {}
#line 5327 "y.tab.c" /* yacc.c:1646  */
    break;

  case 119:
#line 885 "gram.y" /* yacc.c:1646  */
    { 
		setflow = &g[curg].flowt[(int) (yyvsp[0].val)]; 
                elcirc_flowno = (int) (yyvsp[0].val);
	}
#line 5336 "y.tab.c" /* yacc.c:1646  */
    break;

  case 120:
#line 889 "gram.y" /* yacc.c:1646  */
    {}
#line 5342 "y.tab.c" /* yacc.c:1646  */
    break;

  case 121:
#line 891 "gram.y" /* yacc.c:1646  */
    {
	      int fno = (int) (yyvsp[-3].val);
              readbin_adcirc_elev(fno, (char *) (yyvsp[-1].str), (char *) (yyvsp[0].str));
              set_clock(0, flowt[fno].start, 
				flowt[fno].stop, flowt[fno].step,
                                  flowt[fno].nsteps);
              load_clock(ADCIRC, fno);
          }
#line 5355 "y.tab.c" /* yacc.c:1646  */
    break;

  case 122:
#line 899 "gram.y" /* yacc.c:1646  */
    { 
		setflow = &g[curg].flowt[(int) (yyvsp[0].val)]; 
                elcirc_flowno = (int) (yyvsp[0].val);
		g[curg].curadc3d = (int) (yyvsp[0].val); 
	}
#line 5365 "y.tab.c" /* yacc.c:1646  */
    break;

  case 123:
#line 904 "gram.y" /* yacc.c:1646  */
    {}
#line 5371 "y.tab.c" /* yacc.c:1646  */
    break;

  case 124:
#line 905 "gram.y" /* yacc.c:1646  */
    { curadc3d = (int) (yyvsp[0].val); }
#line 5377 "y.tab.c" /* yacc.c:1646  */
    break;

  case 125:
#line 907 "gram.y" /* yacc.c:1646  */
    {
		ReplaceElcircGrid((int) (yyvsp[-1].val), (char *) (yyvsp[0].str));
	}
#line 5385 "y.tab.c" /* yacc.c:1646  */
    break;

  case 126:
#line 911 "gram.y" /* yacc.c:1646  */
    { 
		ReadElcircRegionFile((int) (yyvsp[-2].val), (char *) (yyvsp[0].str));
	}
#line 5393 "y.tab.c" /* yacc.c:1646  */
    break;

  case 127:
#line 915 "gram.y" /* yacc.c:1646  */
    {
		ReadElcirc((int) (yyvsp[-9].val), (char *) (yyvsp[-8].str), (int) (yyvsp[0].val) - 1, (int) (yyvsp[-6].val), (int) (yyvsp[-4].val), (int) (yyvsp[-2].val), 0, 0, 0);
                elcirc_flowno = (int) (yyvsp[-9].val);
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
#line 5406 "y.tab.c" /* yacc.c:1646  */
    break;

  case 128:
#line 924 "gram.y" /* yacc.c:1646  */
    {
		ReadElcirc((int) (yyvsp[-10].val), (char *) (yyvsp[-9].str), (int) (yyvsp[-1].val) - 1, (int) (yyvsp[-7].val), (int) (yyvsp[-5].val), (int) (yyvsp[-3].val), 0, 0, 1);
                elcirc_flowno = (int) (yyvsp[-10].val);
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
#line 5419 "y.tab.c" /* yacc.c:1646  */
    break;

  case 129:
#line 933 "gram.y" /* yacc.c:1646  */
    {
		ReadElcircDepth((int) (yyvsp[-9].val), (char *) (yyvsp[-8].str), (char *) NULL, (double) (yyvsp[0].val), (int) (yyvsp[-6].val), (int) (yyvsp[-4].val), (int) (yyvsp[-2].val), 0, 0, 0);
                elcirc_flowno = (int) (yyvsp[-9].val);
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
#line 5432 "y.tab.c" /* yacc.c:1646  */
    break;

  case 130:
#line 942 "gram.y" /* yacc.c:1646  */
    {
/* read at a given depth relative to the free surface */
		ReadElcircDepthFromFreeSurface((int) (yyvsp[-9].val), (char *) (yyvsp[-8].str), (char *) NULL, (double) (yyvsp[0].val), (int) (yyvsp[-6].val), (int) (yyvsp[-4].val), (int) (yyvsp[-2].val), 0, 0, 0);
                elcirc_flowno = (int) (yyvsp[-9].val);
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
#line 5446 "y.tab.c" /* yacc.c:1646  */
    break;

  case 131:
#line 952 "gram.y" /* yacc.c:1646  */
    {
		ReadElcircSurf((int) (yyvsp[-7].val), (char *) (yyvsp[-6].str), 0, (int) (yyvsp[-4].val), (int) (yyvsp[-2].val), (int) (yyvsp[0].val), 0, 0.0, 0);
                elcirc_flowno = (int) (yyvsp[-7].val);
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
#line 5459 "y.tab.c" /* yacc.c:1646  */
    break;

  case 132:
#line 961 "gram.y" /* yacc.c:1646  */
    {
		ReadElcircSurf((int) (yyvsp[-8].val), (char *) (yyvsp[-7].str), 0, (int) (yyvsp[-5].val), (int) (yyvsp[-3].val), (int) (yyvsp[-1].val), 0, 0.0, 1);
                elcirc_flowno = (int) (yyvsp[-8].val);
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
#line 5472 "y.tab.c" /* yacc.c:1646  */
    break;

  case 133:
#line 970 "gram.y" /* yacc.c:1646  */
    {
		ReadElcircSurf((int) (yyvsp[-7].val), (char *) (yyvsp[-6].str), 1, (int) (yyvsp[-4].val), (int) (yyvsp[-2].val), (int) (yyvsp[0].val), 0, 0.0, 0);
                elcirc_flowno = (int) (yyvsp[-7].val);
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
#line 5485 "y.tab.c" /* yacc.c:1646  */
    break;

  case 134:
#line 979 "gram.y" /* yacc.c:1646  */
    {
		ReadElcircSurf((int) (yyvsp[-8].val), (char *) (yyvsp[-7].str), 1, (int) (yyvsp[-5].val), (int) (yyvsp[-3].val), (int) (yyvsp[-1].val), 0, 0.0, 1);
                elcirc_flowno = (int) (yyvsp[-8].val);
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
#line 5498 "y.tab.c" /* yacc.c:1646  */
    break;

  case 135:
#line 987 "gram.y" /* yacc.c:1646  */
    { 
			curtrans = (int) (yyvsp[0].val); 
			settrans = &(trans[curtrans]); 
	}
#line 5507 "y.tab.c" /* yacc.c:1646  */
    break;

  case 136:
#line 991 "gram.y" /* yacc.c:1646  */
    { curtrans = (int) (yyvsp[0].val); }
#line 5513 "y.tab.c" /* yacc.c:1646  */
    break;

  case 137:
#line 992 "gram.y" /* yacc.c:1646  */
    { setisol = &(g[curg].trans[curtrans].ip); }
#line 5519 "y.tab.c" /* yacc.c:1646  */
    break;

  case 139:
#line 993 "gram.y" /* yacc.c:1646  */
    {strcpy(settrans->uvname, (char *) (yyvsp[0].str));}
#line 5525 "y.tab.c" /* yacc.c:1646  */
    break;

  case 140:
#line 994 "gram.y" /* yacc.c:1646  */
    {strcpy(settrans->vvname, (char *) (yyvsp[0].str));}
#line 5531 "y.tab.c" /* yacc.c:1646  */
    break;

  case 141:
#line 995 "gram.y" /* yacc.c:1646  */
    { settrans->flowno = (int) (yyvsp[0].val);}
#line 5537 "y.tab.c" /* yacc.c:1646  */
    break;

  case 142:
#line 996 "gram.y" /* yacc.c:1646  */
    {strcpy(settrans->salname, (char *) (yyvsp[0].str));}
#line 5543 "y.tab.c" /* yacc.c:1646  */
    break;

  case 143:
#line 997 "gram.y" /* yacc.c:1646  */
    {strcpy(settrans->elevname, (char *) (yyvsp[0].str));}
#line 5549 "y.tab.c" /* yacc.c:1646  */
    break;

  case 144:
#line 998 "gram.y" /* yacc.c:1646  */
    { settrans->gno = (int) (yyvsp[0].val);}
#line 5555 "y.tab.c" /* yacc.c:1646  */
    break;

  case 145:
#line 999 "gram.y" /* yacc.c:1646  */
    { settrans->transgno = (int) (yyvsp[0].val);}
#line 5561 "y.tab.c" /* yacc.c:1646  */
    break;

  case 146:
#line 1000 "gram.y" /* yacc.c:1646  */
    { settrans->display = (int) (yyvsp[0].pset);}
#line 5567 "y.tab.c" /* yacc.c:1646  */
    break;

  case 147:
#line 1001 "gram.y" /* yacc.c:1646  */
    { g[curg].trans[curtrans].display = (int) (yyvsp[0].pset);}
#line 5573 "y.tab.c" /* yacc.c:1646  */
    break;

  case 148:
#line 1002 "gram.y" /* yacc.c:1646  */
    { g[curg].trans[curtrans].display_mag = (int) (yyvsp[0].pset);}
#line 5579 "y.tab.c" /* yacc.c:1646  */
    break;

  case 149:
#line 1003 "gram.y" /* yacc.c:1646  */
    { }
#line 5585 "y.tab.c" /* yacc.c:1646  */
    break;

  case 150:
#line 1005 "gram.y" /* yacc.c:1646  */
    {
		settrans->start = (int) (yyvsp[-4].val);
		settrans->stop = (int) (yyvsp[-2].val);
		settrans->skip = (int) (yyvsp[0].val);
	}
#line 5595 "y.tab.c" /* yacc.c:1646  */
    break;

  case 151:
#line 1010 "gram.y" /* yacc.c:1646  */
    { settrans->npts = (int) (yyvsp[0].val); }
#line 5601 "y.tab.c" /* yacc.c:1646  */
    break;

  case 152:
#line 1011 "gram.y" /* yacc.c:1646  */
    { settrans->transtype = (int) (yyvsp[0].val); }
#line 5607 "y.tab.c" /* yacc.c:1646  */
    break;

  case 153:
#line 1012 "gram.y" /* yacc.c:1646  */
    { AddTransNXY(settrans, (int) (yyvsp[-4].val), (double) (yyvsp[-2].val), (double) (yyvsp[0].val)); }
#line 5613 "y.tab.c" /* yacc.c:1646  */
    break;

  case 154:
#line 1013 "gram.y" /* yacc.c:1646  */
    { AddTransNode(settrans, (int) (yyvsp[-2].val), (int) (yyvsp[0].val)); }
#line 5619 "y.tab.c" /* yacc.c:1646  */
    break;

  case 155:
#line 1015 "gram.y" /* yacc.c:1646  */
    {
		settrans->transtype = 0;
		strcpy(settrans->transname, (char *) (yyvsp[0].str));
	}
#line 5628 "y.tab.c" /* yacc.c:1646  */
    break;

  case 156:
#line 1020 "gram.y" /* yacc.c:1646  */
    {
		settrans->transtype = 1;
		settrans->npts = (int) (yyvsp[-8].val);
		settrans->x1 = (double) (yyvsp[-6].val);
		settrans->y1 = (double) (yyvsp[-4].val);
		settrans->x2 = (double) (yyvsp[-2].val);
		settrans->y2 = (double) (yyvsp[0].val);
	}
#line 5641 "y.tab.c" /* yacc.c:1646  */
    break;

  case 157:
#line 1029 "gram.y" /* yacc.c:1646  */
    { 
		settrans->type = VECTOR;
		ReadNewTrans(settrans, 0);
                trans[curtrans].active = ON;
                set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
                load_clock(TRANSECT, curtrans);
	}
#line 5653 "y.tab.c" /* yacc.c:1646  */
    break;

  case 158:
#line 1037 "gram.y" /* yacc.c:1646  */
    { 
		settrans->type = VECTOR;
		ReadNewTrans(settrans, 1);
                trans[curtrans].active = ON;
                set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
                load_clock(TRANSECT, curtrans);
	}
#line 5665 "y.tab.c" /* yacc.c:1646  */
    break;

  case 159:
#line 1045 "gram.y" /* yacc.c:1646  */
    { 
		settrans->type = VECTOR;
		ReadNewTrans(settrans, 0);
                trans[curtrans].active = ON;
                set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
                load_clock(TRANSECT, curtrans);
	}
#line 5677 "y.tab.c" /* yacc.c:1646  */
    break;

  case 160:
#line 1053 "gram.y" /* yacc.c:1646  */
    { 
		settrans->type = VECTOR;
		ReadNewTrans(settrans, 1);
                trans[curtrans].active = ON;
                set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
                load_clock(TRANSECT, curtrans);
	}
#line 5689 "y.tab.c" /* yacc.c:1646  */
    break;

  case 161:
#line 1061 "gram.y" /* yacc.c:1646  */
    { 
		settrans->type = SCALAR;
		ReadNewTrans(settrans, 0);
                trans[curtrans].active = ON;
                set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
                load_clock(TRANSECT, curtrans);
	}
#line 5701 "y.tab.c" /* yacc.c:1646  */
    break;

  case 162:
#line 1069 "gram.y" /* yacc.c:1646  */
    { 
		settrans->type = SCALAR;
		ReadNewTrans(settrans, 1);
                trans[curtrans].active = ON;
                set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
                load_clock(TRANSECT, curtrans);
	}
#line 5713 "y.tab.c" /* yacc.c:1646  */
    break;

  case 163:
#line 1077 "gram.y" /* yacc.c:1646  */
    {
		elcircmarker = (int) (yyvsp[0].val);
		setadc3d = &adc3d[(int) (yyvsp[0].val)];
		setflow3d = &g[curg].flow3d[(int) (yyvsp[0].val)];
	}
#line 5723 "y.tab.c" /* yacc.c:1646  */
    break;

  case 164:
#line 1083 "gram.y" /* yacc.c:1646  */
    {
		strcpy(setadc3d->datafile, (char *) (yyvsp[0].str));
	}
#line 5731 "y.tab.c" /* yacc.c:1646  */
    break;

  case 165:
#line 1087 "gram.y" /* yacc.c:1646  */
    {
		strcpy(setadc3d->datafile, (char *) (yyvsp[0].str));
	}
#line 5739 "y.tab.c" /* yacc.c:1646  */
    break;

  case 166:
#line 1091 "gram.y" /* yacc.c:1646  */
    {
		strcpy(setadc3d->datafile, (char *) (yyvsp[0].str));
	}
#line 5747 "y.tab.c" /* yacc.c:1646  */
    break;

  case 167:
#line 1095 "gram.y" /* yacc.c:1646  */
    {
		strcpy(setadc3d->datafile, (char *) (yyvsp[0].str));
	}
#line 5755 "y.tab.c" /* yacc.c:1646  */
    break;

  case 168:
#line 1099 "gram.y" /* yacc.c:1646  */
    {
		strcpy(setadc3d->datafile, (char *) (yyvsp[0].str));
	}
#line 5763 "y.tab.c" /* yacc.c:1646  */
    break;

  case 169:
#line 1103 "gram.y" /* yacc.c:1646  */
    {
		strcpy(setadc3d->elevfile, (char *) (yyvsp[0].str));
	}
#line 5771 "y.tab.c" /* yacc.c:1646  */
    break;

  case 170:
#line 1107 "gram.y" /* yacc.c:1646  */
    {
		setadc3d->loctype = NODE;
		setadc3d->loctype = 1;
		setadc3d->node = (int) (yyvsp[-6].val);
		ReadNodeDataNew(elcircmarker, (int) (yyvsp[-4].val), (int) (yyvsp[-2].val), 0);
	}
#line 5782 "y.tab.c" /* yacc.c:1646  */
    break;

  case 171:
#line 1114 "gram.y" /* yacc.c:1646  */
    {
		setadc3d->loctype = NODE;
		setadc3d->loctype = 1;
		setadc3d->node = (int) (yyvsp[-7].val);
		ReadNodeDataNew(elcircmarker, (int) (yyvsp[-5].val), (int) (yyvsp[-3].val), 1);
	}
#line 5793 "y.tab.c" /* yacc.c:1646  */
    break;

  case 172:
#line 1121 "gram.y" /* yacc.c:1646  */
    {
		setadc3d->loctype = XY;
		setadc3d->loctype = 0;
		setadc3d->x = (double) (yyvsp[-8].val);
		setadc3d->y = (double) (yyvsp[-6].val);
		ReadXYDataNew(elcircmarker, (int) (yyvsp[-4].val), (int) (yyvsp[-2].val), 0);
	}
#line 5805 "y.tab.c" /* yacc.c:1646  */
    break;

  case 173:
#line 1129 "gram.y" /* yacc.c:1646  */
    {
		setadc3d->loctype = XY;
		setadc3d->loctype = 0;
		setadc3d->x = (double) (yyvsp[-9].val);
		setadc3d->y = (double) (yyvsp[-7].val);
		ReadXYDataNew(elcircmarker, (int) (yyvsp[-5].val), (int) (yyvsp[-3].val), 1);
	}
#line 5817 "y.tab.c" /* yacc.c:1646  */
    break;

  case 174:
#line 1136 "gram.y" /* yacc.c:1646  */
    { setisol = &(g[curg].salip); }
#line 5823 "y.tab.c" /* yacc.c:1646  */
    break;

  case 176:
#line 1137 "gram.y" /* yacc.c:1646  */
    { setisol = &(g[curg].salip); }
#line 5829 "y.tab.c" /* yacc.c:1646  */
    break;

  case 178:
#line 1138 "gram.y" /* yacc.c:1646  */
    { setisol = &(g[curg].velmagip); }
#line 5835 "y.tab.c" /* yacc.c:1646  */
    break;

  case 180:
#line 1139 "gram.y" /* yacc.c:1646  */
    { setprops = &(setflow3d->p); }
#line 5841 "y.tab.c" /* yacc.c:1646  */
    break;

  case 182:
#line 1141 "gram.y" /* yacc.c:1646  */
    {
		setflow3d->precx = (int) (yyvsp[-2].val);
		setflow3d->precy = (int) (yyvsp[0].val);
	}
#line 5850 "y.tab.c" /* yacc.c:1646  */
    break;

  case 183:
#line 1145 "gram.y" /* yacc.c:1646  */
    { setflow3d->attach = (int) (yyvsp[0].val); }
#line 5856 "y.tab.c" /* yacc.c:1646  */
    break;

  case 184:
#line 1146 "gram.y" /* yacc.c:1646  */
    { setflow3d->loctype = (int) (yyvsp[0].pset); }
#line 5862 "y.tab.c" /* yacc.c:1646  */
    break;

  case 185:
#line 1147 "gram.y" /* yacc.c:1646  */
    { setflow3d->display_marker = (int) (yyvsp[0].pset); }
#line 5868 "y.tab.c" /* yacc.c:1646  */
    break;

  case 186:
#line 1148 "gram.y" /* yacc.c:1646  */
    { setflow3d->display = (int) (yyvsp[0].pset); }
#line 5874 "y.tab.c" /* yacc.c:1646  */
    break;

  case 187:
#line 1149 "gram.y" /* yacc.c:1646  */
    { setflow3d->p.color = (yyvsp[0].val); }
#line 5880 "y.tab.c" /* yacc.c:1646  */
    break;

  case 188:
#line 1150 "gram.y" /* yacc.c:1646  */
    { setflow3d->p.linew = (yyvsp[0].val); }
#line 5886 "y.tab.c" /* yacc.c:1646  */
    break;

  case 189:
#line 1151 "gram.y" /* yacc.c:1646  */
    { setflow3d->p.fillcol = (yyvsp[0].val); }
#line 5892 "y.tab.c" /* yacc.c:1646  */
    break;

  case 190:
#line 1153 "gram.y" /* yacc.c:1646  */
    { 
		setflow3d->wx1 = (double) (yyvsp[-6].val); 
		setflow3d->wy1 = (double) (yyvsp[-4].val); 
		setflow3d->wx2 = (double) (yyvsp[-2].val); 
		setflow3d->wy2 = (double) (yyvsp[0].val); 
	}
#line 5903 "y.tab.c" /* yacc.c:1646  */
    break;

  case 191:
#line 1160 "gram.y" /* yacc.c:1646  */
    { 
		setflow3d->vx = (double) (yyvsp[-2].val); 
		setflow3d->vy = (double) (yyvsp[0].val); 
	}
#line 5912 "y.tab.c" /* yacc.c:1646  */
    break;

  case 192:
#line 1165 "gram.y" /* yacc.c:1646  */
    { 
		setflow3d->locx = (double) (yyvsp[-2].val); 
		setflow3d->locy = (double) (yyvsp[0].val); 
	}
#line 5921 "y.tab.c" /* yacc.c:1646  */
    break;

  case 193:
#line 1170 "gram.y" /* yacc.c:1646  */
    { 
		setflow3d->x = (double) (yyvsp[-2].val); 
		setflow3d->y = (double) (yyvsp[0].val); 
	}
#line 5930 "y.tab.c" /* yacc.c:1646  */
    break;

  case 194:
#line 1177 "gram.y" /* yacc.c:1646  */
    { setflow->display = (yyvsp[0].pset);  }
#line 5936 "y.tab.c" /* yacc.c:1646  */
    break;

  case 195:
#line 1178 "gram.y" /* yacc.c:1646  */
    { setflow->display_elev = (yyvsp[0].pset); }
#line 5942 "y.tab.c" /* yacc.c:1646  */
    break;

  case 196:
#line 1179 "gram.y" /* yacc.c:1646  */
    { setflow->display_elevdepth = (yyvsp[0].pset); }
#line 5948 "y.tab.c" /* yacc.c:1646  */
    break;

  case 197:
#line 1180 "gram.y" /* yacc.c:1646  */
    { setflow->display_maxelevval = (yyvsp[0].pset); }
#line 5954 "y.tab.c" /* yacc.c:1646  */
    break;

  case 198:
#line 1181 "gram.y" /* yacc.c:1646  */
    { setflow->display_maxelev = (yyvsp[0].pset); }
#line 5960 "y.tab.c" /* yacc.c:1646  */
    break;

  case 199:
#line 1182 "gram.y" /* yacc.c:1646  */
    { setflow->display_amp = (yyvsp[0].pset); }
#line 5966 "y.tab.c" /* yacc.c:1646  */
    break;

  case 200:
#line 1183 "gram.y" /* yacc.c:1646  */
    { setflow->display_phase = (yyvsp[0].pset); }
#line 5972 "y.tab.c" /* yacc.c:1646  */
    break;

  case 201:
#line 1184 "gram.y" /* yacc.c:1646  */
    { setflow->display_elevmarkers = (yyvsp[0].pset); }
#line 5978 "y.tab.c" /* yacc.c:1646  */
    break;

  case 202:
#line 1185 "gram.y" /* yacc.c:1646  */
    { setflow->display_mag = (yyvsp[0].pset); }
#line 5984 "y.tab.c" /* yacc.c:1646  */
    break;

  case 203:
#line 1186 "gram.y" /* yacc.c:1646  */
    { setflow->display_wind = (yyvsp[0].pset); }
#line 5990 "y.tab.c" /* yacc.c:1646  */
    break;

  case 204:
#line 1187 "gram.y" /* yacc.c:1646  */
    { setflow->display_inun = (yyvsp[-1].pset); }
#line 5996 "y.tab.c" /* yacc.c:1646  */
    break;

  case 205:
#line 1188 "gram.y" /* yacc.c:1646  */
    { setflow->p.color = (yyvsp[0].val); }
#line 6002 "y.tab.c" /* yacc.c:1646  */
    break;

  case 206:
#line 1189 "gram.y" /* yacc.c:1646  */
    { setisol = &(setflow->elevip); }
#line 6008 "y.tab.c" /* yacc.c:1646  */
    break;

  case 208:
#line 1190 "gram.y" /* yacc.c:1646  */
    { setisol = &(setflow->maxelevip); }
#line 6014 "y.tab.c" /* yacc.c:1646  */
    break;

  case 210:
#line 1191 "gram.y" /* yacc.c:1646  */
    { setisol = &(setflow->ampip); }
#line 6020 "y.tab.c" /* yacc.c:1646  */
    break;

  case 212:
#line 1192 "gram.y" /* yacc.c:1646  */
    { setisol = &(setflow->phaseip); }
#line 6026 "y.tab.c" /* yacc.c:1646  */
    break;

  case 214:
#line 1193 "gram.y" /* yacc.c:1646  */
    { setisol = &(setflow->magip); }
#line 6032 "y.tab.c" /* yacc.c:1646  */
    break;

  case 216:
#line 1194 "gram.y" /* yacc.c:1646  */
    { setflow->flowfreq = (yyvsp[0].val); }
#line 6038 "y.tab.c" /* yacc.c:1646  */
    break;

  case 217:
#line 1195 "gram.y" /* yacc.c:1646  */
    { setflow->freq = (yyvsp[0].val); }
#line 6044 "y.tab.c" /* yacc.c:1646  */
    break;

  case 218:
#line 1196 "gram.y" /* yacc.c:1646  */
    { setelevmarker = &(setflow->em[(int) (yyvsp[0].val)]); }
#line 6050 "y.tab.c" /* yacc.c:1646  */
    break;

  case 219:
#line 1197 "gram.y" /* yacc.c:1646  */
    { setflow->sample = (int) (yyvsp[0].pset); }
#line 6056 "y.tab.c" /* yacc.c:1646  */
    break;

  case 220:
#line 1198 "gram.y" /* yacc.c:1646  */
    { ReadSampleFlow(setflow, (char *) (yyvsp[0].str)); }
#line 6062 "y.tab.c" /* yacc.c:1646  */
    break;

  case 221:
#line 1199 "gram.y" /* yacc.c:1646  */
    { SetMinSampleFlow(elcirc_flowno, (double) (yyvsp[0].val)); }
#line 6068 "y.tab.c" /* yacc.c:1646  */
    break;

  case 222:
#line 1200 "gram.y" /* yacc.c:1646  */
    { ReadSampleFlowXY(setflow, (char *) (yyvsp[0].str)); }
#line 6074 "y.tab.c" /* yacc.c:1646  */
    break;

  case 223:
#line 1201 "gram.y" /* yacc.c:1646  */
    { setflow->samptype = XY; }
#line 6080 "y.tab.c" /* yacc.c:1646  */
    break;

  case 224:
#line 1202 "gram.y" /* yacc.c:1646  */
    { setflow->samptype = NODE; }
#line 6086 "y.tab.c" /* yacc.c:1646  */
    break;

  case 225:
#line 1203 "gram.y" /* yacc.c:1646  */
    { AddSampleFlowNode(setflow, (int) (yyvsp[0].val) - 1); }
#line 6092 "y.tab.c" /* yacc.c:1646  */
    break;

  case 226:
#line 1204 "gram.y" /* yacc.c:1646  */
    { AddSampleFlowElem(setflow, (int) (yyvsp[0].val) - 1); }
#line 6098 "y.tab.c" /* yacc.c:1646  */
    break;

  case 227:
#line 1205 "gram.y" /* yacc.c:1646  */
    { AddSampleFlowXY(setflow, (double) (yyvsp[-2].val), (double) (yyvsp[0].val)); }
#line 6104 "y.tab.c" /* yacc.c:1646  */
    break;

  case 228:
#line 1206 "gram.y" /* yacc.c:1646  */
    { DeleteSampleFlowNode(setflow, (int) (yyvsp[0].val) - 1); }
#line 6110 "y.tab.c" /* yacc.c:1646  */
    break;

  case 229:
#line 1207 "gram.y" /* yacc.c:1646  */
    { DeleteSampleFlowElem(setflow, (int) (yyvsp[0].val) - 1); }
#line 6116 "y.tab.c" /* yacc.c:1646  */
    break;

  case 230:
#line 1208 "gram.y" /* yacc.c:1646  */
    { DeleteSampleFlowXY(setflow, (double) (yyvsp[-2].val), (double) (yyvsp[0].val)); }
#line 6122 "y.tab.c" /* yacc.c:1646  */
    break;

  case 231:
#line 1212 "gram.y" /* yacc.c:1646  */
    { setelevmarker = &(setflow->em[(int) (yyvsp[0].val)]); }
#line 6128 "y.tab.c" /* yacc.c:1646  */
    break;

  case 232:
#line 1213 "gram.y" /* yacc.c:1646  */
    { setelevmarker->active = (int) (yyvsp[0].pset); }
#line 6134 "y.tab.c" /* yacc.c:1646  */
    break;

  case 233:
#line 1214 "gram.y" /* yacc.c:1646  */
    { setelevmarker->type = (int) (yyvsp[0].pset); }
#line 6140 "y.tab.c" /* yacc.c:1646  */
    break;

  case 234:
#line 1215 "gram.y" /* yacc.c:1646  */
    { setelevmarker->display = (int) (yyvsp[0].pset); }
#line 6146 "y.tab.c" /* yacc.c:1646  */
    break;

  case 235:
#line 1216 "gram.y" /* yacc.c:1646  */
    { setprops = &(setelevmarker->p); }
#line 6152 "y.tab.c" /* yacc.c:1646  */
    break;

  case 237:
#line 1217 "gram.y" /* yacc.c:1646  */
    { setelevmarker->loctype = (int) (yyvsp[0].pset); }
#line 6158 "y.tab.c" /* yacc.c:1646  */
    break;

  case 238:
#line 1218 "gram.y" /* yacc.c:1646  */
    { setelevmarker->node = (int) (yyvsp[0].val); }
#line 6164 "y.tab.c" /* yacc.c:1646  */
    break;

  case 239:
#line 1220 "gram.y" /* yacc.c:1646  */
    { 
		setelevmarker->locx = (double) (yyvsp[-2].val); 
		setelevmarker->locy = (double) (yyvsp[0].val); 
	}
#line 6173 "y.tab.c" /* yacc.c:1646  */
    break;

  case 240:
#line 1225 "gram.y" /* yacc.c:1646  */
    { 
		setelevmarker->emin = (double) (yyvsp[-2].val); 
		setelevmarker->emax = (double) (yyvsp[0].val); 
	}
#line 6182 "y.tab.c" /* yacc.c:1646  */
    break;

  case 241:
#line 1236 "gram.y" /* yacc.c:1646  */
    { setgrid = &g[curg].grid[(int) (yyvsp[0].val)]; }
#line 6188 "y.tab.c" /* yacc.c:1646  */
    break;

  case 242:
#line 1237 "gram.y" /* yacc.c:1646  */
    { setgrid = &g[curg].grid[(int) (yyvsp[0].val)]; }
#line 6194 "y.tab.c" /* yacc.c:1646  */
    break;

  case 243:
#line 1238 "gram.y" /* yacc.c:1646  */
    { if (checkptr(setgrid, f_string)) setgrid->display = (yyvsp[0].pset); }
#line 6200 "y.tab.c" /* yacc.c:1646  */
    break;

  case 244:
#line 1239 "gram.y" /* yacc.c:1646  */
    { if (checkptr(setgrid, f_string)) setgrid->display_bath = (yyvsp[0].pset); }
#line 6206 "y.tab.c" /* yacc.c:1646  */
    break;

  case 245:
#line 1240 "gram.y" /* yacc.c:1646  */
    { if (checkptr(setgrid, f_string)) setgrid->display_courant = (yyvsp[-2].pset); }
#line 6212 "y.tab.c" /* yacc.c:1646  */
    break;

  case 246:
#line 1241 "gram.y" /* yacc.c:1646  */
    { if (checkptr(setgrid, f_string)) setgrid->display_courantn = (yyvsp[-2].pset); }
#line 6218 "y.tab.c" /* yacc.c:1646  */
    break;

  case 247:
#line 1242 "gram.y" /* yacc.c:1646  */
    { if (checkptr(setgrid, f_string)) setgrid->display_boundary = (yyvsp[0].pset); }
#line 6224 "y.tab.c" /* yacc.c:1646  */
    break;

  case 248:
#line 1243 "gram.y" /* yacc.c:1646  */
    { if (checkptr(setgrid, f_string)) setprops = &(setgrid->p); }
#line 6230 "y.tab.c" /* yacc.c:1646  */
    break;

  case 250:
#line 1244 "gram.y" /* yacc.c:1646  */
    { if (checkptr(setgrid, f_string)) setprops = &(setgrid->bp); }
#line 6236 "y.tab.c" /* yacc.c:1646  */
    break;

  case 252:
#line 1245 "gram.y" /* yacc.c:1646  */
    { if (checkptr(setgrid, f_string)) setgrid->display_nodes = (yyvsp[0].pset); }
#line 6242 "y.tab.c" /* yacc.c:1646  */
    break;

  case 253:
#line 1246 "gram.y" /* yacc.c:1646  */
    { if (checkptr(setgrid, f_string)) setgrid->display_elements = (yyvsp[0].pset); }
#line 6248 "y.tab.c" /* yacc.c:1646  */
    break;

  case 254:
#line 1247 "gram.y" /* yacc.c:1646  */
    { if (checkptr(setgrid, f_string)) setgrid->display_depths = (yyvsp[0].pset); }
#line 6254 "y.tab.c" /* yacc.c:1646  */
    break;

  case 255:
#line 1248 "gram.y" /* yacc.c:1646  */
    { if (checkptr(setgrid, f_string)) setgrid->display_gridf = (yyvsp[0].pset); }
#line 6260 "y.tab.c" /* yacc.c:1646  */
    break;

  case 256:
#line 1249 "gram.y" /* yacc.c:1646  */
    { if (checkptr(setgrid, f_string)) setisol = &(setgrid->ip); }
#line 6266 "y.tab.c" /* yacc.c:1646  */
    break;

  case 258:
#line 1251 "gram.y" /* yacc.c:1646  */
    { 
		autoscale_grid((int) (yyvsp[-1].pset), g[(int) (yyvsp[-1].pset)].curgrid); 
		set_defaults((int) (yyvsp[-1].pset));
	}
#line 6275 "y.tab.c" /* yacc.c:1646  */
    break;

  case 259:
#line 1256 "gram.y" /* yacc.c:1646  */
    { 
		autoscale_grid(curg, g[curg].curgrid); 
		set_defaults(curg);
	}
#line 6284 "y.tab.c" /* yacc.c:1646  */
    break;

  case 260:
#line 1261 "gram.y" /* yacc.c:1646  */
    {
		extern int readgridfile;
		readgrid((int) (yyvsp[-1].val), (char *) (yyvsp[0].str));
		readgridfile = 1;
	}
#line 6294 "y.tab.c" /* yacc.c:1646  */
    break;

  case 261:
#line 1269 "gram.y" /* yacc.c:1646  */
    { setdrogs = &g[curg].drogues[(int) (yyvsp[0].val)]; }
#line 6300 "y.tab.c" /* yacc.c:1646  */
    break;

  case 262:
#line 1270 "gram.y" /* yacc.c:1646  */
    {  setdrogs->display = (int) (yyvsp[0].pset); }
#line 6306 "y.tab.c" /* yacc.c:1646  */
    break;

  case 263:
#line 1271 "gram.y" /* yacc.c:1646  */
    { setprops = &(setdrogs->p); }
#line 6312 "y.tab.c" /* yacc.c:1646  */
    break;

  case 265:
#line 1273 "gram.y" /* yacc.c:1646  */
    {
    	if (!readdrogues(curdrog, (char *) (yyvsp[0].str), -1, 0, 0)) {
        	fprintf(stderr, "Error reading file %s", (char *) (yyvsp[0].str));
    	} else {
        	set_clock(0, drogues[curdrog].start, drogues[curdrog].stop,
                  drogues[curdrog].step,
                  drogues[curdrog].nsteps);
        	load_clock(DROGUES, curdrog);
    	}
	}
#line 6327 "y.tab.c" /* yacc.c:1646  */
    break;

  case 266:
#line 1284 "gram.y" /* yacc.c:1646  */
    {
		readdrogues(0, (char *) (yyvsp[-6].str), (int) (yyvsp[-4].val), (int) (yyvsp[-2].val), (int) (yyvsp[0].val));
        	set_clock(0, drogues[curdrog].start, drogues[curdrog].stop,
                  drogues[curdrog].step,
                  drogues[curdrog].nsteps);
        	load_clock(DROGUES, curdrog);
	}
#line 6339 "y.tab.c" /* yacc.c:1646  */
    break;

  case 267:
#line 1294 "gram.y" /* yacc.c:1646  */
    { setisol->lactive = (yyvsp[0].pset); }
#line 6345 "y.tab.c" /* yacc.c:1646  */
    break;

  case 268:
#line 1295 "gram.y" /* yacc.c:1646  */
    { setisol->layout = (yyvsp[0].pset); }
#line 6351 "y.tab.c" /* yacc.c:1646  */
    break;

  case 269:
#line 1296 "gram.y" /* yacc.c:1646  */
    { setisol->llabels = (yyvsp[-1].pset); }
#line 6357 "y.tab.c" /* yacc.c:1646  */
    break;

  case 270:
#line 1297 "gram.y" /* yacc.c:1646  */
    { setisol->frame = (int) (yyvsp[0].pset); }
#line 6363 "y.tab.c" /* yacc.c:1646  */
    break;

  case 271:
#line 1298 "gram.y" /* yacc.c:1646  */
    { setisol->framecol = (int) (yyvsp[0].val); }
#line 6369 "y.tab.c" /* yacc.c:1646  */
    break;

  case 272:
#line 1299 "gram.y" /* yacc.c:1646  */
    { setisol->nisol = (int) (yyvsp[0].val); }
#line 6375 "y.tab.c" /* yacc.c:1646  */
    break;

  case 273:
#line 1300 "gram.y" /* yacc.c:1646  */
    { setisol->type = (int) (yyvsp[0].val); }
#line 6381 "y.tab.c" /* yacc.c:1646  */
    break;

  case 274:
#line 1301 "gram.y" /* yacc.c:1646  */
    { setisol->isoltype = (int) (yyvsp[0].val); }
#line 6387 "y.tab.c" /* yacc.c:1646  */
    break;

  case 275:
#line 1302 "gram.y" /* yacc.c:1646  */
    { setisol->visflag = (int) (yyvsp[0].val); }
#line 6393 "y.tab.c" /* yacc.c:1646  */
    break;

  case 276:
#line 1303 "gram.y" /* yacc.c:1646  */
    { 
		int i;
		setisol->writeflag = 1; 
		for (i=0;i<MAXISOLINES;i++) {
    		    setisol->writelevel[i] = 0;
		}
		setisol->wname[0] = 0;
	}
#line 6406 "y.tab.c" /* yacc.c:1646  */
    break;

  case 277:
#line 1311 "gram.y" /* yacc.c:1646  */
    { 
		int i;
		setisol->writeflag = 0; 
		for (i=0;i<MAXISOLINES;i++) {
    		    setisol->writelevel[i] = 0;
		}
		setisol->wname[0] = 0;
	}
#line 6419 "y.tab.c" /* yacc.c:1646  */
    break;

  case 278:
#line 1319 "gram.y" /* yacc.c:1646  */
    { 
		int i;
		setisol->writeflag = 1; 
		for (i=0;i<MAXISOLINES;i++) {
    		    setisol->writelevel[i] = 0;
		}
		strncpy(setisol->wname, (char *) (yyvsp[0].str), 1023);
	}
#line 6432 "y.tab.c" /* yacc.c:1646  */
    break;

  case 279:
#line 1327 "gram.y" /* yacc.c:1646  */
    { 
		int i;
		setisol->writeflag = 1; 
		strncpy(setisol->wname, (char *) (yyvsp[0].str), 1023);
		for (i=0;i<MAXISOLINES;i++) {
    		    setisol->writelevel[i] = 0;
		}
    		setisol->writelevel[(int) (yyvsp[-2].val)] = 1;
	}
#line 6446 "y.tab.c" /* yacc.c:1646  */
    break;

  case 280:
#line 1336 "gram.y" /* yacc.c:1646  */
    { 
    		setisol->writelevel[(int) (yyvsp[0].val)] = 1;
	}
#line 6454 "y.tab.c" /* yacc.c:1646  */
    break;

  case 281:
#line 1339 "gram.y" /* yacc.c:1646  */
    {
            setisol->cis[0] = (double) (yyvsp[-2].val);
            setisol->cint = (double) (yyvsp[0].val);
	    if (setisol->isoltype == 0) {
		int i;
		for (i=1;i< 16;i++) {
		    setisol->cis[i] = setisol->cis[0] + i * setisol->cint;
		}
	    }
        }
#line 6469 "y.tab.c" /* yacc.c:1646  */
    break;

  case 282:
#line 1349 "gram.y" /* yacc.c:1646  */
    { setisol->cmin = (yyvsp[-2].val); setisol->cmax = (yyvsp[0].val); }
#line 6475 "y.tab.c" /* yacc.c:1646  */
    break;

  case 283:
#line 1350 "gram.y" /* yacc.c:1646  */
    { setisol->cis[(int) (yyvsp[-2].val)] = (yyvsp[0].val); }
#line 6481 "y.tab.c" /* yacc.c:1646  */
    break;

  case 284:
#line 1351 "gram.y" /* yacc.c:1646  */
    { setisol->loctype = (yyvsp[0].pset); }
#line 6487 "y.tab.c" /* yacc.c:1646  */
    break;

  case 285:
#line 1352 "gram.y" /* yacc.c:1646  */
    {
            setisol->x = (double) (yyvsp[-2].val);
            setisol->y = (double) (yyvsp[0].val);
        }
#line 6496 "y.tab.c" /* yacc.c:1646  */
    break;

  case 286:
#line 1356 "gram.y" /* yacc.c:1646  */
    { setisol->color[(int) (yyvsp[-2].val)] = (int) (yyvsp[0].val); }
#line 6502 "y.tab.c" /* yacc.c:1646  */
    break;

  case 287:
#line 1357 "gram.y" /* yacc.c:1646  */
    { setisol->linew[(int) (yyvsp[-2].val)] = (int) (yyvsp[0].val); }
#line 6508 "y.tab.c" /* yacc.c:1646  */
    break;

  case 288:
#line 1358 "gram.y" /* yacc.c:1646  */
    { setisol->lines[(int) (yyvsp[-2].val)] = (int) (yyvsp[0].val); }
#line 6514 "y.tab.c" /* yacc.c:1646  */
    break;

  case 289:
#line 1359 "gram.y" /* yacc.c:1646  */
    { setprops = &(setisol->p); }
#line 6520 "y.tab.c" /* yacc.c:1646  */
    break;

  case 291:
#line 1360 "gram.y" /* yacc.c:1646  */
    { setprops = &(setisol->p); }
#line 6526 "y.tab.c" /* yacc.c:1646  */
    break;

  case 293:
#line 1361 "gram.y" /* yacc.c:1646  */
    {
            setisol->xlen = (yyvsp[-2].val);
            setisol->ylen = (yyvsp[0].val);
        }
#line 6535 "y.tab.c" /* yacc.c:1646  */
    break;

  case 294:
#line 1365 "gram.y" /* yacc.c:1646  */
    {
            setisol->xgap = (yyvsp[-2].val);
            setisol->ygap = (yyvsp[0].val);
        }
#line 6544 "y.tab.c" /* yacc.c:1646  */
    break;

  case 295:
#line 1372 "gram.y" /* yacc.c:1646  */
    { sethistbox = &g[curg].hbox[(int) (yyvsp[0].val)]; }
#line 6550 "y.tab.c" /* yacc.c:1646  */
    break;

  case 296:
#line 1373 "gram.y" /* yacc.c:1646  */
    { setprops = &(sethistbox->p); }
#line 6556 "y.tab.c" /* yacc.c:1646  */
    break;

  case 298:
#line 1375 "gram.y" /* yacc.c:1646  */
    {
		sethistbox->precx = (int) (yyvsp[-2].val);
		sethistbox->precy = (int) (yyvsp[0].val);
	}
#line 6565 "y.tab.c" /* yacc.c:1646  */
    break;

  case 299:
#line 1379 "gram.y" /* yacc.c:1646  */
    { sethistbox->attach = (int) (yyvsp[0].val); }
#line 6571 "y.tab.c" /* yacc.c:1646  */
    break;

  case 300:
#line 1380 "gram.y" /* yacc.c:1646  */
    { sethistbox->xtickm = (double ) (yyvsp[-2].val); sethistbox->ytickm = (double ) (yyvsp[0].val);}
#line 6577 "y.tab.c" /* yacc.c:1646  */
    break;

  case 301:
#line 1381 "gram.y" /* yacc.c:1646  */
    { sethistbox->loctype = (int) (yyvsp[0].pset); }
#line 6583 "y.tab.c" /* yacc.c:1646  */
    break;

  case 302:
#line 1382 "gram.y" /* yacc.c:1646  */
    { sethistbox->display_marker = (int) (yyvsp[0].pset); }
#line 6589 "y.tab.c" /* yacc.c:1646  */
    break;

  case 303:
#line 1383 "gram.y" /* yacc.c:1646  */
    { sethistbox->display = (int) (yyvsp[0].pset); }
#line 6595 "y.tab.c" /* yacc.c:1646  */
    break;

  case 304:
#line 1384 "gram.y" /* yacc.c:1646  */
    { sethistbox->adcirc[(int) (yyvsp[-1].val)] = (int) (yyvsp[0].pset) == TRUEP; }
#line 6601 "y.tab.c" /* yacc.c:1646  */
    break;

  case 305:
#line 1385 "gram.y" /* yacc.c:1646  */
    { sethistbox->ap[(int) (yyvsp[-2].val)].color = (int) (yyvsp[0].val); }
#line 6607 "y.tab.c" /* yacc.c:1646  */
    break;

  case 306:
#line 1386 "gram.y" /* yacc.c:1646  */
    { sethistbox->p.color = (yyvsp[0].val); }
#line 6613 "y.tab.c" /* yacc.c:1646  */
    break;

  case 307:
#line 1387 "gram.y" /* yacc.c:1646  */
    { sethistbox->p.linew = (yyvsp[0].val); }
#line 6619 "y.tab.c" /* yacc.c:1646  */
    break;

  case 308:
#line 1388 "gram.y" /* yacc.c:1646  */
    { sethistbox->p.fillcol = (yyvsp[0].val); }
#line 6625 "y.tab.c" /* yacc.c:1646  */
    break;

  case 309:
#line 1389 "gram.y" /* yacc.c:1646  */
    { read_hist((int) (yyvsp[-1].val), TIME, (char *) (yyvsp[0].str)); }
#line 6631 "y.tab.c" /* yacc.c:1646  */
    break;

  case 310:
#line 1390 "gram.y" /* yacc.c:1646  */
    { sethistbox->thist = (int) (yyvsp[0].pset) == TRUEP; }
#line 6637 "y.tab.c" /* yacc.c:1646  */
    break;

  case 311:
#line 1391 "gram.y" /* yacc.c:1646  */
    { sethistbox->hp.color = (int) (yyvsp[0].val); }
#line 6643 "y.tab.c" /* yacc.c:1646  */
    break;

  case 312:
#line 1393 "gram.y" /* yacc.c:1646  */
    { 
		sethistbox->wx1 = (double) (yyvsp[-6].val); 
		sethistbox->wy1 = (double) (yyvsp[-4].val); 
		sethistbox->wx2 = (double) (yyvsp[-2].val); 
		sethistbox->wy2 = (double) (yyvsp[0].val); 
	}
#line 6654 "y.tab.c" /* yacc.c:1646  */
    break;

  case 313:
#line 1400 "gram.y" /* yacc.c:1646  */
    { 
		sethistbox->vx = (double) (yyvsp[-2].val); 
		sethistbox->vy = (double) (yyvsp[0].val); 
	}
#line 6663 "y.tab.c" /* yacc.c:1646  */
    break;

  case 314:
#line 1405 "gram.y" /* yacc.c:1646  */
    { 
		sethistbox->locx = (double) (yyvsp[-2].val); 
		sethistbox->locy = (double) (yyvsp[0].val); 
	}
#line 6672 "y.tab.c" /* yacc.c:1646  */
    break;

  case 315:
#line 1410 "gram.y" /* yacc.c:1646  */
    { 
		sethistbox->x = (double) (yyvsp[-2].val); 
		sethistbox->y = (double) (yyvsp[0].val); 
	}
#line 6681 "y.tab.c" /* yacc.c:1646  */
    break;

  case 316:
#line 1417 "gram.y" /* yacc.c:1646  */
    { setzoombox = &g[curg].zbox[(int) (yyvsp[0].val)]; }
#line 6687 "y.tab.c" /* yacc.c:1646  */
    break;

  case 317:
#line 1418 "gram.y" /* yacc.c:1646  */
    { setprops = &(setzoombox->p); }
#line 6693 "y.tab.c" /* yacc.c:1646  */
    break;

  case 319:
#line 1420 "gram.y" /* yacc.c:1646  */
    {
		setzoombox->precx = (int) (yyvsp[-2].val);
		setzoombox->precy = (int) (yyvsp[0].val);
	}
#line 6702 "y.tab.c" /* yacc.c:1646  */
    break;

  case 320:
#line 1424 "gram.y" /* yacc.c:1646  */
    { setzoombox->attach = (int) (yyvsp[0].val); }
#line 6708 "y.tab.c" /* yacc.c:1646  */
    break;

  case 321:
#line 1425 "gram.y" /* yacc.c:1646  */
    { setzoombox->loctype = (int) (yyvsp[0].pset); }
#line 6714 "y.tab.c" /* yacc.c:1646  */
    break;

  case 322:
#line 1426 "gram.y" /* yacc.c:1646  */
    { setzoombox->display_marker = (int) (yyvsp[0].pset); }
#line 6720 "y.tab.c" /* yacc.c:1646  */
    break;

  case 323:
#line 1427 "gram.y" /* yacc.c:1646  */
    { setzoombox->active = (int) (yyvsp[0].pset); }
#line 6726 "y.tab.c" /* yacc.c:1646  */
    break;

  case 324:
#line 1428 "gram.y" /* yacc.c:1646  */
    { setzoombox->display = (int) (yyvsp[0].pset); }
#line 6732 "y.tab.c" /* yacc.c:1646  */
    break;

  case 325:
#line 1429 "gram.y" /* yacc.c:1646  */
    { setzoombox->expand = (int) (yyvsp[0].val); }
#line 6738 "y.tab.c" /* yacc.c:1646  */
    break;

  case 326:
#line 1430 "gram.y" /* yacc.c:1646  */
    { setzoombox->expand = (int) (yyvsp[0].val); }
#line 6744 "y.tab.c" /* yacc.c:1646  */
    break;

  case 327:
#line 1431 "gram.y" /* yacc.c:1646  */
    { setzoombox->p.color = (yyvsp[0].val); }
#line 6750 "y.tab.c" /* yacc.c:1646  */
    break;

  case 328:
#line 1432 "gram.y" /* yacc.c:1646  */
    { setzoombox->p.linew = (yyvsp[0].val); }
#line 6756 "y.tab.c" /* yacc.c:1646  */
    break;

  case 329:
#line 1433 "gram.y" /* yacc.c:1646  */
    { setzoombox->p.fillcol = (yyvsp[0].val); }
#line 6762 "y.tab.c" /* yacc.c:1646  */
    break;

  case 330:
#line 1435 "gram.y" /* yacc.c:1646  */
    { 
		setzoombox->wx1 = (double) (yyvsp[-6].val); 
		setzoombox->wy1 = (double) (yyvsp[-4].val); 
		setzoombox->wx2 = (double) (yyvsp[-2].val); 
		setzoombox->wy2 = (double) (yyvsp[0].val); 
	}
#line 6773 "y.tab.c" /* yacc.c:1646  */
    break;

  case 331:
#line 1442 "gram.y" /* yacc.c:1646  */
    { 
		setzoombox->vx = (double) (yyvsp[-2].val); 
		setzoombox->vy = (double) (yyvsp[0].val); 
	}
#line 6782 "y.tab.c" /* yacc.c:1646  */
    break;

  case 332:
#line 1447 "gram.y" /* yacc.c:1646  */
    { 
		setzoombox->locx = (double) (yyvsp[-2].val); 
		setzoombox->locy = (double) (yyvsp[0].val); 
	}
#line 6791 "y.tab.c" /* yacc.c:1646  */
    break;

  case 333:
#line 1452 "gram.y" /* yacc.c:1646  */
    { 
		setzoombox->x = (double) (yyvsp[-2].val); 
		setzoombox->y = (double) (yyvsp[0].val); 
	}
#line 6800 "y.tab.c" /* yacc.c:1646  */
    break;

  case 334:
#line 1459 "gram.y" /* yacc.c:1646  */
    { g[curg].wl.active = (yyvsp[0].pset); }
#line 6806 "y.tab.c" /* yacc.c:1646  */
    break;

  case 335:
#line 1460 "gram.y" /* yacc.c:1646  */
    { g[curg].wl.len = (yyvsp[0].val); }
#line 6812 "y.tab.c" /* yacc.c:1646  */
    break;

  case 336:
#line 1461 "gram.y" /* yacc.c:1646  */
    { g[curg].wl.scale = (yyvsp[0].val); }
#line 6818 "y.tab.c" /* yacc.c:1646  */
    break;

  case 337:
#line 1462 "gram.y" /* yacc.c:1646  */
    { g[curg].wl.p.color = (yyvsp[0].val); }
#line 6824 "y.tab.c" /* yacc.c:1646  */
    break;

  case 338:
#line 1463 "gram.y" /* yacc.c:1646  */
    { g[curg].wl.loctype = (yyvsp[0].pset); }
#line 6830 "y.tab.c" /* yacc.c:1646  */
    break;

  case 339:
#line 1465 "gram.y" /* yacc.c:1646  */
    { 
		g[curg].wl.x = (yyvsp[-2].val);
		g[curg].wl.y = (yyvsp[0].val);
	}
#line 6839 "y.tab.c" /* yacc.c:1646  */
    break;

  case 340:
#line 1470 "gram.y" /* yacc.c:1646  */
    { 
		g[curg].wl.units = (yyvsp[0].pset);
		switch (g[curg].wl.units) {
		case MM: g[curg].wl.unitfac = 0.001; break;
		case CM: g[curg].wl.unitfac = 0.01; break;
		case M: g[curg].wl.unitfac = 1.0; break;
		case KM: g[curg].wl.unitfac = 1000.0; break;
		default: fprintf(stderr, "Unknown units for velocity scale\n"); break;
		}
	}
#line 6854 "y.tab.c" /* yacc.c:1646  */
    break;

  case 341:
#line 1480 "gram.y" /* yacc.c:1646  */
    { }
#line 6860 "y.tab.c" /* yacc.c:1646  */
    break;

  case 342:
#line 1485 "gram.y" /* yacc.c:1646  */
    { g[curg].vl.active = (yyvsp[0].pset); }
#line 6866 "y.tab.c" /* yacc.c:1646  */
    break;

  case 343:
#line 1486 "gram.y" /* yacc.c:1646  */
    { g[curg].vl.len = (yyvsp[0].val); }
#line 6872 "y.tab.c" /* yacc.c:1646  */
    break;

  case 344:
#line 1487 "gram.y" /* yacc.c:1646  */
    { g[curg].vl.scale = (yyvsp[0].val); }
#line 6878 "y.tab.c" /* yacc.c:1646  */
    break;

  case 345:
#line 1488 "gram.y" /* yacc.c:1646  */
    { g[curg].vl.p.color = (yyvsp[0].val); }
#line 6884 "y.tab.c" /* yacc.c:1646  */
    break;

  case 346:
#line 1489 "gram.y" /* yacc.c:1646  */
    { g[curg].vl.loctype = (yyvsp[0].pset); }
#line 6890 "y.tab.c" /* yacc.c:1646  */
    break;

  case 347:
#line 1491 "gram.y" /* yacc.c:1646  */
    { 
		g[curg].vl.x = (yyvsp[-2].val);
		g[curg].vl.y = (yyvsp[0].val);
	}
#line 6899 "y.tab.c" /* yacc.c:1646  */
    break;

  case 348:
#line 1496 "gram.y" /* yacc.c:1646  */
    { 
		g[curg].vl.units = (yyvsp[0].pset);
		switch (g[curg].vl.units) {
		case MM: g[curg].vl.unitfac = 0.001; break;
		case CM: g[curg].vl.unitfac = 0.01; break;
		case M: g[curg].vl.unitfac = 1.0; break;
		case KM: g[curg].vl.unitfac = 1000.0; break;
		default: fprintf(stderr, "Unknown units for velocity scale\n"); break;
		}
	}
#line 6914 "y.tab.c" /* yacc.c:1646  */
    break;

  case 349:
#line 1506 "gram.y" /* yacc.c:1646  */
    { }
#line 6920 "y.tab.c" /* yacc.c:1646  */
    break;

  case 350:
#line 1510 "gram.y" /* yacc.c:1646  */
    { g[curg].mapscale.active = (yyvsp[0].pset); }
#line 6926 "y.tab.c" /* yacc.c:1646  */
    break;

  case 351:
#line 1511 "gram.y" /* yacc.c:1646  */
    { g[curg].mapscale.len = (yyvsp[0].val); }
#line 6932 "y.tab.c" /* yacc.c:1646  */
    break;

  case 352:
#line 1512 "gram.y" /* yacc.c:1646  */
    { g[curg].mapscale.p.color = (yyvsp[0].val); }
#line 6938 "y.tab.c" /* yacc.c:1646  */
    break;

  case 353:
#line 1513 "gram.y" /* yacc.c:1646  */
    { g[curg].mapscale.scale = (yyvsp[0].val); }
#line 6944 "y.tab.c" /* yacc.c:1646  */
    break;

  case 354:
#line 1514 "gram.y" /* yacc.c:1646  */
    { g[curg].mapscale.loctype = (yyvsp[0].pset); }
#line 6950 "y.tab.c" /* yacc.c:1646  */
    break;

  case 355:
#line 1516 "gram.y" /* yacc.c:1646  */
    { 
		g[curg].mapscale.x = (yyvsp[-2].val);
		g[curg].mapscale.y = (yyvsp[0].val);
	}
#line 6959 "y.tab.c" /* yacc.c:1646  */
    break;

  case 356:
#line 1521 "gram.y" /* yacc.c:1646  */
    { 
		g[curg].mapscale.units = (yyvsp[0].pset);
		switch (g[curg].mapscale.units) {
		case MM: g[curg].mapscale.unitfac = 0.001; break;
		case CM: g[curg].mapscale.unitfac = 0.01; break;
		case M: g[curg].mapscale.unitfac = 1.0; break;
		case KM: g[curg].mapscale.unitfac = 1000.0; break;
		default: fprintf(stderr, "Unknown units for mapscape scale\n"); break;
		}
	}
#line 6974 "y.tab.c" /* yacc.c:1646  */
    break;

  case 357:
#line 1531 "gram.y" /* yacc.c:1646  */
    { }
#line 6980 "y.tab.c" /* yacc.c:1646  */
    break;

  case 358:
#line 1535 "gram.y" /* yacc.c:1646  */
    { g[curg].tidalclock.active = (yyvsp[0].pset); }
#line 6986 "y.tab.c" /* yacc.c:1646  */
    break;

  case 359:
#line 1536 "gram.y" /* yacc.c:1646  */
    { g[curg].tidalclock.p.color = (yyvsp[0].val); }
#line 6992 "y.tab.c" /* yacc.c:1646  */
    break;

  case 360:
#line 1537 "gram.y" /* yacc.c:1646  */
    { g[curg].tidalclock.p.fillcol = (yyvsp[0].val); }
#line 6998 "y.tab.c" /* yacc.c:1646  */
    break;

  case 361:
#line 1538 "gram.y" /* yacc.c:1646  */
    { g[curg].tidalclock.total_time = (yyvsp[0].val); }
#line 7004 "y.tab.c" /* yacc.c:1646  */
    break;

  case 362:
#line 1539 "gram.y" /* yacc.c:1646  */
    { g[curg].tidalclock.loctype = (yyvsp[0].pset); }
#line 7010 "y.tab.c" /* yacc.c:1646  */
    break;

  case 363:
#line 1541 "gram.y" /* yacc.c:1646  */
    { 
		g[curg].tidalclock.x = (yyvsp[-2].val);
		g[curg].tidalclock.y = (yyvsp[0].val);
	}
#line 7019 "y.tab.c" /* yacc.c:1646  */
    break;

  case 364:
#line 1545 "gram.y" /* yacc.c:1646  */
    { }
#line 7025 "y.tab.c" /* yacc.c:1646  */
    break;

  case 365:
#line 1549 "gram.y" /* yacc.c:1646  */
    { g[curg].timeinfo.active = (yyvsp[0].pset); }
#line 7031 "y.tab.c" /* yacc.c:1646  */
    break;

  case 366:
#line 1551 "gram.y" /* yacc.c:1646  */
    { 
		strcpy(g[curg].timeinfo.start, (char *) (yyvsp[0].str)); 
		time_info_start(curg);
	}
#line 7040 "y.tab.c" /* yacc.c:1646  */
    break;

  case 367:
#line 1556 "gram.y" /* yacc.c:1646  */
    {
	    g[curg].timeinfo.x = (yyvsp[-2].val);
	    g[curg].timeinfo.y = (yyvsp[0].val);
	}
#line 7049 "y.tab.c" /* yacc.c:1646  */
    break;

  case 368:
#line 1560 "gram.y" /* yacc.c:1646  */
    { g[curg].timeinfo.loctype = (yyvsp[0].pset); }
#line 7055 "y.tab.c" /* yacc.c:1646  */
    break;

  case 369:
#line 1561 "gram.y" /* yacc.c:1646  */
    { g[curg].timeinfo.linew = (int) (yyvsp[0].val); }
#line 7061 "y.tab.c" /* yacc.c:1646  */
    break;

  case 370:
#line 1562 "gram.y" /* yacc.c:1646  */
    { g[curg].timeinfo.color = (int) (yyvsp[0].val); }
#line 7067 "y.tab.c" /* yacc.c:1646  */
    break;

  case 371:
#line 1563 "gram.y" /* yacc.c:1646  */
    { g[curg].timeinfo.rot = (int) (yyvsp[0].val); }
#line 7073 "y.tab.c" /* yacc.c:1646  */
    break;

  case 372:
#line 1564 "gram.y" /* yacc.c:1646  */
    { g[curg].timeinfo.font = (int) (yyvsp[0].val); }
#line 7079 "y.tab.c" /* yacc.c:1646  */
    break;

  case 373:
#line 1565 "gram.y" /* yacc.c:1646  */
    { g[curg].timeinfo.just = (int) (yyvsp[0].val); }
#line 7085 "y.tab.c" /* yacc.c:1646  */
    break;

  case 374:
#line 1566 "gram.y" /* yacc.c:1646  */
    { g[curg].timeinfo.charsize = (double) (yyvsp[0].val); }
#line 7091 "y.tab.c" /* yacc.c:1646  */
    break;

  case 375:
#line 1567 "gram.y" /* yacc.c:1646  */
    { strcpy(g[curg].timeinfo.format, (char *) (yyvsp[0].str)); }
#line 7097 "y.tab.c" /* yacc.c:1646  */
    break;

  case 376:
#line 1571 "gram.y" /* yacc.c:1646  */
    { g[curg].timeline.active = (int) (yyvsp[0].pset); }
#line 7103 "y.tab.c" /* yacc.c:1646  */
    break;

  case 377:
#line 1572 "gram.y" /* yacc.c:1646  */
    { g[curg].timeline.len = (int) (yyvsp[0].val); }
#line 7109 "y.tab.c" /* yacc.c:1646  */
    break;

  case 378:
#line 1573 "gram.y" /* yacc.c:1646  */
    { g[curg].timeline.width = (int) (yyvsp[0].val); }
#line 7115 "y.tab.c" /* yacc.c:1646  */
    break;

  case 379:
#line 1574 "gram.y" /* yacc.c:1646  */
    { g[curg].timeline.start = (double) (yyvsp[0].val); }
#line 7121 "y.tab.c" /* yacc.c:1646  */
    break;

  case 380:
#line 1575 "gram.y" /* yacc.c:1646  */
    { g[curg].timeline.stop = (double) (yyvsp[0].val); }
#line 7127 "y.tab.c" /* yacc.c:1646  */
    break;

  case 381:
#line 1576 "gram.y" /* yacc.c:1646  */
    { g[curg].timeline.step = (double) (yyvsp[0].val); }
#line 7133 "y.tab.c" /* yacc.c:1646  */
    break;

  case 382:
#line 1577 "gram.y" /* yacc.c:1646  */
    { g[curg].timeline.p.prec = (int) (yyvsp[0].val); }
#line 7139 "y.tab.c" /* yacc.c:1646  */
    break;

  case 383:
#line 1578 "gram.y" /* yacc.c:1646  */
    { g[curg].timeline.units = (int) (yyvsp[0].val); }
#line 7145 "y.tab.c" /* yacc.c:1646  */
    break;

  case 384:
#line 1579 "gram.y" /* yacc.c:1646  */
    { g[curg].timeline.c1 = g[curg].timeline.c3 = (int) (yyvsp[0].val); }
#line 7151 "y.tab.c" /* yacc.c:1646  */
    break;

  case 385:
#line 1580 "gram.y" /* yacc.c:1646  */
    { g[curg].timeline.c2 = (int) (yyvsp[0].val); }
#line 7157 "y.tab.c" /* yacc.c:1646  */
    break;

  case 386:
#line 1581 "gram.y" /* yacc.c:1646  */
    { g[curg].timeline.loctype = (int) (yyvsp[0].pset); }
#line 7163 "y.tab.c" /* yacc.c:1646  */
    break;

  case 387:
#line 1583 "gram.y" /* yacc.c:1646  */
    { 
		g[curg].timeline.x = (double) (yyvsp[-2].val);
		g[curg].timeline.y = (double) (yyvsp[0].val);
	}
#line 7172 "y.tab.c" /* yacc.c:1646  */
    break;

  case 388:
#line 1587 "gram.y" /* yacc.c:1646  */
    { }
#line 7178 "y.tab.c" /* yacc.c:1646  */
    break;

  case 389:
#line 1591 "gram.y" /* yacc.c:1646  */
    { setslice = &g[curg].sbox[(int) (yyvsp[0].val)]; }
#line 7184 "y.tab.c" /* yacc.c:1646  */
    break;

  case 390:
#line 1592 "gram.y" /* yacc.c:1646  */
    { setprops = &(setslice->p); }
#line 7190 "y.tab.c" /* yacc.c:1646  */
    break;

  case 392:
#line 1594 "gram.y" /* yacc.c:1646  */
    {
		setslice->precx = (int) (yyvsp[-2].val);
		setslice->precy = (int) (yyvsp[0].val);
	}
#line 7199 "y.tab.c" /* yacc.c:1646  */
    break;

  case 393:
#line 1598 "gram.y" /* yacc.c:1646  */
    { setslice->attach = (int) (yyvsp[0].val); }
#line 7205 "y.tab.c" /* yacc.c:1646  */
    break;

  case 394:
#line 1599 "gram.y" /* yacc.c:1646  */
    { setslice->loctype = (int) (yyvsp[0].pset); }
#line 7211 "y.tab.c" /* yacc.c:1646  */
    break;

  case 395:
#line 1600 "gram.y" /* yacc.c:1646  */
    { setslice->display_marker = (int) (yyvsp[0].pset); }
#line 7217 "y.tab.c" /* yacc.c:1646  */
    break;

  case 396:
#line 1601 "gram.y" /* yacc.c:1646  */
    { setslice->active = (int) (yyvsp[0].pset); }
#line 7223 "y.tab.c" /* yacc.c:1646  */
    break;

  case 397:
#line 1602 "gram.y" /* yacc.c:1646  */
    { setslice->display = (int) (yyvsp[0].pset); }
#line 7229 "y.tab.c" /* yacc.c:1646  */
    break;

  case 398:
#line 1603 "gram.y" /* yacc.c:1646  */
    { setslice->p.color = (yyvsp[0].val); }
#line 7235 "y.tab.c" /* yacc.c:1646  */
    break;

  case 399:
#line 1604 "gram.y" /* yacc.c:1646  */
    { setslice->p.linew = (yyvsp[0].val); }
#line 7241 "y.tab.c" /* yacc.c:1646  */
    break;

  case 400:
#line 1605 "gram.y" /* yacc.c:1646  */
    { setslice->p.fillcol = (yyvsp[0].val); }
#line 7247 "y.tab.c" /* yacc.c:1646  */
    break;

  case 401:
#line 1607 "gram.y" /* yacc.c:1646  */
    { 
		setslice->wx1 = (double) (yyvsp[-6].val); 
		setslice->wy1 = (double) (yyvsp[-4].val); 
		setslice->wx2 = (double) (yyvsp[-2].val); 
		setslice->wy2 = (double) (yyvsp[0].val); 
	}
#line 7258 "y.tab.c" /* yacc.c:1646  */
    break;

  case 402:
#line 1614 "gram.y" /* yacc.c:1646  */
    { 
		setslice->vx = (double) (yyvsp[-2].val); 
		setslice->vy = (double) (yyvsp[0].val); 
	}
#line 7267 "y.tab.c" /* yacc.c:1646  */
    break;

  case 403:
#line 1619 "gram.y" /* yacc.c:1646  */
    { 
		setslice->locx = (double) (yyvsp[-2].val); 
		setslice->locy = (double) (yyvsp[0].val); 
	}
#line 7276 "y.tab.c" /* yacc.c:1646  */
    break;

  case 404:
#line 1624 "gram.y" /* yacc.c:1646  */
    { 
		setslice->x = (double) (yyvsp[-2].val); 
		setslice->y = (double) (yyvsp[0].val); 
	}
#line 7285 "y.tab.c" /* yacc.c:1646  */
    break;

  case 405:
#line 1631 "gram.y" /* yacc.c:1646  */
    { setprops->color = (yyvsp[0].val); }
#line 7291 "y.tab.c" /* yacc.c:1646  */
    break;

  case 406:
#line 1632 "gram.y" /* yacc.c:1646  */
    { setprops->linew = (yyvsp[0].val); }
#line 7297 "y.tab.c" /* yacc.c:1646  */
    break;

  case 407:
#line 1633 "gram.y" /* yacc.c:1646  */
    { setprops->lines = (yyvsp[0].val); }
#line 7303 "y.tab.c" /* yacc.c:1646  */
    break;

  case 408:
#line 1634 "gram.y" /* yacc.c:1646  */
    { setprops->format = (yyvsp[0].pset); }
#line 7309 "y.tab.c" /* yacc.c:1646  */
    break;

  case 409:
#line 1635 "gram.y" /* yacc.c:1646  */
    { setprops->font = (yyvsp[0].val); }
#line 7315 "y.tab.c" /* yacc.c:1646  */
    break;

  case 410:
#line 1636 "gram.y" /* yacc.c:1646  */
    { setprops->prec = (yyvsp[0].val); }
#line 7321 "y.tab.c" /* yacc.c:1646  */
    break;

  case 411:
#line 1637 "gram.y" /* yacc.c:1646  */
    { setprops->charsize = (yyvsp[0].val); }
#line 7327 "y.tab.c" /* yacc.c:1646  */
    break;

  case 412:
#line 1638 "gram.y" /* yacc.c:1646  */
    { setprops->symbol = (yyvsp[0].val); }
#line 7333 "y.tab.c" /* yacc.c:1646  */
    break;

  case 413:
#line 1639 "gram.y" /* yacc.c:1646  */
    { setprops->symsize = (yyvsp[0].val); }
#line 7339 "y.tab.c" /* yacc.c:1646  */
    break;

  case 414:
#line 1640 "gram.y" /* yacc.c:1646  */
    { setprops->fill = (yyvsp[0].pset); }
#line 7345 "y.tab.c" /* yacc.c:1646  */
    break;

  case 415:
#line 1641 "gram.y" /* yacc.c:1646  */
    { setprops->fillusing = (yyvsp[0].pset); }
#line 7351 "y.tab.c" /* yacc.c:1646  */
    break;

  case 416:
#line 1642 "gram.y" /* yacc.c:1646  */
    { setprops->fillcol = (yyvsp[0].val); }
#line 7357 "y.tab.c" /* yacc.c:1646  */
    break;

  case 417:
#line 1643 "gram.y" /* yacc.c:1646  */
    { setprops->fillpat = (yyvsp[0].val); }
#line 7363 "y.tab.c" /* yacc.c:1646  */
    break;

  case 418:
#line 1644 "gram.y" /* yacc.c:1646  */
    { setprops->arrow = (yyvsp[0].val); }
#line 7369 "y.tab.c" /* yacc.c:1646  */
    break;

  case 419:
#line 1645 "gram.y" /* yacc.c:1646  */
    { setprops->atype = (yyvsp[0].val); }
#line 7375 "y.tab.c" /* yacc.c:1646  */
    break;

  case 420:
#line 1646 "gram.y" /* yacc.c:1646  */
    { setprops->asize = (yyvsp[0].val); }
#line 7381 "y.tab.c" /* yacc.c:1646  */
    break;

  case 421:
#line 1650 "gram.y" /* yacc.c:1646  */
    { curg = (int) (yyvsp[0].pset); }
#line 7387 "y.tab.c" /* yacc.c:1646  */
    break;

  case 422:
#line 1651 "gram.y" /* yacc.c:1646  */
    { curg = (int) (yyvsp[0].val); }
#line 7393 "y.tab.c" /* yacc.c:1646  */
    break;

  case 423:
#line 1652 "gram.y" /* yacc.c:1646  */
    { kill_graph((yyvsp[0].pset)); }
#line 7399 "y.tab.c" /* yacc.c:1646  */
    break;

  case 424:
#line 1653 "gram.y" /* yacc.c:1646  */
    { kill_graph(maxgraph); }
#line 7405 "y.tab.c" /* yacc.c:1646  */
    break;

  case 425:
#line 1655 "gram.y" /* yacc.c:1646  */
    {
	    extern int go_locateflag;
	    go_locateflag = ((yyvsp[0].pset) == ON);
	}
#line 7414 "y.tab.c" /* yacc.c:1646  */
    break;

  case 426:
#line 1660 "gram.y" /* yacc.c:1646  */
    {
	    cg = curg = (int) (yyvsp[0].pset);
	    draw_focus(curg);
	    defineworld(g[curg].w.xg1, g[curg].w.yg1, g[curg].w.xg2, g[curg].w.yg2, 
			islogx(curg), islogy(curg));
	    viewport(g[curg].v.xv1, g[curg].v.yv1, g[curg].v.xv2, g[curg].v.yv2);
	    draw_focus(curg);
	    update_all(curg);
	}
#line 7428 "y.tab.c" /* yacc.c:1646  */
    break;

  case 427:
#line 1669 "gram.y" /* yacc.c:1646  */
    { draw_focus_flag = (yyvsp[0].pset); }
#line 7434 "y.tab.c" /* yacc.c:1646  */
    break;

  case 428:
#line 1670 "gram.y" /* yacc.c:1646  */
    { focus_policy = (yyvsp[0].pset); }
#line 7440 "y.tab.c" /* yacc.c:1646  */
    break;

  case 429:
#line 1671 "gram.y" /* yacc.c:1646  */
    { focus_policy = (yyvsp[0].pset); }
#line 7446 "y.tab.c" /* yacc.c:1646  */
    break;

  case 430:
#line 1672 "gram.y" /* yacc.c:1646  */
    { focus_policy = (yyvsp[0].pset); }
#line 7452 "y.tab.c" /* yacc.c:1646  */
    break;

  case 431:
#line 1673 "gram.y" /* yacc.c:1646  */
    { cursource = (yyvsp[0].pset); }
#line 7458 "y.tab.c" /* yacc.c:1646  */
    break;

  case 432:
#line 1674 "gram.y" /* yacc.c:1646  */
    { push_world(); }
#line 7464 "y.tab.c" /* yacc.c:1646  */
    break;

  case 433:
#line 1675 "gram.y" /* yacc.c:1646  */
    { pop_world(); }
#line 7470 "y.tab.c" /* yacc.c:1646  */
    break;

  case 434:
#line 1676 "gram.y" /* yacc.c:1646  */
    { cycle_world_stack(); }
#line 7476 "y.tab.c" /* yacc.c:1646  */
    break;

  case 435:
#line 1677 "gram.y" /* yacc.c:1646  */
    {
	    if ((int) (yyvsp[0].val) > 0)
		show_world_stack((int) (yyvsp[0].val) - 1);
	}
#line 7485 "y.tab.c" /* yacc.c:1646  */
    break;

  case 436:
#line 1682 "gram.y" /* yacc.c:1646  */
    {
	    add_world(curg, (yyvsp[-14].val), (yyvsp[-12].val), (yyvsp[-10].val), (yyvsp[-8].val), (yyvsp[-6].val), (yyvsp[-4].val), (yyvsp[-2].val), (yyvsp[0].val));
	}
#line 7493 "y.tab.c" /* yacc.c:1646  */
    break;

  case 437:
#line 1685 "gram.y" /* yacc.c:1646  */
    { clear_world_stack(); }
#line 7499 "y.tab.c" /* yacc.c:1646  */
    break;

  case 438:
#line 1687 "gram.y" /* yacc.c:1646  */
    {
	    g[curg].w.xg1 = (yyvsp[-6].val);
	    g[curg].w.yg1 = (yyvsp[-4].val);
	    g[curg].w.xg2 = (yyvsp[-2].val);
	    g[curg].w.yg2 = (yyvsp[0].val);
	}
#line 7510 "y.tab.c" /* yacc.c:1646  */
    break;

  case 439:
#line 1693 "gram.y" /* yacc.c:1646  */
    { g[curg].w.xg1 = (yyvsp[0].val); }
#line 7516 "y.tab.c" /* yacc.c:1646  */
    break;

  case 440:
#line 1694 "gram.y" /* yacc.c:1646  */
    { g[curg].w.xg2 = (yyvsp[0].val); }
#line 7522 "y.tab.c" /* yacc.c:1646  */
    break;

  case 441:
#line 1695 "gram.y" /* yacc.c:1646  */
    { g[curg].w.yg1 = (yyvsp[0].val); }
#line 7528 "y.tab.c" /* yacc.c:1646  */
    break;

  case 442:
#line 1696 "gram.y" /* yacc.c:1646  */
    { g[curg].w.yg2 = (yyvsp[0].val); }
#line 7534 "y.tab.c" /* yacc.c:1646  */
    break;

  case 443:
#line 1698 "gram.y" /* yacc.c:1646  */
    {
	    g[curg].v.xv1 = (yyvsp[-6].val);
	    g[curg].v.yv1 = (yyvsp[-4].val);
	    g[curg].v.xv2 = (yyvsp[-2].val);
	    g[curg].v.yv2 = (yyvsp[0].val);
	}
#line 7545 "y.tab.c" /* yacc.c:1646  */
    break;

  case 444:
#line 1704 "gram.y" /* yacc.c:1646  */
    { g[curg].v.xv1 = (yyvsp[0].val); }
#line 7551 "y.tab.c" /* yacc.c:1646  */
    break;

  case 445:
#line 1705 "gram.y" /* yacc.c:1646  */
    { g[curg].v.xv2 = (yyvsp[0].val); }
#line 7557 "y.tab.c" /* yacc.c:1646  */
    break;

  case 446:
#line 1706 "gram.y" /* yacc.c:1646  */
    { g[curg].v.yv1 = (yyvsp[0].val); }
#line 7563 "y.tab.c" /* yacc.c:1646  */
    break;

  case 447:
#line 1707 "gram.y" /* yacc.c:1646  */
    { g[curg].v.yv2 = (yyvsp[0].val); }
#line 7569 "y.tab.c" /* yacc.c:1646  */
    break;

  case 448:
#line 1708 "gram.y" /* yacc.c:1646  */
    {
	    set_plotstr_string(&g[curg].labs.title, (char *) (yyvsp[0].str));
	}
#line 7577 "y.tab.c" /* yacc.c:1646  */
    break;

  case 449:
#line 1711 "gram.y" /* yacc.c:1646  */
    {
	    g[curg].labs.title.font = checkon(FONTP, g[curg].labs.title.font, (int) (yyvsp[0].val));
	}
#line 7585 "y.tab.c" /* yacc.c:1646  */
    break;

  case 450:
#line 1714 "gram.y" /* yacc.c:1646  */
    { g[curg].labs.title.charsize = (yyvsp[0].val); }
#line 7591 "y.tab.c" /* yacc.c:1646  */
    break;

  case 451:
#line 1715 "gram.y" /* yacc.c:1646  */
    {
	    g[curg].labs.title.color = checkon(COLOR, g[curg].labs.title.color, (int) (yyvsp[0].val));
	}
#line 7599 "y.tab.c" /* yacc.c:1646  */
    break;

  case 452:
#line 1719 "gram.y" /* yacc.c:1646  */
    {
	    g[curg].labs.title.linew = checkon(LINEWIDTH, g[curg].labs.title.linew, (int) (yyvsp[0].val));
	}
#line 7607 "y.tab.c" /* yacc.c:1646  */
    break;

  case 453:
#line 1722 "gram.y" /* yacc.c:1646  */
    {
	    set_plotstr_string(&g[curg].labs.stitle, (char *) (yyvsp[0].str));
	}
#line 7615 "y.tab.c" /* yacc.c:1646  */
    break;

  case 454:
#line 1725 "gram.y" /* yacc.c:1646  */
    {
	    g[curg].labs.stitle.font = checkon(FONTP, g[curg].labs.stitle.font, (int) (yyvsp[0].val));
	}
#line 7623 "y.tab.c" /* yacc.c:1646  */
    break;

  case 455:
#line 1728 "gram.y" /* yacc.c:1646  */
    { g[curg].labs.stitle.charsize = (yyvsp[0].val); }
#line 7629 "y.tab.c" /* yacc.c:1646  */
    break;

  case 456:
#line 1730 "gram.y" /* yacc.c:1646  */
    {
	    g[curg].labs.stitle.color = checkon(COLOR, g[curg].labs.stitle.color, (int) (yyvsp[0].val));
	}
#line 7637 "y.tab.c" /* yacc.c:1646  */
    break;

  case 457:
#line 1734 "gram.y" /* yacc.c:1646  */
    {
	    g[curg].labs.stitle.linew = checkon(LINEWIDTH, g[curg].labs.stitle.color, (int) (yyvsp[0].val));
	}
#line 7645 "y.tab.c" /* yacc.c:1646  */
    break;

  case 458:
#line 1737 "gram.y" /* yacc.c:1646  */
    { g[curg].f.active = (yyvsp[0].pset); }
#line 7651 "y.tab.c" /* yacc.c:1646  */
    break;

  case 459:
#line 1738 "gram.y" /* yacc.c:1646  */
    { g[curg].f.type = (int) (yyvsp[0].val); }
#line 7657 "y.tab.c" /* yacc.c:1646  */
    break;

  case 460:
#line 1739 "gram.y" /* yacc.c:1646  */
    { g[curg].f.lines = checkon(LINESTYLE, g[curg].f.lines, (int) (yyvsp[0].val)); }
#line 7663 "y.tab.c" /* yacc.c:1646  */
    break;

  case 461:
#line 1740 "gram.y" /* yacc.c:1646  */
    { g[curg].f.linew = checkon(LINEWIDTH, g[curg].f.linew, (int) (yyvsp[0].val)); }
#line 7669 "y.tab.c" /* yacc.c:1646  */
    break;

  case 462:
#line 1741 "gram.y" /* yacc.c:1646  */
    { g[curg].f.color = checkon(COLOR, g[curg].f.color, (int) (yyvsp[0].val)); }
#line 7675 "y.tab.c" /* yacc.c:1646  */
    break;

  case 463:
#line 1742 "gram.y" /* yacc.c:1646  */
    { g[curg].f.fillbg = (yyvsp[0].pset); }
#line 7681 "y.tab.c" /* yacc.c:1646  */
    break;

  case 464:
#line 1743 "gram.y" /* yacc.c:1646  */
    { g[curg].f.bgcolor = (int) (yyvsp[0].val); }
#line 7687 "y.tab.c" /* yacc.c:1646  */
    break;

  case 465:
#line 1744 "gram.y" /* yacc.c:1646  */
    { g[(yyvsp[-1].pset)].active = (yyvsp[0].pset); }
#line 7693 "y.tab.c" /* yacc.c:1646  */
    break;

  case 466:
#line 1745 "gram.y" /* yacc.c:1646  */
    { g[(yyvsp[-2].pset)].label = (yyvsp[0].pset); }
#line 7699 "y.tab.c" /* yacc.c:1646  */
    break;

  case 467:
#line 1746 "gram.y" /* yacc.c:1646  */
    { g[(yyvsp[-3].pset)].auto_type = (yyvsp[0].pset); }
#line 7705 "y.tab.c" /* yacc.c:1646  */
    break;

  case 468:
#line 1747 "gram.y" /* yacc.c:1646  */
    { g[(yyvsp[-3].pset)].auto_type = (yyvsp[0].pset); }
#line 7711 "y.tab.c" /* yacc.c:1646  */
    break;

  case 469:
#line 1748 "gram.y" /* yacc.c:1646  */
    { g[(yyvsp[-2].pset)].parmsread = ((yyvsp[0].pset) == FALSEP); }
#line 7717 "y.tab.c" /* yacc.c:1646  */
    break;

  case 470:
#line 1749 "gram.y" /* yacc.c:1646  */
    { g[(yyvsp[-2].pset)].hidden = ((yyvsp[0].pset) == TRUEP); }
#line 7723 "y.tab.c" /* yacc.c:1646  */
    break;

  case 471:
#line 1750 "gram.y" /* yacc.c:1646  */
    { g[(yyvsp[-2].pset)].type = (yyvsp[0].pset); }
#line 7729 "y.tab.c" /* yacc.c:1646  */
    break;

  case 472:
#line 1751 "gram.y" /* yacc.c:1646  */
    { g[(yyvsp[-2].pset)].pointset = ((yyvsp[0].pset) == ON); }
#line 7735 "y.tab.c" /* yacc.c:1646  */
    break;

  case 473:
#line 1753 "gram.y" /* yacc.c:1646  */
    {
	    g[(yyvsp[-4].pset)].fx = (yyvsp[-1].pset);
	    g[(yyvsp[-4].pset)].fy = (yyvsp[0].pset);
	}
#line 7744 "y.tab.c" /* yacc.c:1646  */
    break;

  case 474:
#line 1758 "gram.y" /* yacc.c:1646  */
    {
	    g[(yyvsp[-5].pset)].px = (yyvsp[-2].val);
	    g[(yyvsp[-5].pset)].py = (yyvsp[0].val);
	}
#line 7753 "y.tab.c" /* yacc.c:1646  */
    break;

  case 475:
#line 1763 "gram.y" /* yacc.c:1646  */
    {
	    g[(yyvsp[-5].pset)].dsx = (yyvsp[-2].val);
	    g[(yyvsp[-5].pset)].dsy = (yyvsp[0].val);
	}
#line 7762 "y.tab.c" /* yacc.c:1646  */
    break;

  case 476:
#line 1767 "gram.y" /* yacc.c:1646  */
    { g[(yyvsp[-3].pset)].pt_type = (int) (yyvsp[0].val); }
#line 7768 "y.tab.c" /* yacc.c:1646  */
    break;

  case 479:
#line 1773 "gram.y" /* yacc.c:1646  */
    {}
#line 7774 "y.tab.c" /* yacc.c:1646  */
    break;

  case 480:
#line 1774 "gram.y" /* yacc.c:1646  */
    {}
#line 7780 "y.tab.c" /* yacc.c:1646  */
    break;

  case 481:
#line 1775 "gram.y" /* yacc.c:1646  */
    {}
#line 7786 "y.tab.c" /* yacc.c:1646  */
    break;

  case 482:
#line 1779 "gram.y" /* yacc.c:1646  */
    {}
#line 7792 "y.tab.c" /* yacc.c:1646  */
    break;

  case 483:
#line 1780 "gram.y" /* yacc.c:1646  */
    {}
#line 7798 "y.tab.c" /* yacc.c:1646  */
    break;

  case 484:
#line 1781 "gram.y" /* yacc.c:1646  */
    {}
#line 7804 "y.tab.c" /* yacc.c:1646  */
    break;

  case 485:
#line 1782 "gram.y" /* yacc.c:1646  */
    {}
#line 7810 "y.tab.c" /* yacc.c:1646  */
    break;

  case 486:
#line 1783 "gram.y" /* yacc.c:1646  */
    {}
#line 7816 "y.tab.c" /* yacc.c:1646  */
    break;

  case 487:
#line 1784 "gram.y" /* yacc.c:1646  */
    {}
#line 7822 "y.tab.c" /* yacc.c:1646  */
    break;

  case 488:
#line 1788 "gram.y" /* yacc.c:1646  */
    {}
#line 7828 "y.tab.c" /* yacc.c:1646  */
    break;

  case 489:
#line 1789 "gram.y" /* yacc.c:1646  */
    {}
#line 7834 "y.tab.c" /* yacc.c:1646  */
    break;

  case 490:
#line 1790 "gram.y" /* yacc.c:1646  */
    {}
#line 7840 "y.tab.c" /* yacc.c:1646  */
    break;

  case 491:
#line 1794 "gram.y" /* yacc.c:1646  */
    { set_axis_prop(whichgraph, naxis, (yyvsp[0].pset), 0.0); }
#line 7846 "y.tab.c" /* yacc.c:1646  */
    break;

  case 492:
#line 1795 "gram.y" /* yacc.c:1646  */
    { set_axis_prop(whichgraph, naxis, (yyvsp[-1].pset), (yyvsp[0].val)); }
#line 7852 "y.tab.c" /* yacc.c:1646  */
    break;

  case 493:
#line 1796 "gram.y" /* yacc.c:1646  */
    { set_axis_prop(whichgraph, naxis, (yyvsp[-1].pset), (yyvsp[0].val)); }
#line 7858 "y.tab.c" /* yacc.c:1646  */
    break;

  case 494:
#line 1797 "gram.y" /* yacc.c:1646  */
    { set_axis_prop(whichgraph, naxis, (yyvsp[-1].pset), (yyvsp[0].val)); }
#line 7864 "y.tab.c" /* yacc.c:1646  */
    break;

  case 495:
#line 1798 "gram.y" /* yacc.c:1646  */
    { set_axis_prop(whichgraph, naxis, (yyvsp[-1].pset), (yyvsp[0].val)); }
#line 7870 "y.tab.c" /* yacc.c:1646  */
    break;

  case 496:
#line 1799 "gram.y" /* yacc.c:1646  */
    { set_axis_prop(whichgraph, naxis, (yyvsp[-2].pset), (yyvsp[0].val)); }
#line 7876 "y.tab.c" /* yacc.c:1646  */
    break;

  case 497:
#line 1800 "gram.y" /* yacc.c:1646  */
    { set_axis_prop(whichgraph, naxis, (yyvsp[-1].pset), (yyvsp[0].pset)); }
#line 7882 "y.tab.c" /* yacc.c:1646  */
    break;

  case 498:
#line 1805 "gram.y" /* yacc.c:1646  */
    {}
#line 7888 "y.tab.c" /* yacc.c:1646  */
    break;

  case 499:
#line 1806 "gram.y" /* yacc.c:1646  */
    {}
#line 7894 "y.tab.c" /* yacc.c:1646  */
    break;

  case 500:
#line 1807 "gram.y" /* yacc.c:1646  */
    {}
#line 7900 "y.tab.c" /* yacc.c:1646  */
    break;

  case 501:
#line 1808 "gram.y" /* yacc.c:1646  */
    {}
#line 7906 "y.tab.c" /* yacc.c:1646  */
    break;

  case 502:
#line 1809 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].active = (yyvsp[0].pset); }
#line 7912 "y.tab.c" /* yacc.c:1646  */
    break;

  case 503:
#line 1819 "gram.y" /* yacc.c:1646  */
    {
	    g[curg].t[naxis].t_flag = (yyvsp[0].pset);
	    g[curg].t[naxis].t_mflag = (yyvsp[0].pset);
	}
#line 7921 "y.tab.c" /* yacc.c:1646  */
    break;

  case 504:
#line 1823 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_flag = (yyvsp[0].pset); }
#line 7927 "y.tab.c" /* yacc.c:1646  */
    break;

  case 505:
#line 1824 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_mflag = (yyvsp[0].pset); }
#line 7933 "y.tab.c" /* yacc.c:1646  */
    break;

  case 506:
#line 1825 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tmajor = (yyvsp[0].val); }
#line 7939 "y.tab.c" /* yacc.c:1646  */
    break;

  case 507:
#line 1826 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tminor = (yyvsp[0].val); }
#line 7945 "y.tab.c" /* yacc.c:1646  */
    break;

  case 508:
#line 1827 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].offsx = (yyvsp[0].val); }
#line 7951 "y.tab.c" /* yacc.c:1646  */
    break;

  case 509:
#line 1828 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].offsy = (yyvsp[0].val); }
#line 7957 "y.tab.c" /* yacc.c:1646  */
    break;

  case 510:
#line 1829 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].alt = (yyvsp[0].pset); }
#line 7963 "y.tab.c" /* yacc.c:1646  */
    break;

  case 511:
#line 1830 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tmin = (yyvsp[0].val); }
#line 7969 "y.tab.c" /* yacc.c:1646  */
    break;

  case 512:
#line 1831 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tmax = (yyvsp[0].val); }
#line 7975 "y.tab.c" /* yacc.c:1646  */
    break;

  case 513:
#line 1832 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_num = (int) (yyvsp[0].val); }
#line 7981 "y.tab.c" /* yacc.c:1646  */
    break;

  case 514:
#line 1833 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_inout = (yyvsp[0].pset); }
#line 7987 "y.tab.c" /* yacc.c:1646  */
    break;

  case 515:
#line 1834 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_log = (yyvsp[0].pset); }
#line 7993 "y.tab.c" /* yacc.c:1646  */
    break;

  case 516:
#line 1835 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_size = (yyvsp[0].val); }
#line 7999 "y.tab.c" /* yacc.c:1646  */
    break;

  case 517:
#line 1836 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_size = (yyvsp[0].val); }
#line 8005 "y.tab.c" /* yacc.c:1646  */
    break;

  case 518:
#line 1837 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_msize = (yyvsp[0].val); }
#line 8011 "y.tab.c" /* yacc.c:1646  */
    break;

  case 519:
#line 1838 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_color = g[curg].t[naxis].t_mcolor = (int) (yyvsp[0].val); }
#line 8017 "y.tab.c" /* yacc.c:1646  */
    break;

  case 520:
#line 1839 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_linew = g[curg].t[naxis].t_mlinew = (int) (yyvsp[0].val); }
#line 8023 "y.tab.c" /* yacc.c:1646  */
    break;

  case 521:
#line 1840 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_color = (int) (yyvsp[0].val); }
#line 8029 "y.tab.c" /* yacc.c:1646  */
    break;

  case 522:
#line 1841 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_mcolor = (int) (yyvsp[0].val); }
#line 8035 "y.tab.c" /* yacc.c:1646  */
    break;

  case 523:
#line 1842 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_linew = (int) (yyvsp[0].val); }
#line 8041 "y.tab.c" /* yacc.c:1646  */
    break;

  case 524:
#line 1843 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_mlinew = (int) (yyvsp[0].val); }
#line 8047 "y.tab.c" /* yacc.c:1646  */
    break;

  case 525:
#line 1844 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_lines = (int) (yyvsp[0].val); }
#line 8053 "y.tab.c" /* yacc.c:1646  */
    break;

  case 526:
#line 1845 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_mlines = (int) (yyvsp[0].val); }
#line 8059 "y.tab.c" /* yacc.c:1646  */
    break;

  case 527:
#line 1846 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_gridflag = (yyvsp[0].pset); }
#line 8065 "y.tab.c" /* yacc.c:1646  */
    break;

  case 528:
#line 1847 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_mgridflag = (yyvsp[0].pset); }
#line 8071 "y.tab.c" /* yacc.c:1646  */
    break;

  case 529:
#line 1848 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_op = (yyvsp[0].pset); }
#line 8077 "y.tab.c" /* yacc.c:1646  */
    break;

  case 530:
#line 1849 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_type = AUTO; }
#line 8083 "y.tab.c" /* yacc.c:1646  */
    break;

  case 531:
#line 1850 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_type = SPEC; }
#line 8089 "y.tab.c" /* yacc.c:1646  */
    break;

  case 532:
#line 1851 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_spec = (int) (yyvsp[0].val); }
#line 8095 "y.tab.c" /* yacc.c:1646  */
    break;

  case 533:
#line 1852 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_specloc[(int) (yyvsp[-2].val)] = (yyvsp[0].val); }
#line 8101 "y.tab.c" /* yacc.c:1646  */
    break;

  case 534:
#line 1856 "gram.y" /* yacc.c:1646  */
    {}
#line 8107 "y.tab.c" /* yacc.c:1646  */
    break;

  case 535:
#line 1857 "gram.y" /* yacc.c:1646  */
    {}
#line 8113 "y.tab.c" /* yacc.c:1646  */
    break;

  case 536:
#line 1861 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_flag = (yyvsp[0].pset); }
#line 8119 "y.tab.c" /* yacc.c:1646  */
    break;

  case 537:
#line 1862 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_type = AUTO; }
#line 8125 "y.tab.c" /* yacc.c:1646  */
    break;

  case 538:
#line 1863 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_type = SPEC; }
#line 8131 "y.tab.c" /* yacc.c:1646  */
    break;

  case 539:
#line 1864 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_prec = (int) (yyvsp[0].val); }
#line 8137 "y.tab.c" /* yacc.c:1646  */
    break;

  case 540:
#line 1865 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_format = (yyvsp[0].pset); }
#line 8143 "y.tab.c" /* yacc.c:1646  */
    break;

  case 541:
#line 1866 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_format = (yyvsp[0].val); }
#line 8149 "y.tab.c" /* yacc.c:1646  */
    break;

  case 542:
#line 1867 "gram.y" /* yacc.c:1646  */
    { strcpy(g[curg].t[naxis].tl_appstr, (char *) (yyvsp[0].str)); }
#line 8155 "y.tab.c" /* yacc.c:1646  */
    break;

  case 543:
#line 1868 "gram.y" /* yacc.c:1646  */
    { strcpy(g[curg].t[naxis].tl_prestr, (char *) (yyvsp[0].str)); }
#line 8161 "y.tab.c" /* yacc.c:1646  */
    break;

  case 544:
#line 1869 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_layout = HORIZONTAL; }
#line 8167 "y.tab.c" /* yacc.c:1646  */
    break;

  case 545:
#line 1870 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_layout = VERTICAL; }
#line 8173 "y.tab.c" /* yacc.c:1646  */
    break;

  case 546:
#line 1871 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_layout = SPEC; }
#line 8179 "y.tab.c" /* yacc.c:1646  */
    break;

  case 547:
#line 1872 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_angle = (int) (yyvsp[0].val); }
#line 8185 "y.tab.c" /* yacc.c:1646  */
    break;

  case 548:
#line 1873 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_just = (int) (yyvsp[0].pset); }
#line 8191 "y.tab.c" /* yacc.c:1646  */
    break;

  case 549:
#line 1874 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_skip = (int) (yyvsp[0].val); }
#line 8197 "y.tab.c" /* yacc.c:1646  */
    break;

  case 550:
#line 1875 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_staggered = (int) (yyvsp[0].val); }
#line 8203 "y.tab.c" /* yacc.c:1646  */
    break;

  case 551:
#line 1876 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_op = (yyvsp[0].pset); }
#line 8209 "y.tab.c" /* yacc.c:1646  */
    break;

  case 552:
#line 1877 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_sign = (yyvsp[0].pset); }
#line 8215 "y.tab.c" /* yacc.c:1646  */
    break;

  case 553:
#line 1878 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_start = (yyvsp[0].val); }
#line 8221 "y.tab.c" /* yacc.c:1646  */
    break;

  case 554:
#line 1879 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_stop = (yyvsp[0].val); }
#line 8227 "y.tab.c" /* yacc.c:1646  */
    break;

  case 555:
#line 1880 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_starttype = (int) (yyvsp[0].pset); }
#line 8233 "y.tab.c" /* yacc.c:1646  */
    break;

  case 556:
#line 1881 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_starttype = (int) (yyvsp[0].pset); }
#line 8239 "y.tab.c" /* yacc.c:1646  */
    break;

  case 557:
#line 1882 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_stoptype = (int) (yyvsp[0].pset); }
#line 8245 "y.tab.c" /* yacc.c:1646  */
    break;

  case 558:
#line 1883 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_stoptype = (int) (yyvsp[0].pset); }
#line 8251 "y.tab.c" /* yacc.c:1646  */
    break;

  case 559:
#line 1884 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_vgap = (yyvsp[0].val); }
#line 8257 "y.tab.c" /* yacc.c:1646  */
    break;

  case 560:
#line 1885 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_hgap = (yyvsp[0].val); }
#line 8263 "y.tab.c" /* yacc.c:1646  */
    break;

  case 561:
#line 1886 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_charsize = (yyvsp[0].val); }
#line 8269 "y.tab.c" /* yacc.c:1646  */
    break;

  case 562:
#line 1887 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_font = (int) (yyvsp[0].val); }
#line 8275 "y.tab.c" /* yacc.c:1646  */
    break;

  case 563:
#line 1888 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_color = (int) (yyvsp[0].val); }
#line 8281 "y.tab.c" /* yacc.c:1646  */
    break;

  case 564:
#line 1889 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].tl_linew = (int) (yyvsp[0].val); }
#line 8287 "y.tab.c" /* yacc.c:1646  */
    break;

  case 565:
#line 1890 "gram.y" /* yacc.c:1646  */
    { set_plotstr_string(&g[curg].t[naxis].t_speclab[(int) (yyvsp[-2].val)], (char *) (yyvsp[0].str)); }
#line 8293 "y.tab.c" /* yacc.c:1646  */
    break;

  case 566:
#line 1894 "gram.y" /* yacc.c:1646  */
    { set_plotstr_string(&g[curg].t[naxis].label, (char *) (yyvsp[0].str)); }
#line 8299 "y.tab.c" /* yacc.c:1646  */
    break;

  case 567:
#line 1895 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].label_layout = PERP; }
#line 8305 "y.tab.c" /* yacc.c:1646  */
    break;

  case 568:
#line 1896 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].label_layout = PARA; }
#line 8311 "y.tab.c" /* yacc.c:1646  */
    break;

  case 569:
#line 1897 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].label_place = (yyvsp[0].pset); }
#line 8317 "y.tab.c" /* yacc.c:1646  */
    break;

  case 570:
#line 1898 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].label_place = (yyvsp[0].pset); }
#line 8323 "y.tab.c" /* yacc.c:1646  */
    break;

  case 571:
#line 1899 "gram.y" /* yacc.c:1646  */
    {
	    g[curg].t[naxis].label.x = (yyvsp[-2].val);
	    g[curg].t[naxis].label.y = (yyvsp[0].val);
	}
#line 8332 "y.tab.c" /* yacc.c:1646  */
    break;

  case 572:
#line 1903 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].label.just = (int) (yyvsp[0].pset); }
#line 8338 "y.tab.c" /* yacc.c:1646  */
    break;

  case 573:
#line 1904 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].label.charsize = (yyvsp[0].val); }
#line 8344 "y.tab.c" /* yacc.c:1646  */
    break;

  case 574:
#line 1905 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].label.font = (int) (yyvsp[0].val); }
#line 8350 "y.tab.c" /* yacc.c:1646  */
    break;

  case 575:
#line 1906 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].label.color = (int) (yyvsp[0].val); }
#line 8356 "y.tab.c" /* yacc.c:1646  */
    break;

  case 576:
#line 1907 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].label.linew = (int) (yyvsp[0].val); }
#line 8362 "y.tab.c" /* yacc.c:1646  */
    break;

  case 577:
#line 1911 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_drawbar = (yyvsp[0].pset); }
#line 8368 "y.tab.c" /* yacc.c:1646  */
    break;

  case 578:
#line 1912 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_drawbarcolor = (int) (yyvsp[0].val); }
#line 8374 "y.tab.c" /* yacc.c:1646  */
    break;

  case 579:
#line 1913 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_drawbarlines = (int) (yyvsp[0].val); }
#line 8380 "y.tab.c" /* yacc.c:1646  */
    break;

  case 580:
#line 1914 "gram.y" /* yacc.c:1646  */
    { g[curg].t[naxis].t_drawbarlinew = (int) (yyvsp[0].val); }
#line 8386 "y.tab.c" /* yacc.c:1646  */
    break;

  case 581:
#line 1926 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = DISK; }
#line 8392 "y.tab.c" /* yacc.c:1646  */
    break;

  case 582:
#line 1927 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = PIPE; }
#line 8398 "y.tab.c" /* yacc.c:1646  */
    break;

  case 583:
#line 1931 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = RIGHT; }
#line 8404 "y.tab.c" /* yacc.c:1646  */
    break;

  case 584:
#line 1932 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = LEFT; }
#line 8410 "y.tab.c" /* yacc.c:1646  */
    break;

  case 585:
#line 1933 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = CENTER; }
#line 8416 "y.tab.c" /* yacc.c:1646  */
    break;

  case 586:
#line 1942 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = (yyvsp[0].pset); }
#line 8422 "y.tab.c" /* yacc.c:1646  */
    break;

  case 587:
#line 1943 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = (yyvsp[0].pset); }
#line 8428 "y.tab.c" /* yacc.c:1646  */
    break;

  case 588:
#line 1944 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = (yyvsp[0].pset); }
#line 8434 "y.tab.c" /* yacc.c:1646  */
    break;

  case 589:
#line 1945 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = (yyvsp[0].pset); }
#line 8440 "y.tab.c" /* yacc.c:1646  */
    break;

  case 590:
#line 1946 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = XYFIXED; }
#line 8446 "y.tab.c" /* yacc.c:1646  */
    break;

  case 591:
#line 1950 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = IN; }
#line 8452 "y.tab.c" /* yacc.c:1646  */
    break;

  case 592:
#line 1951 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = OUT; }
#line 8458 "y.tab.c" /* yacc.c:1646  */
    break;

  case 593:
#line 1952 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = BOTH; }
#line 8464 "y.tab.c" /* yacc.c:1646  */
    break;

  case 594:
#line 1956 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = NORMAL; }
#line 8470 "y.tab.c" /* yacc.c:1646  */
    break;

  case 595:
#line 1957 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = ABSOLUTE; }
#line 8476 "y.tab.c" /* yacc.c:1646  */
    break;

  case 596:
#line 1958 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = NEGATE; }
#line 8482 "y.tab.c" /* yacc.c:1646  */
    break;

  case 597:
#line 1971 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = DECIMAL; }
#line 8488 "y.tab.c" /* yacc.c:1646  */
    break;

  case 598:
#line 1972 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = EXPONENTIAL; }
#line 8494 "y.tab.c" /* yacc.c:1646  */
    break;

  case 599:
#line 1973 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = POWER; }
#line 8500 "y.tab.c" /* yacc.c:1646  */
    break;

  case 600:
#line 1974 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = GENERAL; }
#line 8506 "y.tab.c" /* yacc.c:1646  */
    break;

  case 601:
#line 1975 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = DDMMYY; }
#line 8512 "y.tab.c" /* yacc.c:1646  */
    break;

  case 602:
#line 1976 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = MMDDYY; }
#line 8518 "y.tab.c" /* yacc.c:1646  */
    break;

  case 603:
#line 1977 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = MMYY; }
#line 8524 "y.tab.c" /* yacc.c:1646  */
    break;

  case 604:
#line 1978 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = MMDD; }
#line 8530 "y.tab.c" /* yacc.c:1646  */
    break;

  case 605:
#line 1979 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = MONTHDAY; }
#line 8536 "y.tab.c" /* yacc.c:1646  */
    break;

  case 606:
#line 1980 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = DAYMONTH; }
#line 8542 "y.tab.c" /* yacc.c:1646  */
    break;

  case 607:
#line 1981 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = DDMONTHSYYHHMMSS; }
#line 8548 "y.tab.c" /* yacc.c:1646  */
    break;

  case 608:
#line 1982 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = MONTHS; }
#line 8554 "y.tab.c" /* yacc.c:1646  */
    break;

  case 609:
#line 1983 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = MONTHL; }
#line 8560 "y.tab.c" /* yacc.c:1646  */
    break;

  case 610:
#line 1984 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = DAYOFWEEKS; }
#line 8566 "y.tab.c" /* yacc.c:1646  */
    break;

  case 611:
#line 1985 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = DAYOFWEEKL; }
#line 8572 "y.tab.c" /* yacc.c:1646  */
    break;

  case 612:
#line 1986 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = DAYOFYEAR; }
#line 8578 "y.tab.c" /* yacc.c:1646  */
    break;

  case 613:
#line 1987 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = HMS; }
#line 8584 "y.tab.c" /* yacc.c:1646  */
    break;

  case 614:
#line 1988 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = MMDDHMS; }
#line 8590 "y.tab.c" /* yacc.c:1646  */
    break;

  case 615:
#line 1989 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = MMDDYYHMS; }
#line 8596 "y.tab.c" /* yacc.c:1646  */
    break;

  case 616:
#line 1990 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = DEGREESLON; }
#line 8602 "y.tab.c" /* yacc.c:1646  */
    break;

  case 617:
#line 1991 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = DEGREESMMLON; }
#line 8608 "y.tab.c" /* yacc.c:1646  */
    break;

  case 618:
#line 1992 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = DEGREESMMSSLON; }
#line 8614 "y.tab.c" /* yacc.c:1646  */
    break;

  case 619:
#line 1993 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = MMSSLON; }
#line 8620 "y.tab.c" /* yacc.c:1646  */
    break;

  case 620:
#line 1994 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = DEGREESLAT; }
#line 8626 "y.tab.c" /* yacc.c:1646  */
    break;

  case 621:
#line 1995 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = DEGREESMMLAT; }
#line 8632 "y.tab.c" /* yacc.c:1646  */
    break;

  case 622:
#line 1996 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = DEGREESMMSSLAT; }
#line 8638 "y.tab.c" /* yacc.c:1646  */
    break;

  case 623:
#line 1997 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = MMSSLAT; }
#line 8644 "y.tab.c" /* yacc.c:1646  */
    break;

  case 624:
#line 2000 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = (yyvsp[0].pset); }
#line 8650 "y.tab.c" /* yacc.c:1646  */
    break;

  case 625:
#line 2001 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = (yyvsp[0].pset); }
#line 8656 "y.tab.c" /* yacc.c:1646  */
    break;

  case 626:
#line 2005 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = (yyvsp[0].pset); }
#line 8662 "y.tab.c" /* yacc.c:1646  */
    break;

  case 627:
#line 2006 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = (yyvsp[0].pset); }
#line 8668 "y.tab.c" /* yacc.c:1646  */
    break;

  case 628:
#line 2007 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = (yyvsp[0].pset); }
#line 8674 "y.tab.c" /* yacc.c:1646  */
    break;

  case 629:
#line 2008 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = (yyvsp[0].pset); }
#line 8680 "y.tab.c" /* yacc.c:1646  */
    break;

  case 630:
#line 2012 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = ON; }
#line 8686 "y.tab.c" /* yacc.c:1646  */
    break;

  case 631:
#line 2013 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = OFF; }
#line 8692 "y.tab.c" /* yacc.c:1646  */
    break;

  case 632:
#line 2016 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = WORLD; }
#line 8698 "y.tab.c" /* yacc.c:1646  */
    break;

  case 633:
#line 2017 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = VIEW; }
#line 8704 "y.tab.c" /* yacc.c:1646  */
    break;

  case 634:
#line 2028 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = TRUEP; }
#line 8710 "y.tab.c" /* yacc.c:1646  */
    break;

  case 635:
#line 2029 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = FALSEP; }
#line 8716 "y.tab.c" /* yacc.c:1646  */
    break;

  case 636:
#line 2032 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = PATTERN; }
#line 8722 "y.tab.c" /* yacc.c:1646  */
    break;

  case 637:
#line 2033 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = COLOR; }
#line 8728 "y.tab.c" /* yacc.c:1646  */
    break;

  case 638:
#line 2034 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = NONE; }
#line 8734 "y.tab.c" /* yacc.c:1646  */
    break;

  case 639:
#line 2037 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = ABOVE; }
#line 8740 "y.tab.c" /* yacc.c:1646  */
    break;

  case 640:
#line 2038 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = BELOW; }
#line 8746 "y.tab.c" /* yacc.c:1646  */
    break;

  case 641:
#line 2039 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = LEFT; }
#line 8752 "y.tab.c" /* yacc.c:1646  */
    break;

  case 642:
#line 2040 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = RIGHT; }
#line 8758 "y.tab.c" /* yacc.c:1646  */
    break;

  case 643:
#line 2041 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = TOP; }
#line 8764 "y.tab.c" /* yacc.c:1646  */
    break;

  case 644:
#line 2042 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = BOTTOM ; }
#line 8770 "y.tab.c" /* yacc.c:1646  */
    break;

  case 645:
#line 2043 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = BOTH ; }
#line 8776 "y.tab.c" /* yacc.c:1646  */
    break;

  case 646:
#line 2047 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = MM; }
#line 8782 "y.tab.c" /* yacc.c:1646  */
    break;

  case 647:
#line 2048 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = CM; }
#line 8788 "y.tab.c" /* yacc.c:1646  */
    break;

  case 648:
#line 2049 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = M; }
#line 8794 "y.tab.c" /* yacc.c:1646  */
    break;

  case 649:
#line 2050 "gram.y" /* yacc.c:1646  */
    { (yyval.pset) = KM; }
#line 8800 "y.tab.c" /* yacc.c:1646  */
    break;

  case 650:
#line 2065 "gram.y" /* yacc.c:1646  */
    {
	    int itmp = (int) (yyvsp[-3].val) - 1;
	    if (itmp >= ls) {
		yyerror("subscript out of range");
		return 1;
	    } else {
		(yyvsp[-5].ptr)[itmp] = (yyvsp[0].val);
		result = (yyvsp[0].val);
	    }
	}
#line 8815 "y.tab.c" /* yacc.c:1646  */
    break;

  case 651:
#line 2079 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    for (i = 0; i < lxy; i++) {
		(yyvsp[-2].ptr)[i] = (yyvsp[0].ptr)[i];
	    }
	    result = (yyvsp[0].ptr)[0];
	}
#line 8827 "y.tab.c" /* yacc.c:1646  */
    break;

  case 652:
#line 2087 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    for (i = 0; i < lxy; i++) {
		(yyvsp[-2].ptr)[i] = (yyvsp[0].val);
	    }
	    result = (yyvsp[0].val);
	}
#line 8839 "y.tab.c" /* yacc.c:1646  */
    break;

  case 653:
#line 2098 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[0].ptr)[i];
	    }
	}
#line 8852 "y.tab.c" /* yacc.c:1646  */
    break;

  case 654:
#line 2107 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[0].val);
	    }
	}
#line 8865 "y.tab.c" /* yacc.c:1646  */
    break;

  case 655:
#line 2116 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].val) + (yyvsp[0].val);
	    }
	}
#line 8878 "y.tab.c" /* yacc.c:1646  */
    break;

  case 656:
#line 2125 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] + (yyvsp[0].ptr)[i];
	    }
	}
#line 8891 "y.tab.c" /* yacc.c:1646  */
    break;

  case 657:
#line 2134 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].val) + (yyvsp[0].ptr)[i];
	    }
	}
#line 8904 "y.tab.c" /* yacc.c:1646  */
    break;

  case 658:
#line 2143 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] + (yyvsp[0].val);
	    }
	}
#line 8917 "y.tab.c" /* yacc.c:1646  */
    break;

  case 659:
#line 2152 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].val) - (yyvsp[0].val);
	    }
	}
#line 8930 "y.tab.c" /* yacc.c:1646  */
    break;

  case 660:
#line 2161 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] - (yyvsp[0].ptr)[i];
	    }
	}
#line 8943 "y.tab.c" /* yacc.c:1646  */
    break;

  case 661:
#line 2170 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].val) - (yyvsp[0].ptr)[i];
	    }
	}
#line 8956 "y.tab.c" /* yacc.c:1646  */
    break;

  case 662:
#line 2179 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] - (yyvsp[0].val);
	    }
	}
#line 8969 "y.tab.c" /* yacc.c:1646  */
    break;

  case 663:
#line 2188 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].val) * (yyvsp[0].val);
	    }
	}
#line 8982 "y.tab.c" /* yacc.c:1646  */
    break;

  case 664:
#line 2197 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] * (yyvsp[0].ptr)[i];
	    }
	}
#line 8995 "y.tab.c" /* yacc.c:1646  */
    break;

  case 665:
#line 2206 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].val) * (yyvsp[0].ptr)[i];
	    }
	}
#line 9008 "y.tab.c" /* yacc.c:1646  */
    break;

  case 666:
#line 2215 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] * (yyvsp[0].val);
	    }
	}
#line 9021 "y.tab.c" /* yacc.c:1646  */
    break;

  case 667:
#line 2224 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    if ((yyvsp[0].val) == 0.0) {
		yyerror("Divide by Zero");
		return 1;
	    }
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].val) / (yyvsp[0].val);
	    }
	}
#line 9038 "y.tab.c" /* yacc.c:1646  */
    break;

  case 668:
#line 2237 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		if ((yyvsp[0].ptr)[i] == 0.0) {
		    yyerror("Divide by Zero");
		    return 1;
		}
	    }
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] / (yyvsp[0].ptr)[i];
	    }
	}
#line 9057 "y.tab.c" /* yacc.c:1646  */
    break;

  case 669:
#line 2252 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		if ((yyvsp[0].ptr)[i] == 0.0) {
		    yyerror("Divide by Zero");
		    return 1;
		}
	    }
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].val) / (yyvsp[0].ptr)[i];
	    }
	}
#line 9076 "y.tab.c" /* yacc.c:1646  */
    break;

  case 670:
#line 2267 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    if ((yyvsp[0].val) == 0.0) {
		yyerror("Divide by Zero");
		return 1;
	    }
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] / (yyvsp[0].val);
	    }
	}
#line 9093 "y.tab.c" /* yacc.c:1646  */
    break;

  case 671:
#line 2280 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = pow((yyvsp[-2].val), (yyvsp[0].val));
	    }
	}
#line 9106 "y.tab.c" /* yacc.c:1646  */
    break;

  case 672:
#line 2289 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = pow((yyvsp[-2].val), (yyvsp[0].ptr)[i]);
	    }
	}
#line 9119 "y.tab.c" /* yacc.c:1646  */
    break;

  case 673:
#line 2298 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = pow((yyvsp[-2].ptr)[i], (yyvsp[0].val));
	    }
	}
#line 9132 "y.tab.c" /* yacc.c:1646  */
    break;

  case 674:
#line 2307 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = pow((yyvsp[-2].ptr)[i], (yyvsp[0].ptr)[i]);
	    }
	}
#line 9145 "y.tab.c" /* yacc.c:1646  */
    break;

  case 675:
#line 2316 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = fabs((yyvsp[-1].val));
	    }
	}
#line 9158 "y.tab.c" /* yacc.c:1646  */
    break;

  case 676:
#line 2325 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = fabs((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9171 "y.tab.c" /* yacc.c:1646  */
    break;

  case 677:
#line 2334 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = acos((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9184 "y.tab.c" /* yacc.c:1646  */
    break;

  case 678:
#line 2343 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = asin((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9197 "y.tab.c" /* yacc.c:1646  */
    break;

  case 679:
#line 2352 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = atan((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9210 "y.tab.c" /* yacc.c:1646  */
    break;

  case 680:
#line 2361 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = atan2((yyvsp[-3].ptr)[i], (yyvsp[-1].ptr)[i]);
	    }
	}
#line 9223 "y.tab.c" /* yacc.c:1646  */
    break;

  case 681:
#line 2370 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = ceil((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9236 "y.tab.c" /* yacc.c:1646  */
    break;

  case 682:
#line 2379 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = cos((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9249 "y.tab.c" /* yacc.c:1646  */
    break;

  case 683:
#line 2388 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] *= M_PI / 180.0;
	    }
	}
#line 9262 "y.tab.c" /* yacc.c:1646  */
    break;

  case 684:
#line 2397 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = xx[i];
	    }
	}
#line 9275 "y.tab.c" /* yacc.c:1646  */
    break;

  case 685:
#line 2406 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = yy[i];
	    }
	}
#line 9288 "y.tab.c" /* yacc.c:1646  */
    break;

  case 686:
#line 2415 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = erf((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9301 "y.tab.c" /* yacc.c:1646  */
    break;

  case 687:
#line 2424 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = erfc((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9314 "y.tab.c" /* yacc.c:1646  */
    break;

  case 688:
#line 2433 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = exp((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9327 "y.tab.c" /* yacc.c:1646  */
    break;

  case 689:
#line 2442 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = floor((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9340 "y.tab.c" /* yacc.c:1646  */
    break;

  case 690:
#line 2451 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = hypot((yyvsp[-3].ptr)[i], (yyvsp[-1].ptr)[i]);
	    }
	}
#line 9353 "y.tab.c" /* yacc.c:1646  */
    break;

  case 691:
#line 2460 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = hypot((yyvsp[-3].val), (yyvsp[-1].ptr)[i]);
	    }
	}
#line 9366 "y.tab.c" /* yacc.c:1646  */
    break;

  case 692:
#line 2469 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = hypot((yyvsp[-3].ptr)[i], (yyvsp[-1].val));
	    }
	}
#line 9379 "y.tab.c" /* yacc.c:1646  */
    break;

  case 693:
#line 2478 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = hypot((yyvsp[-3].val), (yyvsp[-1].val));
	    }
	}
#line 9392 "y.tab.c" /* yacc.c:1646  */
    break;

  case 694:
#line 2487 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = i + 1;
	    }
	}
#line 9405 "y.tab.c" /* yacc.c:1646  */
    break;

  case 695:
#line 2496 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[0].func);
	    }
	}
#line 9418 "y.tab.c" /* yacc.c:1646  */
    break;

  case 696:
#line 2505 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (int) (yyvsp[-1].ptr)[i];
	    }
	}
#line 9431 "y.tab.c" /* yacc.c:1646  */
    break;

  case 697:
#line 2514 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = lrand48() % (long) ((yyvsp[-1].val));
	    }
	}
#line 9444 "y.tab.c" /* yacc.c:1646  */
    break;

  case 698:
#line 2523 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = lgamma((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9457 "y.tab.c" /* yacc.c:1646  */
    break;

  case 699:
#line 2532 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = log((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9470 "y.tab.c" /* yacc.c:1646  */
    break;

  case 700:
#line 2541 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = log10((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9483 "y.tab.c" /* yacc.c:1646  */
    break;

  case 701:
#line 2550 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = 1.0 / (1.0 + exp(-((yyvsp[-5].ptr)[i] - (yyvsp[-3].val))/ (yyvsp[-1].val)));
	    }
	}
#line 9496 "y.tab.c" /* yacc.c:1646  */
    break;

  case 702:
#line 2559 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-3].ptr)[i] >= (yyvsp[-1].ptr)[i] ? (yyvsp[-3].ptr)[i] : (yyvsp[-1].ptr)[i];
	    }
	}
#line 9509 "y.tab.c" /* yacc.c:1646  */
    break;

  case 703:
#line 2568 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-3].ptr)[i] <= (yyvsp[-1].ptr)[i] ? (yyvsp[-3].ptr)[i] : (yyvsp[-1].ptr)[i];
	    }
	}
#line 9522 "y.tab.c" /* yacc.c:1646  */
    break;

  case 704:
#line 2577 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = fmod((yyvsp[-3].ptr)[i], (yyvsp[-1].ptr)[i]);
	    }
	}
#line 9535 "y.tab.c" /* yacc.c:1646  */
    break;

  case 705:
#line 2586 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = fx((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9548 "y.tab.c" /* yacc.c:1646  */
    break;

  case 706:
#line 2595 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    double tmp;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = normp((yyvsp[-1].ptr)[i], &tmp);
	    }
	}
#line 9562 "y.tab.c" /* yacc.c:1646  */
    break;

  case 707:
#line 2605 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = M_PI;
	    }
	}
#line 9575 "y.tab.c" /* yacc.c:1646  */
    break;

  case 708:
#line 2614 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = M_PI / 180.0;
	    }
	}
#line 9588 "y.tab.c" /* yacc.c:1646  */
    break;

  case 709:
#line 2623 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (double) drand48();
	    }
	}
#line 9601 "y.tab.c" /* yacc.c:1646  */
    break;

  case 710:
#line 2632 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = sin((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9614 "y.tab.c" /* yacc.c:1646  */
    break;

  case 711:
#line 2641 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-1].ptr)[i] * (yyvsp[-1].ptr)[i];
	    }
	}
#line 9627 "y.tab.c" /* yacc.c:1646  */
    break;

  case 712:
#line 2650 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = sqrt((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9640 "y.tab.c" /* yacc.c:1646  */
    break;

  case 713:
#line 2659 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = tan((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9653 "y.tab.c" /* yacc.c:1646  */
    break;

  case 714:
#line 2668 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] > (yyvsp[0].ptr)[i];
	    }
	}
#line 9666 "y.tab.c" /* yacc.c:1646  */
    break;

  case 715:
#line 2677 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] < (yyvsp[0].ptr)[i];
	    }
	}
#line 9679 "y.tab.c" /* yacc.c:1646  */
    break;

  case 716:
#line 2686 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] <= (yyvsp[0].ptr)[i];
	    }
	}
#line 9692 "y.tab.c" /* yacc.c:1646  */
    break;

  case 717:
#line 2695 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] >= (yyvsp[0].ptr)[i];
	    }
	}
#line 9705 "y.tab.c" /* yacc.c:1646  */
    break;

  case 718:
#line 2704 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] == (yyvsp[0].ptr)[i];
	    }
	}
#line 9718 "y.tab.c" /* yacc.c:1646  */
    break;

  case 719:
#line 2713 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] != (yyvsp[0].ptr)[i];
	    }
	}
#line 9731 "y.tab.c" /* yacc.c:1646  */
    break;

  case 720:
#line 2722 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] && (yyvsp[0].ptr)[i];
	    }
	}
#line 9744 "y.tab.c" /* yacc.c:1646  */
    break;

  case 721:
#line 2731 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] || (yyvsp[0].ptr)[i];
	    }
	}
#line 9757 "y.tab.c" /* yacc.c:1646  */
    break;

  case 722:
#line 2740 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = !((yyvsp[0].ptr)[i]);
	    }
	}
#line 9770 "y.tab.c" /* yacc.c:1646  */
    break;

  case 723:
#line 2749 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-1].ptr)[i];
	    }
	}
#line 9783 "y.tab.c" /* yacc.c:1646  */
    break;

  case 724:
#line 2757 "gram.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = -(yyvsp[0].ptr)[i];
	    }
	}
#line 9796 "y.tab.c" /* yacc.c:1646  */
    break;

  case 726:
#line 2768 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = (yyvsp[0].val);
	}
#line 9804 "y.tab.c" /* yacc.c:1646  */
    break;

  case 727:
#line 2771 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = (yyvsp[-3].ptr)[(int) (yyvsp[-1].val)];
	}
#line 9812 "y.tab.c" /* yacc.c:1646  */
    break;

  case 728:
#line 2774 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = (yyvsp[-2].val) + (yyvsp[0].val);
	}
#line 9820 "y.tab.c" /* yacc.c:1646  */
    break;

  case 729:
#line 2777 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = (yyvsp[-2].val) - (yyvsp[0].val);
	}
#line 9828 "y.tab.c" /* yacc.c:1646  */
    break;

  case 730:
#line 2780 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = (yyvsp[-2].val) * (yyvsp[0].val);
	}
#line 9836 "y.tab.c" /* yacc.c:1646  */
    break;

  case 731:
#line 2784 "gram.y" /* yacc.c:1646  */
    {
	    if ((yyvsp[0].val) != 0.0) {
		(yyval.val) = (yyvsp[-2].val) / (yyvsp[0].val);
	    } else {
		yyerror("Divide by Zero");
		return 1;
	    }
	}
#line 9849 "y.tab.c" /* yacc.c:1646  */
    break;

  case 732:
#line 2792 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = fmod((yyvsp[-2].val), (yyvsp[0].val));
	}
#line 9857 "y.tab.c" /* yacc.c:1646  */
    break;

  case 733:
#line 2795 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = pow((yyvsp[-2].val), (yyvsp[0].val));
	}
#line 9865 "y.tab.c" /* yacc.c:1646  */
    break;

  case 734:
#line 2798 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = fabs((yyvsp[-1].val));
	}
#line 9873 "y.tab.c" /* yacc.c:1646  */
    break;

  case 735:
#line 2801 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = acos((yyvsp[-1].val));
	}
#line 9881 "y.tab.c" /* yacc.c:1646  */
    break;

  case 736:
#line 2804 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = asin((yyvsp[-1].val));
	}
#line 9889 "y.tab.c" /* yacc.c:1646  */
    break;

  case 737:
#line 2807 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = atan((yyvsp[-1].val));
	}
#line 9897 "y.tab.c" /* yacc.c:1646  */
    break;

  case 738:
#line 2810 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = atan2((yyvsp[-3].val), (yyvsp[-1].val));
	}
#line 9905 "y.tab.c" /* yacc.c:1646  */
    break;

  case 739:
#line 2813 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = ceil((yyvsp[-1].val));
	}
#line 9913 "y.tab.c" /* yacc.c:1646  */
    break;

  case 740:
#line 2816 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = cos((yyvsp[-1].val));
	}
#line 9921 "y.tab.c" /* yacc.c:1646  */
    break;

  case 741:
#line 2819 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = 180.0 / M_PI;
	}
#line 9929 "y.tab.c" /* yacc.c:1646  */
    break;

  case 742:
#line 2822 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = *xx;
	}
#line 9937 "y.tab.c" /* yacc.c:1646  */
    break;

  case 743:
#line 2825 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = *yy;
	}
#line 9945 "y.tab.c" /* yacc.c:1646  */
    break;

  case 744:
#line 2828 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = erf((yyvsp[-1].val));
	}
#line 9953 "y.tab.c" /* yacc.c:1646  */
    break;

  case 745:
#line 2831 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = erfc((yyvsp[-1].val));
	}
#line 9961 "y.tab.c" /* yacc.c:1646  */
    break;

  case 746:
#line 2834 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = exp((yyvsp[-1].val));
	}
#line 9969 "y.tab.c" /* yacc.c:1646  */
    break;

  case 747:
#line 2837 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = floor((yyvsp[-1].val));
	}
#line 9977 "y.tab.c" /* yacc.c:1646  */
    break;

  case 748:
#line 2840 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = hypot((yyvsp[-3].val), (yyvsp[-1].val));
	}
#line 9985 "y.tab.c" /* yacc.c:1646  */
    break;

  case 749:
#line 2843 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = g[(yyvsp[-2].pset)].v.xv1;
	}
#line 9993 "y.tab.c" /* yacc.c:1646  */
    break;

  case 750:
#line 2846 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = g[(yyvsp[-2].pset)].v.xv2;
	}
#line 10001 "y.tab.c" /* yacc.c:1646  */
    break;

  case 751:
#line 2849 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = g[(yyvsp[-2].pset)].v.yv1;
	}
#line 10009 "y.tab.c" /* yacc.c:1646  */
    break;

  case 752:
#line 2852 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = g[(yyvsp[-2].pset)].v.yv2;
	}
#line 10017 "y.tab.c" /* yacc.c:1646  */
    break;

  case 753:
#line 2855 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = g[(yyvsp[-2].pset)].w.xg1;
	}
#line 10025 "y.tab.c" /* yacc.c:1646  */
    break;

  case 754:
#line 2858 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = g[(yyvsp[-2].pset)].w.xg2;
	}
#line 10033 "y.tab.c" /* yacc.c:1646  */
    break;

  case 755:
#line 2861 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = g[(yyvsp[-2].pset)].w.yg1;
	}
#line 10041 "y.tab.c" /* yacc.c:1646  */
    break;

  case 756:
#line 2864 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = g[(yyvsp[-2].pset)].w.yg2;
	}
#line 10049 "y.tab.c" /* yacc.c:1646  */
    break;

  case 757:
#line 2867 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = g[curg].v.xv1;
	}
#line 10057 "y.tab.c" /* yacc.c:1646  */
    break;

  case 758:
#line 2870 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = g[curg].v.xv2;
	}
#line 10065 "y.tab.c" /* yacc.c:1646  */
    break;

  case 759:
#line 2873 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = g[curg].v.yv1;
	}
#line 10073 "y.tab.c" /* yacc.c:1646  */
    break;

  case 760:
#line 2876 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = g[curg].v.yv2;
	}
#line 10081 "y.tab.c" /* yacc.c:1646  */
    break;

  case 761:
#line 2879 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = g[curg].w.xg1;
	}
#line 10089 "y.tab.c" /* yacc.c:1646  */
    break;

  case 762:
#line 2882 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = g[curg].w.xg2;
	}
#line 10097 "y.tab.c" /* yacc.c:1646  */
    break;

  case 763:
#line 2885 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = g[curg].w.yg1;
	}
#line 10105 "y.tab.c" /* yacc.c:1646  */
    break;

  case 764:
#line 2888 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = g[curg].w.yg2;
	}
#line 10113 "y.tab.c" /* yacc.c:1646  */
    break;

  case 765:
#line 2891 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = setindex;
	}
#line 10121 "y.tab.c" /* yacc.c:1646  */
    break;

  case 766:
#line 2894 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = setsetno;
	}
#line 10129 "y.tab.c" /* yacc.c:1646  */
    break;

  case 767:
#line 2897 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = (long) (yyvsp[-1].val);
	}
#line 10137 "y.tab.c" /* yacc.c:1646  */
    break;

  case 768:
#line 2900 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = lrand48() % (long) ((yyvsp[-1].val));
	}
#line 10145 "y.tab.c" /* yacc.c:1646  */
    break;

  case 769:
#line 2903 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = lgamma((yyvsp[-1].val));
	}
#line 10153 "y.tab.c" /* yacc.c:1646  */
    break;

  case 770:
#line 2906 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = log((yyvsp[-1].val));
	}
#line 10161 "y.tab.c" /* yacc.c:1646  */
    break;

  case 771:
#line 2909 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = log10((yyvsp[-1].val));
	}
#line 10169 "y.tab.c" /* yacc.c:1646  */
    break;

  case 772:
#line 2913 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = 1.0 / (1.0 + exp(-((yyvsp[-5].val) - (yyvsp[-3].val))/ (yyvsp[-1].val)));
	}
#line 10177 "y.tab.c" /* yacc.c:1646  */
    break;

  case 773:
#line 2916 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = (yyvsp[-3].val) >= (yyvsp[-1].val) ? (yyvsp[-3].val) : (yyvsp[-1].val);
	}
#line 10185 "y.tab.c" /* yacc.c:1646  */
    break;

  case 774:
#line 2919 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = (yyvsp[-3].val) <= (yyvsp[-1].val) ? (yyvsp[-3].val) : (yyvsp[-1].val);
	}
#line 10193 "y.tab.c" /* yacc.c:1646  */
    break;

  case 775:
#line 2922 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = fmod((yyvsp[-3].val), (yyvsp[-1].val));
	}
#line 10201 "y.tab.c" /* yacc.c:1646  */
    break;

  case 776:
#line 2925 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = fx((yyvsp[-1].val));
	}
#line 10209 "y.tab.c" /* yacc.c:1646  */
    break;

  case 777:
#line 2928 "gram.y" /* yacc.c:1646  */
    {
	    double tmp;
	    (yyval.val) = normp((yyvsp[-1].val), &tmp);
	}
#line 10218 "y.tab.c" /* yacc.c:1646  */
    break;

  case 778:
#line 2932 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = M_PI;
	}
#line 10226 "y.tab.c" /* yacc.c:1646  */
    break;

  case 779:
#line 2935 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = M_PI / 180.0;
	}
#line 10234 "y.tab.c" /* yacc.c:1646  */
    break;

  case 780:
#line 2938 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = (double) drand48();
	}
#line 10242 "y.tab.c" /* yacc.c:1646  */
    break;

  case 781:
#line 2941 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = sin((yyvsp[-1].val));
	}
#line 10250 "y.tab.c" /* yacc.c:1646  */
    break;

  case 782:
#line 2944 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = pow((yyvsp[-1].val), 2.0);
	}
#line 10258 "y.tab.c" /* yacc.c:1646  */
    break;

  case 783:
#line 2947 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = sqrt((yyvsp[-1].val));
	}
#line 10266 "y.tab.c" /* yacc.c:1646  */
    break;

  case 784:
#line 2950 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = tan((yyvsp[-1].val));
	}
#line 10274 "y.tab.c" /* yacc.c:1646  */
    break;

  case 785:
#line 2953 "gram.y" /* yacc.c:1646  */
    {
	    if ((int) (yyvsp[-2].val))
		(yyval.val) = (yyvsp[0].val);
	}
#line 10283 "y.tab.c" /* yacc.c:1646  */
    break;

  case 786:
#line 2957 "gram.y" /* yacc.c:1646  */
    {
	    if ((int) (yyvsp[-4].val)) {
		(yyval.val) = (yyvsp[-2].val);
	    } else {
		(yyval.val) = (yyvsp[0].val);
	    }
	}
#line 10295 "y.tab.c" /* yacc.c:1646  */
    break;

  case 787:
#line 2964 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = (yyvsp[-2].val) > (yyvsp[0].val);
	}
#line 10303 "y.tab.c" /* yacc.c:1646  */
    break;

  case 788:
#line 2967 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = (yyvsp[-2].val) < (yyvsp[0].val);
	}
#line 10311 "y.tab.c" /* yacc.c:1646  */
    break;

  case 789:
#line 2970 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = (yyvsp[-2].val) <= (yyvsp[0].val);
	}
#line 10319 "y.tab.c" /* yacc.c:1646  */
    break;

  case 790:
#line 2973 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = (yyvsp[-2].val) >= (yyvsp[0].val);
	}
#line 10327 "y.tab.c" /* yacc.c:1646  */
    break;

  case 791:
#line 2976 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = (yyvsp[-2].val) == (yyvsp[0].val);
	}
#line 10335 "y.tab.c" /* yacc.c:1646  */
    break;

  case 792:
#line 2979 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = (yyvsp[-2].val) != (yyvsp[0].val);
	}
#line 10343 "y.tab.c" /* yacc.c:1646  */
    break;

  case 793:
#line 2982 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = (yyvsp[-2].val) && (yyvsp[0].val);
	}
#line 10351 "y.tab.c" /* yacc.c:1646  */
    break;

  case 794:
#line 2985 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = (yyvsp[-2].val) || (yyvsp[0].val);
	}
#line 10359 "y.tab.c" /* yacc.c:1646  */
    break;

  case 795:
#line 2988 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = !((yyvsp[0].val));
	}
#line 10367 "y.tab.c" /* yacc.c:1646  */
    break;

  case 796:
#line 2991 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = (yyvsp[-1].val);
	}
#line 10375 "y.tab.c" /* yacc.c:1646  */
    break;

  case 797:
#line 2994 "gram.y" /* yacc.c:1646  */
    {
	    (yyval.val) = -(yyvsp[0].val);
	}
#line 10383 "y.tab.c" /* yacc.c:1646  */
    break;


#line 10387 "y.tab.c" /* yacc.c:1646  */
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYTERROR;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  return yyresult;
}
#line 2998 "gram.y" /* yacc.c:1906  */


void fixupstr(char *val)
{
    int vl = strlen(val);
    lowtoupper(val);
    val[vl + 1] = 0;
    val[vl] = '\n';
}

void scanner(char *s, double *x, double *y, int len, double *a, double *b, double *c, double *d, int lenscr, int i, int setno, int *errpos)
{
    interr = 0;
    whichgraph = cg;
    whichset = setno;
    if (s[0] == '#') {
	return;
    }
    pos = 0;
    aa = a;
    bb = b;
    cc = c;
    dd = d;
    xx = x;
    yy = y;
    lxy = len;
    ls = lenscr;
    setindex = i + 1;
    curset = setsetno = setno;
    strcpy(f_string, s);
    strcpy(statusstr, s);
    fcnt = 0;
    yyparse();
    *errpos = interr;
    for (i = 0; i < fcnt; i++) {
	free(freelist[i]);
	freelist[i] = NULL;
    }
}

void runbatch(char *bfile)
{
    double x, y, a, b, c, d;
    int i, setno, errpos, lcnt = 1;
    char stext[256];
    FILE *fp;
    if (strcmp("stdin", bfile)) {
	fp = fopen(bfile, "r");
    }
    else {
	fp = stdin;
    }
    if (fp == NULL) {
        fprintf(stderr, "Error opening batch file \"%s\"\n", bfile);
        exit(1);
    }
    while(fgets(stext, 255, fp) != NULL) {
        if (stext[0] == '#') {
            continue;
        }
        lowtoupper(stext);
/* TODO check on 0, 0 here for index and setno */
        scanner(stext, &x, &y, 1, ax, bx, cx, dx, 1, 0, 0, &errpos);
        stext[0] = 0;
        if (gotparams && paramfile[0]) {
            if (!getparms(cg, paramfile)) {
            }
            gotparams = 0;
        } else if (gotread && readfile[0]) {
            if (getdata(cg, readfile, readsrc, readtype)) {
                drawgraph();
            }
            gotread = 0;
        }
    }
    if (fp != stdin) {
	fclose(fp);
    }
}

symtab_entry key[] = {
"ABORT", ABORT,
"ABOVE", ABOVE,
"ABSOLUTE", ABSOLUTE,
"ACTIVATE", ACTIVATE,
"ACTIVE", ACTIVE,
"ADCIRC", ADCIRC,
"ADCIRC3DFLOW", ADCIRC3DFLOW,
"ALL", ALL,
"ALT", ALT,
"ALTERNATE", ALTERNATE,
"ALTXAXIS", ALTXAXIS,
"ALTYAXIS", ALTYAXIS,
"AMP", AMP,
"ANGLE", ANGLE,
"ANNOTATE", ANNOTATE,
"APPEND", APPEND,
"AREA", AREA,
"ARROW", ARROW,
"ASCEND", ASCEND,
"ASCENDING", ASCEND,
"AT", AT,
"ATTACH", ATTACH,
"AUTO", AUTO,
"AUTOSCALE", AUTOSCALE,
"AUTOTICKS", AUTOTICKS,
"AVERAGE", AVERAGE,
"AVG", AVG,
"AXES", AXES,
"AXIS", AXIS,
"BACKBUFFER", BACKBUFFER,
"BACKGROUND", BACKGROUND,
"BAR", BAR,
"BATCH", BATCH,
"BATH", BATH,
"BATHYMETRY", BATHYMETRY,
"BELOW", BELOW,
"BIN", BIN,
"BINARY", BINARY,
"BOTH", BOTH,
"BOTTOM", BOTTOM,
"BOUNDARY", BOUNDARY,
"BOX", BOX,
"CELLS", CELLS,
"CENTER", CENTER,
"CH3D", CH3D,
"CHAR", CHAR,
"CHDIR", CHDIR,
"CIRCLE", CIRCLE,
"CLEAR", CLEAR,
"CLICK", CLICK,
"CLOCK", CLOCK,
"CLOSE", CLOSE,
"CM", CM,
"CMAP", CMAP,
"COLOR", COLOR,
"COLORMAP", COLORMAP,
"COMMENT", COMMENT,
"CONC", CONC,
"CONCENTRATION", CONCENTRATION,
"CONCENTRATIONS", CONCENTRATIONS,
"COPY", COPY,
"COURANT", COURANT,
"CROSS", CROSS,
"CYCLE", CYCLE,
"DAYMONTH", DAYMONTH,
"DAYOFWEEKL", DAYOFWEEKL,
"DAYOFWEEKS", DAYOFWEEKS,
"DAYOFYEAR", DAYOFYEAR,
"DDMMYY", DDMMYY,
"DECIMAL", DECIMAL,
"DEF", DEF,
"DEFAULT", DEFAULT,
"DEGREESLAT", DEGREESLAT,
"DEGREESLON", DEGREESLON,
"DEGREESMMLAT", DEGREESMMLAT,
"DEGREESMMLON", DEGREESMMLON,
"DEGREESMMSSLAT", DEGREESMMSSLAT,
"DEGREESMMSSLON", DEGREESMMSSLON,
"DELAYP", DELAYP,
"DELETE", DELETE,
"DEPTH", DEPTH,
"DEPTHS", DEPTHS,
"DESCEND", DESCEND,
"DESCENDING", DESCEND,
"DEVICE", DEVICE,
"DEVXY", DEVXY,
"DFT", DFT,
"DIAMOND", DIAMOND,
"DIFFERENCE", DIFFERENCE,
"DISK", DISK,
"DISPLAY", DISPLAY,
"DOT", DOT,
"DOUBLEBUFFER", DOUBLEBUFFER,
"DOWN", DOWN,
"DRAW2", DRAW2,
"DROGUE", DROGUE,
"DROGUES", DROGUES,
"DT", DT,
"DXDX", DXDX,
"DXP", DXP,
"DYDY", DYDY,
"DYP", DYP,
"ECHO", ECHO,
"EDIT", EDIT,
"ELA", ELA,
"ELCIRC", ELCIRC,
"ELEMENT", ELEMENT,
"ELEMENTS", ELEMENTS,
"ELEV", ELEV,
"ELEVATION", ELEVATION,
"ELEVATIONS", ELEVATIONS,
"ELEVMARKER", ELEVMARKER,
"ELSE", ELSE,
"END", END,
"ERRORBAR", ERRORBAR,
"EXIT", EXIT,
"EXPAND", EXPAND,
"EXPONENTIAL", EXPONENTIAL,
"FACTOR", FACTOR,
"FALSE", FALSEP,
"FAST", FAST,
"FEET", FEET,
"FFT", FFT,
"FILE", FILEP,
"FILL", FILL,
"FIND", FIND,
"FIXEDPOINT", FIXEDPOINT,
"FLOW", FLOW,
"FLUSH", FLUSH,
"FLUX", FLUX,
"FOCUS", FOCUS,
"FOLLOWS", FOLLOWS,
"FONT", FONTP,
"FOREGROUND", FOREGROUND,
"FORMAT", FORMAT,
"FORT14", FORT14,
"FORT63", FORT63,
"FORT64", FORT64,
"FORWARD", FORWARD,
"FRAME", FRAMEP,
"FREQ", FREQ,
"FRONTBUFFER", FRONTBUFFER,
"GENERAL", GENERAL,
"GETP", GETP,
"GOTO", GOTO,
"GRAPH", GRAPH,
"GRAPHNO", GRAPHNO,
"GRAPHS", GRAPHS,
"GRAPHTYPE", GRAPHTYPE,
"GRID", GRID,
"HARDCOPY", HARDCOPY,
"HBAR", HBAR,
"HELP", HELP,
"HGAP", HGAP,
"HIDDEN", HIDDEN,
"HISTBOX", HISTBOX,
"HISTO", HISTO,
"HISTORY", HISTORY,
"HMS", HMS,
"HORIZONTAL", HORIZONTAL,
"HPGLL", HPGLL,
"HPGLP", HPGLP,
"IF", IF,
"IGNORE", IGNORE,
"IHL", IHL,
"IMAGE", IMAGE,
"IMAGES", IMAGES,
"IN", IN,
"INCLUDE", INCLUDE,
"INFO", INFO,
"INIT", INIT,
"INITGRAPHICS", INITGRAPHICS,
"INOUT", INOUT,
"INTEGRATE", INTEGRATE,
"INTERP", INTERP,
"INUNDATION", INUNDATION,
"INVDFT", INVDFT,
"INVFFT", INVFFT,
"ISOLINE", ISOLINE,
"ISOLINES", ISOLINES,
"JUST", JUST,
"KILL", KILL,
"KM", KM,
"LABEL", LABEL,
"LAYOUT", LAYOUT,
"LEAVE", LEAVE,
"LEAVEGRAPHICS", LEAVEGRAPHICS,
"LEFT", LEFT,
"LEGEND", LEGEND,
"LENGTH", LENGTH,
"LEVEL", LEVEL,
"LEVELS", LEVELS,
"LIMITS", LIMITS,
"LINE", LINE,
"LINES", LINES,
"LINESTYLE", LINESTYLE,
"LINETO", LINETO,
"LINEW", LINEW,
"LINEWIDTH", LINEWIDTH,
"LINK", LINK,
"LOAD", LOAD,
"LOC", LOC,
"LOCATE", LOCATE,
"LOCATOR", LOCATOR,
"LOCTYPE", LOCTYPE,
"LOG", LOG,
"LOGX", LOGX,
"LOGXY", LOGXY,
"LOGY", LOGY,
"M", M,
"MAG", MAG,
"MAGNITUDE", MAGNITUDE,
"MAJOR", MAJOR,
"MAPSCALE", MAPSCALE,
"MARKER", MARKER,
"MARKERS", MARKERS,
"MAX", MAXP,
"MAXLEVELS", MAXLEVELS,
"METHOD", METHOD,
"MIFL", MIFL,
"MIFP", MIFP,
"MILES", MILES,
"MIN", MINP,
"MINOR", MINOR,
"MISSINGP", MISSINGP,
"MM", MM,
"MMDD", MMDD,
"MMDDHMS", MMDDHMS,
"MMDDYY", MMDDYY,
"MMDDYYHMS", MMDDYYHMS,
"MMSSLAT", MMSSLAT,
"MMSSLON", MMSSLON,
"MMYY", MMYY,
"MONTHDAY", MONTHDAY,
"MONTHL", MONTHL,
"MONTHS", MONTHS,
"MOVE", MOVE,
"MOVE2", MOVE2,
"MOVETO", MOVETO,
"NEGATE", NEGATE,
"NO", NO,
"NODE", NODE,
"NODES", NODES,
"NONE", NONE,
"NORMAL", NORMAL,
"NORTH", NORTH,
"NXY", NXY,
"OFF", OFF,
"OFFSETX", OFFSETX,
"OFFSETY", OFFSETY,
"ON", ON,
"OP", OP,
"OPEN", OPEN,
"ORIENT", ORIENT,
"OUT", OUT,
"PAGE", PAGE,
"PARA", PARA,
"PARALLEL", PARALLEL,
"PARAMETERS", PARAMETERS,
"PARAMS", PARAMS,
"PARMS", PARMS,
"PATTERN", PATTERN,
"PER", PER,
"PERIMETER", PERIMETER,
"PERP", PERP,
"PERPENDICULAR", PERPENDICULAR,
"PHASE", PHASE,
"PIE", PIE,
"PIPE", PIPE,
"PLACE", PLACE,
"PLUS", PLUS,
"POINT", POINT,
"POLAR", POLAR,
"POLY", POLY,
"POLYI", POLYI,
"POLYO", POLYO,
"POP", POP,
"POWER", POWER,
"PREC", PREC,
"PREFIX", PREFIX,
"PREPEND", PREPEND,
"PRINT", PRINT,
"PROP", PROP,
"PS", PS,
"PSCOLORL", PSCOLORL,
"PSCOLORP", PSCOLORP,
"PSMONOL", PSMONOL,
"PSMONOP", PSMONOP,
"PUSH", PUSH,
"PUTP", PUTP,
"QUIT", QUIT,
"READ", READ,
"READBIN", READBIN,
"REDRAW", REDRAW,
"REGION", REGION,
"REGIONS", REGIONS,
"REGNUM", REGNUM,
"REGRESS", REGRESS,
"REMOVE", REMOVE,
"RENDER", RENDER,
"REPORT", REPORT,
"RESET", RESET,
"REVERSE", REVERSE,
"REWIND", REWIND,
"RIGHT", RIGHT,
"RISER", RISER,
"ROT", ROT,
"RUN", RUN,
"SALINITY", SALINITY,
"SAMPLE", SAMPLE,
"SAVE", SAVE,
"SCALAR", SCALAR,
"SCALE", SCALE,
"SCIENTIFIC", SCIENTIFIC,
"SET", SET,
"SETS", SETS,
"SHOW", SHOW,
"SHRINK", SHRINK,
"SIGMA", SIGMA,
"SIGN", SIGN,
"SIZE", SIZE,
"SKIP", SKIP,
"SLEEP", SLEEP,
"SLICE", SLICE,
"SOURCE", SOURCE,
"SPEC", SPEC,
"SPECIFIED", SPECIFIED,
"SPECTRUM", SPECTRUM,
"SPLITS", SPLITS,
"SQUARE", SQUARE,
"STACK", STACK,
"STACKEDBAR", STACKEDBAR,
"STACKEDHBAR", STACKEDHBAR,
"STACKEDLINE", STACKEDLINE,
"STAGGER", STAGGER,
"STAR", STAR,
"START", START,
"STARTSTEP", STARTSTEP,
"STARTTYPE", STARTTYPE,
"STATUS", STATUS,
"STEP", STEP,
"STOP", STOP,
"STREAMLINES", STREAMLINES,
"STRING", STRING,
"STRINGS", STRINGS,
"SUBTITLE", SUBTITLE,
"SURFACE", SURFACE,
"SWAPBUFFER", SWAPBUFFER,
"SYMBOL", SYMBOL,
"SYSTEM", SYSTEM,
"TEANL", TEANL,
"TEXT", TEXT,
"TICK", TICKP,
"TICKLABEL", TICKLABEL,
"TICKMARKS", TICKMARKS,
"TIDALCLOCK", TIDALCLOCK,
"TIME", TIME,
"TIMEINFO", TIMEINFO,
"TIMELINE", TIMELINE,
"TITLE", TITLE,
"TO", TO,
"TOP", TOP,
"TOTAL", TOTAL,
"TRACK", TRACK,
"TRANSECT", TRANSECT,
"TRIANGLE1", TRIANGLE1,
"TRIANGLE2", TRIANGLE2,
"TRIANGLE3", TRIANGLE3,
"TRIANGLE4", TRIANGLE4,
"TRUE", TRUEP,
"TYPE", TYPE,
"UNITS", UNITS,
"UP", UP,
"VALUE", VALUE,
"VECTOR", VECTOR,
"VEL", VEL,
"VELOCITY", VELOCITY,
"VERTICAL", VERTICAL,
"VGAP", VGAP,
"VIEW", VIEW,
"VSCALE", VSCALE,
"VX1", VX1,
"VX2", VX2,
"VY1", VY1,
"VY2", VY2,
"WIDTH", WIDTH,
"WIND", WIND,
"WITH", WITH,
"WORLD", WORLD,
"WRAP", WRAP,
"WRITE", WRITE,
"WSCALE", WSCALE,
"WX1", WX1,
"WX2", WX2,
"WY1", WY1,
"WY2", WY2,
"X1", X1,
"X2", X2,
"X3", X3,
"X4", X4,
"X5", X5,
"XAXES", XAXES,
"XAXIS", XAXIS,
"XCOR", XCOR,
"XMAX", XMAX,
"XMIN", XMIN,
"XY", XY,
"XYARC", XYARC,
"XYBOX", XYBOX,
"XYDX", XYDX,
"XYDXDX", XYDXDX,
"XYDXDY", XYDXDY,
"XYDY", XYDY,
"XYDYDY", XYDYDY,
"XYFIXED", XYFIXED,
"XYHILO", XYHILO,
"XYRT", XYRT,
"XYSEG", XYSEG,
"XYSTRING", XYSTRING,
"XYUV", XYUV,
"XYX2Y2", XYX2Y2,
"XYXX", XYXX,
"XYYY", XYYY,
"XYZ", XYZ,
"XYZW", XYZW,
"Y1", Y1,
"Y2", Y2,
"Y3", Y3,
"Y4", Y4,
"Y5", Y5,
"YAXES", YAXES,
"YAXIS", YAXIS,
"YES", YES,
"YMAX", YMAX,
"YMIN", YMIN,
"ZEROXAXIS", ZEROXAXIS,
"ZEROYAXIS", ZEROYAXIS,
"ZOOM", ZOOM,
"ZOOMBOX", ZOOMBOX
};

int maxparms = sizeof(key) / sizeof(symtab_entry);
int maxfunc = sizeof(key) / sizeof(symtab_entry);

int findf(symtab_entry *key, char *s, int tlen)
{

    int low, high, mid;

    low = 0;
    high = tlen - 1;
    while (low <= high) {
	mid = (low + high) / 2;
	if (strcmp(s, key[mid].s) < 0) {
	    high = mid - 1;
	} else {
	    if (strcmp(s, key[mid].s) > 0) {
		low = mid + 1;
	    } else {
		return (mid);
	    }
	}
    }
    return (-1);
}

int getcharstr(void)
{
    if (pos >= strlen(f_string))
	 return EOF;
    return (f_string[pos++]);
}

void ungetchstr(void)
{
    if (pos > 0)
	pos--;
}

int yylex(void)
{
    int c, i;
    int found;
    static char s[256];
    char sbuf[256];

    while ((c = getcharstr()) == ' ' || c == '\t');
    if (c == EOF) {
	return (0);
    }
    if (c == '"') {
	i = 0;
	while ((c = getcharstr()) != '"' && c != EOF) {
	    if (c == '\\') {
		int ctmp;
		ctmp = getcharstr();
		if (ctmp != '"') {
		    ungetchstr();
		}
		else {
		    c = ctmp;
		}
	    }
	    s[i] = c;
	    i++;
	}
	if (c == EOF) {
	    sprintf(sbuf, "Nonterminating string\n");
	    yyerror(sbuf);
	    return 0;
	}
	s[i] = '\0';
	yylval.str = s;
	return CHRSTR;
    }
    if (c == '.' || isdigit(c)) {
	char stmp[80];
	double d;
	int i, gotdot = 0;

	i = 0;
	while (c == '.' || isdigit(c)) {
	    if (c == '.') {
		if (gotdot) {
		    yyerror("Reading number, too many dots");
	    	    return 0;
		} else {
		    gotdot = 1;
		}
	    }
	    stmp[i++] = c;
	    c = getcharstr();
	}
	if (c == 'E' || c == 'e') {
	    stmp[i++] = c;
	    c = getcharstr();
	    if (c == '+' || c == '-') {
		stmp[i++] = c;
		c = getcharstr();
	    }
	    while (isdigit(c)) {
		stmp[i++] = c;
		c = getcharstr();
	    }
	}
	if (gotdot && i == 1) {
	    ungetchstr();
	    return '.';
	}
	stmp[i] = '\0';
	ungetchstr();
	sscanf(stmp, "%lf", &d);
	yylval.val = d;
	return NUMBER;
    }
/* graphs, sets, regions resp. */
    if (c == 'G' || c == 'S' || c == 'R') {
	char stmp[80];
	double d;
	int i = 0, ctmp = c, gn, sn, rn;
	c = getcharstr();
	while (isdigit(c)) {
	    stmp[i++] = c;
	    c = getcharstr();
	}
	if (i == 0) {
	    c = ctmp;
	    ungetchstr();
	} else {
	    ungetchstr();
	    if (ctmp == 'G') {
	        stmp[i] = '\0';
		gn = atoi(stmp);
		if (gn >= 0 && gn < maxgraph) {
		    yylval.ival = gn;
		    whichgraph = gn;
		    return GRAPHNO;
		}
	    }
	}
    }
    if (isalpha(c)) {
	char *p = sbuf;
	int gno = -1, setno = -1, xy = -1, elno = -1;

	do {
	    *p++ = c;
	} while ((c = getcharstr()) != EOF && isalnum(c));
	ungetchstr();
	*p = '\0';
        if (debuglevel == 2) {
	    printf("->%s<-\n", sbuf);
	}
	if ((found = findf(key, sbuf, maxfunc)) >= 0) {
	    if (key[found].type == VAR) {
		switch (sbuf[0]) {
		case 'A':
		    yylval.ptr = aa;
		    return VAR;
		case 'B':
		    yylval.ptr = bb;
		    return VAR;
		case 'C':
		    yylval.ptr = cc;
		    return VAR;
		case 'D':
		    yylval.ptr = dd;
		    return VAR;
		}
	    }
	    else if (key[found].type == FITPARM) {
		int index = sbuf[1] - '0';
		yylval.val = nonl_parms[index];
		return FITPARM;
	    }
	    else { /* set up special cases */
		switch (key[found].type) {
		case XAXIS:
		    naxis = 0;
		    break;
		case YAXIS:
		    naxis = 1;
		    break;
		case ZEROXAXIS:
		    naxis = 2;
		    break;
		case ZEROYAXIS:
		    naxis = 3;
		    break;
		case ALTXAXIS:
		    naxis = 4;
		    break;
		case ALTYAXIS:
		    naxis = 5;
		    break;
		case AXES:
		    naxis = 6;
		    break;
		case XAXES:
		    naxis = 7;
		    break;
		case YAXES:
		    naxis = 8;
		    break;
		case GRAPHS:
		    yylval.ival = -1;
		    whichgraph = -1;
		    return GRAPHS;
		case SETS:
		    yylval.ival = -1;
		    whichset = -1;
		    return SETS;
		default:
		    break;
		}
	    }
	    yylval.func = key[found].type;
	    return key[found].type;
	} else {
	    strcat(sbuf, ": No such function or variable");
	    yyerror(sbuf);
	    return 0;
	}
    }
    switch (c) {
    case '>':
	return follow('=', GE, GT);
    case '<':
	return follow('=', LE, LT);
    case '=':
	return follow('=', EQ, '=');
    case '!':
	return follow('=', NE, NOT);
    case '|':
	return follow('|', OR, '|');
    case '&':
	return follow('&', AND, '&');
    case '\n':
	return '\n';
    default:
	return c;
    }
}

int follow(int expect, int ifyes, int ifno)
{
    int c = getcharstr();

    if (c == expect) {
	return ifyes;
    }
    ungetchstr();
    return ifno;
}

void yyerror(char *s)
{
    int i;
    char buf[256];
    sprintf(buf, "Error: %s: %s", s, f_string);
    i = strlen(buf);
    buf[i - 1] = 0;
    errwin(buf);
    interr = 1;
}

#define C1 0.1978977093962766
#define C2 0.1352915131768107

double rnorm(double mean, double sdev)
{
    double u = drand48();

    return mean + sdev * (pow(u, C2) - pow(1.0 - u, C2)) / C1;
}

double fx(double x)
{
    return 1.0 / sqrt(2.0 * M_PI) * exp(-x * x * 0.5);
}

double normp(double b, double *s)
{
    double sum, dx, a = -8.0, fx(double x);
    int i, n = 48;

    sum = fx(a) + fx(b);
    dx = (b - a) / n;
    for (i = 1; i <= ((n - 1) / 2); i++)
	sum = sum + 4.0 * fx(a + (2.0 * i - 1.0) * dx) + 2.0 * fx(a + 2.0 * i * dx);
    sum = sum + 4.0 * fx(b - dx);
    *s = fx(b);
    return sum * dx / 3.0;
}
