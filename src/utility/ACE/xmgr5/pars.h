<<<<<<< .mine
/* A Bison parser, made by GNU Bison 2.4.2.  */
=======
>>>>>>> .r440

<<<<<<< .mine
/* Skeleton interface for Bison's Yacc-like parsers in C
   
      Copyright (C) 1984, 1989-1990, 2000-2006, 2009-2010 Free Software
   Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
=======
/* A Bison parser, made by GNU Bison 2.4.1.  */

/* Skeleton interface for Bison's Yacc-like parsers in C
   
      Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
>>>>>>> .r440
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


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     ABS = 258,
     ACOS = 259,
     ALL = 260,
     ASIN = 261,
     ATAN = 262,
     ATAN2 = 263,
     CEIL = 264,
     COS = 265,
     DEG = 266,
     DX = 267,
     DY = 268,
     ERF = 269,
     ERFC = 270,
     EXP = 271,
     FLOOR = 272,
     HYPOT = 273,
     INDEX = 274,
     INT = 275,
     INVN = 276,
     INVT = 277,
     IRAND = 278,
     LGAMMA = 279,
     LN = 280,
     LOG = 281,
     LOGISTIC = 282,
     MAXP = 283,
     MINP = 284,
     MOD = 285,
     NORM = 286,
     NORMP = 287,
     PI = 288,
     RAD = 289,
     RAND = 290,
     RNORM = 291,
     SETNO = 292,
     SIN = 293,
     SQR = 294,
     SQRT = 295,
     TAN = 296,
     INUM = 297,
     CTD = 298,
     ADP = 299,
     TPC = 300,
     VX1 = 301,
     VX2 = 302,
     VY1 = 303,
     VY2 = 304,
     WX1 = 305,
     WX2 = 306,
     WY1 = 307,
     WY2 = 308,
     DELAYP = 309,
     DOUBLEBUFFER = 310,
     DOWN = 311,
     ABSOLUTE = 312,
     ABORT = 313,
     ACTIVATE = 314,
     ACTIVE = 315,
     ALT = 316,
     ALTERNATE = 317,
     ALTXAXIS = 318,
     ALTYAXIS = 319,
     ANGLE = 320,
     ANNOTATE = 321,
     APPEND = 322,
     AREA = 323,
     ARROW = 324,
     AUTO = 325,
     AUTOSCALE = 326,
     AUTOTICKS = 327,
     AVG = 328,
     AXIS = 329,
     AXES = 330,
     BACKBUFFER = 331,
     BACKGROUND = 332,
     BAR = 333,
     BATCH = 334,
     BLOCK = 335,
     BIN = 336,
     BOTH = 337,
     BOTTOM = 338,
     BOX = 339,
     BOXPLOT = 340,
     BP = 341,
     CD = 342,
     CELLS = 343,
     CENTER = 344,
     CHAR = 345,
     CHRSTR = 346,
     CLEAR = 347,
     CLICK = 348,
     CMAP = 349,
     COLOR = 350,
     COMMENT = 351,
     COPY = 352,
     CORIE = 353,
     CYCLE = 354,
     DB = 355,
     DECIMAL = 356,
     DEF = 357,
     DEFAULT = 358,
     DELETE = 359,
     DEVICE = 360,
     DFT = 361,
     DIFFERENCE = 362,
     DISK = 363,
     DRAW2 = 364,
     DROP = 365,
     DXDX = 366,
     DXP = 367,
     DYDY = 368,
     DYP = 369,
     ECHO = 370,
     EDIT = 371,
     ELSE = 372,
     END = 373,
     ERRORBAR = 374,
     EXIT = 375,
     EXPONENTIAL = 376,
     FALSEP = 377,
     FFT = 378,
     FILEP = 379,
     FILL = 380,
     FIND = 381,
     FIXEDPOINT = 382,
     FLUSH = 383,
     FOCUS = 384,
     FOLLOWS = 385,
     FONTP = 386,
     FOREGROUND = 387,
     FORMAT = 388,
     FRONTBUFFER = 389,
     FRAMEP = 390,
     GETP = 391,
     GIFL = 392,
     GIFP = 393,
     GRAPH = 394,
     GRAPHNO = 395,
     GRAPHS = 396,
     GRAPHTYPE = 397,
     GRID = 398,
     HARDCOPY = 399,
     HBAR = 400,
     HBOXPLOT = 401,
     HGAP = 402,
     HIDDEN = 403,
     HORIZONTAL = 404,
     HISTO = 405,
     IF = 406,
     IGNORE = 407,
     IHL = 408,
     IMAGE = 409,
     IN = 410,
     INIT = 411,
     INITGRAPHICS = 412,
     INOUT = 413,
     INTEGRATE = 414,
     INTERP = 415,
     INVDFT = 416,
     INVFFT = 417,
     ISOLINE = 418,
     ISOLINES = 419,
     JUST = 420,
     KILL = 421,
     LABEL = 422,
     LAYOUT = 423,
     LEAVE = 424,
     LEAVEGRAPHICS = 425,
     LEFT = 426,
     LEGEND = 427,
     LENGTH = 428,
     LEVEL = 429,
     LEVELS = 430,
     LINE = 431,
     LINESTYLE = 432,
     LINETO = 433,
     LINEWIDTH = 434,
     LINK = 435,
     LOAD = 436,
     LOCATOR = 437,
     LOCTYPE = 438,
     LOGX = 439,
     LOGY = 440,
     LOGXY = 441,
     MAJOR = 442,
     MIFL = 443,
     MIFP = 444,
     MINOR = 445,
     MISSINGP = 446,
     MONITOR = 447,
     MOVE = 448,
     MOVE2 = 449,
     MOVETO = 450,
     NEGATE = 451,
     NO = 452,
     NONE = 453,
     NORMAL = 454,
     NXY = 455,
     OFF = 456,
     OFFSETX = 457,
     OFFSETY = 458,
     ON = 459,
     OP = 460,
     ORIENT = 461,
     OUT = 462,
     PAGE = 463,
     PARA = 464,
     PARALLEL = 465,
     PARAMETERS = 466,
     PARAMS = 467,
     PATTERN = 468,
     PERIMETER = 469,
     PERP = 470,
     PERPENDICULAR = 471,
     PIE = 472,
     PIPE = 473,
     PLACE = 474,
     POINT = 475,
     POLAR = 476,
     POWER = 477,
     PREC = 478,
     PREPEND = 479,
     PRINT = 480,
     PS = 481,
     PSCOLORP = 482,
     PSMONOP = 483,
     PSCOLORL = 484,
     PSMONOL = 485,
     PUSH = 486,
     POP = 487,
     PUTP = 488,
     READ = 489,
     REDRAW = 490,
     REGRESS = 491,
     REGNUM = 492,
     REGIONS = 493,
     RENDER = 494,
     REVERSE = 495,
     RIGHT = 496,
     RISER = 497,
     ROT = 498,
     RUNAVG = 499,
     RUNMED = 500,
     RUNSTD = 501,
     RUNMIN = 502,
     RUNMAX = 503,
     SAMPLE = 504,
     SAVEALL = 505,
     SCALE = 506,
     SCIENTIFIC = 507,
     SET = 508,
     SETNUM = 509,
     SETS = 510,
     SIGN = 511,
     SIZE = 512,
     SKIP = 513,
     SLEEP = 514,
     SLICE = 515,
     SOURCE = 516,
     SPEC = 517,
     SPECIFIED = 518,
     SPECTRUM = 519,
     STACK = 520,
     STACKEDBAR = 521,
     STACKEDHBAR = 522,
     STACKEDLINE = 523,
     STAGGER = 524,
     START = 525,
     STARTTYPE = 526,
     STATUS = 527,
     STEP = 528,
     STOP = 529,
     STRING = 530,
     SUBTITLE = 531,
     SWAPBUFFER = 532,
     SYMBOL = 533,
     TICKP = 534,
     TICKLABEL = 535,
     TICKMARKS = 536,
     TITLE = 537,
     TO = 538,
     TOP = 539,
     TRUEP = 540,
     TYPE = 541,
     UP = 542,
     VELOCITY = 543,
     VERTICAL = 544,
     VGAP = 545,
     VIEW = 546,
     WITH = 547,
     WORLD = 548,
     WRITE = 549,
     X = 550,
     X0 = 551,
     X1 = 552,
     XAXES = 553,
     XAXIS = 554,
     XCOR = 555,
     XMAX = 556,
     XMIN = 557,
     FEGRID = 558,
     RECTGRID = 559,
     XY = 560,
     XYARC = 561,
     XYBOX = 562,
     XYBOXPLOT = 563,
     XYFIXED = 564,
     XYHILO = 565,
     XYRT = 566,
     XYSEG = 567,
     XYSTRING = 568,
     XYDX = 569,
     XYDY = 570,
     XYDXDX = 571,
     XYDYDY = 572,
     XYDXDY = 573,
     XYX2Y2 = 574,
     XYXX = 575,
     XYYY = 576,
     XYZ = 577,
     XYZW = 578,
     XYUV = 579,
     Y = 580,
     Y0 = 581,
     Y1 = 582,
     Y2 = 583,
     Y3 = 584,
     Y4 = 585,
     Y5 = 586,
     YAXES = 587,
     YAXIS = 588,
     YES = 589,
     YMAX = 590,
     YMIN = 591,
     ZEROXAXIS = 592,
     ZEROYAXIS = 593,
     ABOVE = 594,
     BELOW = 595,
     POLYI = 596,
     POLYO = 597,
     GENERAL = 598,
     DDMMYY = 599,
     MMDDYY = 600,
     MMYY = 601,
     MMDD = 602,
     MONTHDAY = 603,
     DAYMONTH = 604,
     MONTHS = 605,
     MONTHL = 606,
     DDMONTHSYYHHMMSS = 607,
     DDMONTHSYY = 608,
     DAYOFWEEKS = 609,
     DAYOFWEEKL = 610,
     DAYOFYEAR = 611,
     HMS = 612,
     HH = 613,
     MMDDHMS = 614,
     MMDDYYHMS = 615,
     DEGREESLON = 616,
     DEGREESMMLON = 617,
     DEGREESMMSSLON = 618,
     MMSSLON = 619,
     DEGREESLAT = 620,
     DEGREESMMLAT = 621,
     DEGREESMMSSLAT = 622,
     MMSSLAT = 623,
     DOT = 624,
     STAR = 625,
     PLUS = 626,
     CROSS = 627,
     CIRCLE = 628,
     SQUARE = 629,
     DIAMOND = 630,
     TRIANGLE1 = 631,
     TRIANGLE2 = 632,
     TRIANGLE3 = 633,
     TRIANGLE4 = 634,
     SPLINE = 635,
     LANDSCAPE = 636,
     PORTRAIT = 637,
     FREE = 638,
     FIXED = 639,
     STATUSBAR = 640,
     LOCATORBAR = 641,
     TOOLBAR = 642,
     ELCIRC = 643,
     SCALAR = 644,
     VECTOR = 645,
     HEAT = 646,
     HISTORY = 647,
     PROFILE = 648,
     NODE = 649,
     VAR = 650,
     NUMBER = 651,
     OR = 652,
     AND = 653,
     NE = 654,
     EQ = 655,
     GE = 656,
     LE = 657,
     LT = 658,
     GT = 659,
     NOT = 660,
     UMINUS = 661
   };
#endif
/* Tokens.  */
#define ABS 258
#define ACOS 259
#define ALL 260
#define ASIN 261
#define ATAN 262
#define ATAN2 263
#define CEIL 264
#define COS 265
#define DEG 266
#define DX 267
#define DY 268
#define ERF 269
#define ERFC 270
#define EXP 271
#define FLOOR 272
#define HYPOT 273
#define INDEX 274
#define INT 275
#define INVN 276
#define INVT 277
#define IRAND 278
#define LGAMMA 279
#define LN 280
#define LOG 281
#define LOGISTIC 282
#define MAXP 283
#define MINP 284
#define MOD 285
#define NORM 286
#define NORMP 287
#define PI 288
#define RAD 289
#define RAND 290
#define RNORM 291
#define SETNO 292
#define SIN 293
#define SQR 294
#define SQRT 295
#define TAN 296
#define INUM 297
#define CTD 298
#define ADP 299
#define TPC 300
#define VX1 301
#define VX2 302
#define VY1 303
#define VY2 304
#define WX1 305
#define WX2 306
#define WY1 307
#define WY2 308
#define DELAYP 309
#define DOUBLEBUFFER 310
#define DOWN 311
#define ABSOLUTE 312
#define ABORT 313
#define ACTIVATE 314
#define ACTIVE 315
#define ALT 316
#define ALTERNATE 317
#define ALTXAXIS 318
#define ALTYAXIS 319
#define ANGLE 320
#define ANNOTATE 321
#define APPEND 322
#define AREA 323
#define ARROW 324
#define AUTO 325
#define AUTOSCALE 326
#define AUTOTICKS 327
#define AVG 328
#define AXIS 329
#define AXES 330
#define BACKBUFFER 331
#define BACKGROUND 332
#define BAR 333
#define BATCH 334
#define BLOCK 335
#define BIN 336
#define BOTH 337
#define BOTTOM 338
#define BOX 339
#define BOXPLOT 340
#define BP 341
#define CD 342
#define CELLS 343
#define CENTER 344
#define CHAR 345
#define CHRSTR 346
#define CLEAR 347
#define CLICK 348
#define CMAP 349
#define COLOR 350
#define COMMENT 351
#define COPY 352
#define CORIE 353
#define CYCLE 354
#define DB 355
#define DECIMAL 356
#define DEF 357
#define DEFAULT 358
#define DELETE 359
#define DEVICE 360
#define DFT 361
#define DIFFERENCE 362
#define DISK 363
#define DRAW2 364
#define DROP 365
#define DXDX 366
#define DXP 367
#define DYDY 368
#define DYP 369
#define ECHO 370
#define EDIT 371
#define ELSE 372
#define END 373
#define ERRORBAR 374
#define EXIT 375
#define EXPONENTIAL 376
#define FALSEP 377
#define FFT 378
#define FILEP 379
#define FILL 380
#define FIND 381
#define FIXEDPOINT 382
#define FLUSH 383
#define FOCUS 384
#define FOLLOWS 385
#define FONTP 386
#define FOREGROUND 387
#define FORMAT 388
#define FRONTBUFFER 389
#define FRAMEP 390
#define GETP 391
#define GIFL 392
#define GIFP 393
#define GRAPH 394
#define GRAPHNO 395
#define GRAPHS 396
#define GRAPHTYPE 397
#define GRID 398
#define HARDCOPY 399
#define HBAR 400
#define HBOXPLOT 401
#define HGAP 402
#define HIDDEN 403
#define HORIZONTAL 404
#define HISTO 405
#define IF 406
#define IGNORE 407
#define IHL 408
#define IMAGE 409
#define IN 410
#define INIT 411
#define INITGRAPHICS 412
#define INOUT 413
#define INTEGRATE 414
#define INTERP 415
#define INVDFT 416
#define INVFFT 417
#define ISOLINE 418
#define ISOLINES 419
#define JUST 420
#define KILL 421
#define LABEL 422
#define LAYOUT 423
#define LEAVE 424
#define LEAVEGRAPHICS 425
#define LEFT 426
#define LEGEND 427
#define LENGTH 428
#define LEVEL 429
#define LEVELS 430
#define LINE 431
#define LINESTYLE 432
#define LINETO 433
#define LINEWIDTH 434
#define LINK 435
#define LOAD 436
#define LOCATOR 437
#define LOCTYPE 438
#define LOGX 439
#define LOGY 440
#define LOGXY 441
#define MAJOR 442
#define MIFL 443
#define MIFP 444
#define MINOR 445
#define MISSINGP 446
#define MONITOR 447
#define MOVE 448
#define MOVE2 449
#define MOVETO 450
#define NEGATE 451
#define NO 452
#define NONE 453
#define NORMAL 454
#define NXY 455
#define OFF 456
#define OFFSETX 457
#define OFFSETY 458
#define ON 459
#define OP 460
#define ORIENT 461
#define OUT 462
#define PAGE 463
#define PARA 464
#define PARALLEL 465
#define PARAMETERS 466
#define PARAMS 467
#define PATTERN 468
#define PERIMETER 469
#define PERP 470
#define PERPENDICULAR 471
#define PIE 472
#define PIPE 473
#define PLACE 474
#define POINT 475
#define POLAR 476
#define POWER 477
#define PREC 478
#define PREPEND 479
#define PRINT 480
#define PS 481
#define PSCOLORP 482
#define PSMONOP 483
#define PSCOLORL 484
#define PSMONOL 485
#define PUSH 486
#define POP 487
#define PUTP 488
#define READ 489
#define REDRAW 490
#define REGRESS 491
#define REGNUM 492
#define REGIONS 493
#define RENDER 494
#define REVERSE 495
#define RIGHT 496
#define RISER 497
#define ROT 498
#define RUNAVG 499
#define RUNMED 500
#define RUNSTD 501
#define RUNMIN 502
#define RUNMAX 503
#define SAMPLE 504
#define SAVEALL 505
#define SCALE 506
#define SCIENTIFIC 507
#define SET 508
#define SETNUM 509
#define SETS 510
#define SIGN 511
#define SIZE 512
#define SKIP 513
#define SLEEP 514
#define SLICE 515
#define SOURCE 516
#define SPEC 517
#define SPECIFIED 518
#define SPECTRUM 519
#define STACK 520
#define STACKEDBAR 521
#define STACKEDHBAR 522
#define STACKEDLINE 523
#define STAGGER 524
#define START 525
#define STARTTYPE 526
#define STATUS 527
#define STEP 528
#define STOP 529
#define STRING 530
#define SUBTITLE 531
#define SWAPBUFFER 532
#define SYMBOL 533
#define TICKP 534
#define TICKLABEL 535
#define TICKMARKS 536
#define TITLE 537
#define TO 538
#define TOP 539
#define TRUEP 540
#define TYPE 541
#define UP 542
#define VELOCITY 543
#define VERTICAL 544
#define VGAP 545
#define VIEW 546
#define WITH 547
#define WORLD 548
#define WRITE 549
#define X 550
#define X0 551
#define X1 552
#define XAXES 553
#define XAXIS 554
#define XCOR 555
#define XMAX 556
#define XMIN 557
#define FEGRID 558
#define RECTGRID 559
#define XY 560
#define XYARC 561
#define XYBOX 562
#define XYBOXPLOT 563
#define XYFIXED 564
#define XYHILO 565
#define XYRT 566
#define XYSEG 567
#define XYSTRING 568
#define XYDX 569
#define XYDY 570
#define XYDXDX 571
#define XYDYDY 572
#define XYDXDY 573
#define XYX2Y2 574
#define XYXX 575
#define XYYY 576
#define XYZ 577
#define XYZW 578
#define XYUV 579
#define Y 580
#define Y0 581
#define Y1 582
#define Y2 583
#define Y3 584
#define Y4 585
#define Y5 586
#define YAXES 587
#define YAXIS 588
#define YES 589
#define YMAX 590
#define YMIN 591
#define ZEROXAXIS 592
#define ZEROYAXIS 593
#define ABOVE 594
#define BELOW 595
#define POLYI 596
#define POLYO 597
#define GENERAL 598
#define DDMMYY 599
#define MMDDYY 600
#define MMYY 601
#define MMDD 602
#define MONTHDAY 603
#define DAYMONTH 604
#define MONTHS 605
#define MONTHL 606
#define DDMONTHSYYHHMMSS 607
#define DDMONTHSYY 608
#define DAYOFWEEKS 609
#define DAYOFWEEKL 610
#define DAYOFYEAR 611
#define HMS 612
#define HH 613
#define MMDDHMS 614
#define MMDDYYHMS 615
#define DEGREESLON 616
#define DEGREESMMLON 617
#define DEGREESMMSSLON 618
#define MMSSLON 619
#define DEGREESLAT 620
#define DEGREESMMLAT 621
#define DEGREESMMSSLAT 622
#define MMSSLAT 623
#define DOT 624
#define STAR 625
#define PLUS 626
#define CROSS 627
#define CIRCLE 628
#define SQUARE 629
#define DIAMOND 630
#define TRIANGLE1 631
#define TRIANGLE2 632
#define TRIANGLE3 633
#define TRIANGLE4 634
#define SPLINE 635
#define LANDSCAPE 636
#define PORTRAIT 637
#define FREE 638
#define FIXED 639
#define STATUSBAR 640
#define LOCATORBAR 641
#define TOOLBAR 642
#define ELCIRC 643
#define SCALAR 644
#define VECTOR 645
#define HEAT 646
#define HISTORY 647
#define PROFILE 648
#define NODE 649
#define VAR 650
#define NUMBER 651
#define OR 652
#define AND 653
#define NE 654
#define EQ 655
#define GE 656
#define LE 657
#define LT 658
#define GT 659
#define NOT 660
#define UMINUS 661




<<<<<<< .mine
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 1685 of yacc.c  */
=======
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 1676 of yacc.c  */
>>>>>>> .r440
#line 77 "pars.y"

    double val;
    long ival;
    double *ptr;
    long func;
    long pset;
    char *str;
<<<<<<< .mine



/* Line 1685 of yacc.c  */
#line 874 "y.tab.h"
=======



/* Line 1676 of yacc.c  */
#line 875 "y.tab.h"
>>>>>>> .r440
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

extern YYSTYPE yylval;


