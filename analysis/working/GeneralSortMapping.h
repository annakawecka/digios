/***********************************************************************
 * 
 * This is the Mapping for Armory/GeneralSort.C
 * 
 * idDetMap = 0XX ; PSD
 *          = 1XX ; Recoil
 *          = 2XX ; ELUM
 *          = 3XX ; EZERO
 *          = 4XX ; TAC & RF-Timing
 * 
***********************************************************************/

#define NARRAY  40
#define NRDT    8
#define NELUM   2
#define NEZERO  5
#define NTAC    3

#define MWIN 100 //M value for energy filter from digi setting


Int_t idConst = 1010;

Int_t idDetMap[160] = {100, 101, 102, 103, 104, 105, 106, 107,  -1,  -1,   //VME1-DIG1-Recoil
                       400, 401, 402,  -1,  -1,  -1,  -1,  -1,  -1,  -1,   //VME1-DIG2-TAC
                       300, 301,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,   //VME1-DIG3-IC
                        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,   //VME1-DIG4
                   
                         1,   0,   5,   4,   3,   2,   1,   0,  -1,  -1,/*1*/
                         3,   2,   1,   0,   5,   4,   3,   2,  -1,  -1,/*2*/
                       240,  10,   9,   8,   7,   6,   5,   4,  -1,  -1,/*3*/
                        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,//
                        
                         7,   6,  11,  10,   9,   8,   7,   6,  -1,  -1,/*4*/
                        15,  14,  13,  12,  11,  10,   9,   8,  -1,  -1,/*5*/
                        17,  16,  15,  14,  13,  12,  17,  16,  -1,  -1,/*6*/
                        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,//
                        
                        19,  18,  17,  16,  15,  14,  13,  12,  -1,  -1,/*7*/
                        21,  20,  19,  18,  23,  22,  21,  20,  -1,  -1,/*8*/
                        23,  22,  21,  20,  19,  18,  23,  22,  -1,  -1,/*9*/      
                        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1};//


Int_t idKindMap[160] = {-1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
                        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
                        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
                        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
                         1,   1,   0,   0,   0,   0,   0,   0,  -1,  -1,//1
                         2,   2,   2,   2,   1,   1,   1,   1,  -1,  -1,//2
                         0,   0,   0,   0,   0,   0,   2,   2,  -1,  -1,//3
                        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
                         2,   2,   1,   1,   1,   1,   1,   1,  -1,  -1,//4
                         0,   0,   0,   0,   2,   2,   2,   2,  -1,  -1,//5
                         2,   2,   2,   2,   2,   2,   0,   0,  -1,  -1,//6
                        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
                         0,   0,   1,   1,   1,   1,   1,   1,  -1,  -1,//7
                         1,   1,   1,   1,   0,   0,   0,   0,  -1,  -1,//8
                         2,   2,   2,   2,   2,   2,   1,   1,  -1,  -1,//9
                        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1};
