/***********************************************************************
 * 
 * This is the Mapping for Armory/GeneralSort.C
 * 
 * idDetMap = 0XX ; PSD
 *          = 1XX ; Recoil,  dE of the Recoil MUST be odd
 *          = 2XX ; ELUM
 *          = 3XX ; EZERO
 *          = 4XX ; TAC & RF-Timing
 *          = 5XX ; Circular Recoil
 * 
***********************************************************************/

#define NARRAY  24
#define NRDT    8
#define NELUM   0
#define NEZERO  0
#define NTAC    0
#define NCRDT   0

#define MWIN 100 //M value for energy filter from digi setting

Int_t idConst = 1010;

Int_t idDetMap[160] = { -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 400,   //VME1-DIG1-TAC
                       100, 101, 102, 103, 104, 105, 106, 107,  -1,  -1,   //VME1-DIG2-Recoil
                       200, 201, 202,  -1,  -1,  -1,  -1, 207,  -1,  -1,   //VME1-DIG3-Elum, 11
                        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,   //VME1-DIG4
                   
                         1,   0,   5,   4,   3,   2,   1,   0,  -1,  -1,/*1*/
                         3,   2,   1,   0,   5,   4,   3,   2,  -1,  -1,/*2*/
                        22,  10,   9,   8,   7,   6,   5,   4,  -1,  -1,/*3*/
                        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,//
                        
                         7,   6,  16,  10,   9,   8,   7,   6,  -1,  -1,/*4*/
                        15,  14,  13,  12,  16,  10,   9,   8,  -1,  -1,/*5*/
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
                         1,   1,   1,   1,   1,   1,   0,   0,  -1,  -1,//6
                        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
                        
                         0,   0,   2,   2,   2,   2,   2,   2,  -1,  -1,//7
                         1,   1,   1,   1,   0,   0,   0,   0,  -1,  -1,//8
                         2,   2,   2,   2,   2,   2,   1,   1,  -1,  -1,//9
                        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1};
