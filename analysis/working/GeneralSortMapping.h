/***********************************************************************
 * 
 * This is the Mapping for Armory/GeneralSort.C
 * 
 * idDetMap = [  0, 100) ; PSD
 *          = [100, 110) ; Recoil
 *          = (200, 240] ; ELUM
 *          = [300, 310) ; EZERO
 *          = [400, 450] ; TAC & RF-Timing
 * 
***********************************************************************/

//Arrays for mapping things...
//id = ((VME_ID -1 )x4 + (Dig_ID-1)) x 10 + CH_ID + 1010

//detID =   0 ...  99 = PSD array
//detID = 100 ... 110 = recoil
//detID = 200 ... 240 = ELUM
//detID = 300 ... 310 = EZERO
//detID = 400 ... 450 = TAC timing


int idConst = 1010;
Int_t idDetMap[160] = {  4,   3,   2,   1,   0,   4,   3,   2,   1,   0,   //VME1-DIG1-TL-E,R
                         4,   3,   2,   1,   0,   4,   3,   2,   1,   0,   //VME1-DIG2-TL-F,N
                        29,  28,  27,  26,  25,  29,  28,  27,  26,  25,   //VME1-DIG3-BL-E,R
                       100, 101, 102, 103, 104, 105, 106, 107,  -1,  -1,   //VME1-DIG4
                        29,  28,  27,  26,  25,  29,  28,  27,  26,  25,   //VME2-DIG1-BL-F,N
                        24,  23,  22,  21,  20,  24,  23,  22,  21,  20,   //VME2-DIG2-BB-E,R
                        24,  23,  22,  21,  20,  24,  23,  22,  21,  20,   //VME2-DIG3-BB-F,N
                       400,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,   //VME2-DIG4
                        19,  18,  17,  16,  15,  19,  18,  17,  16,  15,   //VME3-DIG1-BR-E,R
                        19,  18,  17,  16,  15,  19,  18,  17,  16,  15,   //VME3-DIG2-BR-F,N
                        14,  13,  12,  11,  10,  14,  13,  12,  11,  10,   //VME3-DIG3-TR-E,R
                       200, 201, 202, 203, 204, 205, 206, 207,  -1,  -1,   //VME3-DIG4
                        14,  13,  12,  11,  10,  14,  13,  12,  11,  10,   //VME4-DIG1-TR-F,N
                         9,   8,   7,   6,   5,   9,   8,   7,   6,   5,   //VME4-DIG2-TT-E,R
                         9,   8,   7,   6,   5,   9,   8,   7,   6,   5,   //VME4-DIG3-TT-F,N
                       208, 209, 210, 211, 212, 213, 214, 215,  -1,  -1};  //VME4-DIG4

//Kind Map, 0 = energy, 1 = Far, 2 = Near, 3 = Ring
//kind map for other detectors is not used.
Int_t idKindMap[160] = { 0,   0,   0,   0,   0,   3,   3,   3,   3,   3,//VME1
                         2,   2,   2,   2,   2,   1,   1,   1,   1,   1,
                         0,   0,   0,   0,   0,   3,   3,   3,   3,   3,
                        -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
                         2,   2,   2,   2,   2,   1,   1,   1,   1,   1,//VME2
                         0,   0,   0,   0,   0,   3,   3,   3,   3,   3,
                         2,   2,   2,   2,   2,   1,   1,   1,   1,   1,
                        -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
                         0,   0,   0,   0,   0,   3,   3,   3,   3,   3,//VME3
                         2,   2,   2,   2,   2,   1,   1,   1,   1,   1,
                         0,   0,   0,   0,   0,   3,   3,   3,   3,   3,
                        -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
                         2,   2,   2,   2,   2,   1,   1,   1,   1,   1,//VME4
                         0,   0,   0,   0,   0,   3,   3,   3,   3,   3,
                         2,   2,   2,   2,   2,   1,   1,   1,   1,   1,
                        -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2};
