#!/usr/bin/python

from epics import caget
from epics import caput
import os
import time

print("======== DataBase for PSD Array ========")

loop = 1
tOld = 0

while loop == 1:
    f = open("tempDB.txt", 'w')
    t1 = int(round(time.time() * 1000 ))
    
    fexp=open("/home/helios/digios/expName.sh", 'r')
    line=fexp.readline()
    line=fexp.readline()
    expName=line[8:-1]
    line=fexp.readline()
    runNum=line[11:-1]
    fexp.close()
    
    string = "RunNum value=%s\n" % (runNum)
    f.write(string)

    bashCommand="du -c ~/digios/analysis/data/%s_run_%03d* | tail -1 | awk '{print $1}'" % (expName, int(runNum))
    result=os.popen(bashCommand).read().rstrip()
    print("ExpName : %s, runNum : %03d, file size : %s" % (expName, int(runNum),  result))
    string="fileSize value=%s\n" % (result)
    f.write(string)

    for VME in range (1,5):
        pv="DAQC%d_CV_BuffersAvail" % (VME)
        result=caget(pv)
        print("VME%d Buffer" % VME, result)
        if ( float(result) < 300 ) :
            caput("Online_CS_StartStop", "Stop")
            print("=================== VME%d Buffer < 300!!!!!" % VME)
            time.sleep(10)
            caput("Online_CS_StartStop", "Start")
            print("=================== resume")
            f.close()
            break

        string="buffer,VME=%d value=%s\n" % (VME, result)
        f.write(string)
        for DIG in range (1,5):
            for CH in range(0, 10):
                #timeStamp = int(round(time.time() * 1000 ))
                pv="VME0%d:MDIG%d:disc_count%d_RBV" % (VME, DIG, CH)
                result=caget(pv)
                string="hitRate,VME=%d,DIG=%d,CH=%d value=%s\n" % (VME, DIG, CH, result)
                f.write(string)
                
                pv="VME0%d:MDIG%d:led_threshold%d" % (VME, DIG, CH)
                result=caget(pv)
                string="threshold,VME=%d,DIG=%d,CH=%d value=%s\n" % (VME, DIG, CH, result)
                f.write(string)

    result=caget("VME32:MTRG:RAW_TRIG_RATE_COUNTER_2_RBV")
    string="sumHitX value=%s\n" % (result)
    f.write(string)                

    result=caget("VME32:MTRG:RAW_TRIG_RATE_COUNTER_3_RBV")
    string="sumHitY value=%s\n" % (result)
    f.write(string)

    f.close()
    #lineNum = len(open("tempDB.txt").readlines())
    #print "=== file Length : ", lineNum
    
    #fh = open('tempDB.txt','r')
    #for i in fh.readlines(): print(i[:-1])

    t2 = int(round(time.time() * 1000 ))
    tDiff = t2-t1

    print "=======================", t2-tOld
    #print t2, ", diff = ", tDiff , " mses"
    
    #usually take 4000 msec for all channels
    if( t2-tOld > 2000 ) :
        #bashCommand='curl -sS -i -XPOST "http://heliosdb.onenet:8086/write?db=testing" --data-binary @tempDB.txt --max-time 10'
        bashCommand='curl -sS -i -XPOST "http://heliosdb.onenet:8086/write?db=testing" --data-binary @tempDB.txt --speed-time 5 --speed-limit 1000'
        os.system(bashCommand)
        tOld = t2
    
    time.sleep(3)


