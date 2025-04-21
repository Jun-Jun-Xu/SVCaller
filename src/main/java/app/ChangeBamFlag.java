package app;

import htsjdk.samtools.*;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import utils.Methods;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * @author xujun on 2022-12-05-4:00 PM
 * @project SVCaller
 */
public class ChangeBamFlag {
    String inputFile = null;
    String outputFile = null;
    String inforFile = null;
    HashMap hm = null;
    int[] flags = null;

    public ChangeBamFlag(String[] parameters) {
        inputFile = parameters[0];
        outputFile = parameters[1];
        this.running();
    }
    public void running (){
        String temp = null;List<String> tem = new ArrayList<>();
        int flag = 0; int tempFlag = 0;
        try{
            SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(new File(inputFile).getAbsolutePath()));
            SAMRecordIterator r = sr.iterator();
            SAMFileHeader header=sr.getFileHeader();
            SAMFileWriterFactory samWriterFactory = new SAMFileWriterFactory();
            SAMFileWriter samWriter = samWriterFactory.makeBAMWriter(header, false, new File(outputFile));
            while(r.hasNext()){
                SAMRecord samtem=r.next();
                flag = samtem.getFlags();
                tempFlag = Methods.binarySearch(flag,flags);
                boolean bl = true;
                while (flag>0){
                    if (tempFlag==256){
                        bl=false;break;
                    }
                    flag = flag - tempFlag;
                    tempFlag = Methods.binarySearch(flag,flags);
                }
                if (!bl){
                    samtem.setFlags(samtem.getFlags()+1792);
                }
                samWriter.addAlignment(samtem);
            }
            samWriter.close();
        }
        catch (Exception ex){
            ex.getStackTrace();
        }
    }
}
