package app.genotypeDel;

/**
 * @author xujun on 2022-07-13-10:11 AM
 * @project SVCaller
 */
import htsjdk.samtools.*;
import pgl.infra.utils.IOUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.util.*;
import utils.Methods;

public class FindURPE {
    int[] flags = null;HashSet hs = null;//if the flag of reads is contained in hs, the read is not unreasonable PE read.
    HashSet readName = null; HashMap namePosMap = null;
    String inputFile = null;
    String outputDirS = null;
    String quality = null ;

    public FindURPE(String[] parameters){
        flags = new int[12];
        hs = new HashSet();
        hs.add((int)Math.pow(2,1));hs.add((int)Math.pow(2,2));hs.add((int)Math.pow(2,3));hs.add((int)Math.pow(2,8));hs.add((int)Math.pow(2,10));
        for (int i = 0 ; i < 12 ; i++){
            flags[i]= (int)Math.pow(2,i);
        }
        inputFile = parameters[0];
        outputDirS = parameters[1];
        quality = parameters[2];
        this.find();
    }

    private void find () {
        try{
            readName = new HashSet(); namePosMap = new HashMap();
            BufferedWriter bwPE = IOUtils.getTextGzipWriter((new File(outputDirS+"_PE.txt.gz").getAbsolutePath()));
            SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(new File(inputFile).getAbsolutePath()));
            SAMRecordIterator r = sr.iterator();
            int flag =0; int tempFlag = 0;
            String cigar = null; String cigarSA = null; String SA = null;
            String [] value = null; String [] type = null;
            String [] value1 = null; String [] type1 = null;
            String name = null ; String chr = null; int start = 0; int end = 0;
            String chr1 = null; int start1 =0;HashSet output = new HashSet();
            while(r.hasNext()) {
                SAMRecord tem=r.next();
                if (tem.getMappingQuality() < Integer.parseInt(this.quality))continue;
                if (!tem.getReferenceName().equals(tem.getMateReferenceName()))continue;
                flag = tem.getFlags();
                tempFlag = Methods.binarySearch(flag,flags);
                boolean blPE = false;
                //|| !tem.getAttributes().get(0).tag.equals("SA")
                while (flag>0){
                    if (hs.contains(tempFlag)){
                        blPE=true;break;
                    }
                    flag = flag - tempFlag;
                    tempFlag = Methods.binarySearch(flag,flags);
                }
                if (!blPE){
                    bwPE.write(tem.getSAMString());//bwPE.newLine();
                }
            }
            sr.close();
            bwPE.flush();bwPE.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
