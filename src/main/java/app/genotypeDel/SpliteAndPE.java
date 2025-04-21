package app.genotypeDel;

/**
 * @author xujun on 2022-06-08-3:30 PM
 * @project SVCaller
 */

import htsjdk.samtools.*;
import pgl.infra.utils.IOUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.util.*;
import utils.Methods;

public class SpliteAndPE {
    int[] flags = null;HashSet hs = null;
    HashSet readName = null; HashMap namePosMap = null;
    String inputFile = null;
    String outputDirS = null;
    String quality = null ;

    public SpliteAndPE(String[] parameters) {
        flags = new int[12];
        hs = new HashSet();
        hs.add((int)Math.pow(2,1));hs.add((int)Math.pow(2,2));hs.add((int)Math.pow(2,3));hs.add((int)Math.pow(2,8));hs.add((int)Math.pow(2,10));
        for (int i = 0 ; i < 12 ; i++){
            flags[i]= (int)Math.pow(2,i);
        }
        inputFile = parameters[0];
        outputDirS = parameters[1];
        quality = parameters[2];
        this.findBreakPoint();
    }

    private void findBreakPoint () {
        try{
            readName = new HashSet(); namePosMap = new HashMap();
            BufferedWriter bwS = IOUtils.getTextGzipWriter((new File(outputDirS+"_splite.txt.gz").getAbsolutePath()));
            BufferedWriter bwPE = IOUtils.getTextGzipWriter((new File(outputDirS+"_PE.txt.gz").getAbsolutePath()));
            bwS.write("chr1\tstart1\tend1\tchr2\tstart2\tend2");bwS.newLine();
            bwPE.write("chr1\tstart1\tend1\tchr2\tstart2\tend2");bwPE.newLine();
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
                boolean blS = false;boolean blPE = false;
                if (tem.getFlags() < 2047 || tem.getAttributes().get(0).value.toString().split(";").length!=1)blS=true;
                //|| !tem.getAttributes().get(0).tag.equals("SA")
                while (flag>0){
                    if (hs.contains(tempFlag)){
                        blPE=true;break;
                    }
                    flag = flag - tempFlag;
                    tempFlag = Methods.binarySearch(flag,flags);
                }
                if (blS && blPE)continue;
                name = tem.getReadName();chr = tem.getContig();start=tem.getStart();
                StringBuilder sb = new StringBuilder();
                if (!blS){
                    int sum = 0 ; int s = 0; int e = 0;
                    SA = tem.getAttributes().get(0).value.toString();
                    if (Integer.valueOf(SA.split(",")[4]) < Integer.parseInt(this.quality))continue;
                    String chrSA = SA.split(",")[0].replaceAll("SA:Z:","");
                    if (chr.equals(chrSA)){
                        cigar =  tem.getCigarString();
                        value = cigar.replaceAll("[a-zA-Z=]"," ").split(" ");
                        type = cigar.replaceAll("[0-9]+"," ").split(" ");
                        ArrayList typeList = new ArrayList<String>(Arrays.asList(type));
                        sum = 0;s=0; e = value.length;
                        if (type[1].equals("S") || type[1].equals("H")){s=1;};
                        if (type[e].equals("S") || type[e].equals("H")){e=e-1;};
                        for (int i=s ;i< e;i++){
                            if(!typeList.get(i+1).equals("I") && !typeList.get(i+1).equals("P")){
                                sum+=Integer.parseInt(value[i]);
                            }
                        }
                        cigarSA = SA.split(",")[3];
                        chr1 = SA.split(",")[0]; start1 = Integer.parseInt(SA.split(",")[1]);
                        value1 = cigarSA.replaceAll("[a-zA-Z=]"," ").split(" ");
                        type1 = cigarSA.replaceAll("[0-9]+"," ").split(" ");
                        ArrayList typeList1 = new ArrayList<String>(Arrays.asList(type1));
                        sum = 0;s=0; e = value1.length;
                        if (type1[1].equals("S") || type1[1].equals("H")){s=1;};
                        if (type1[e].equals("S") || type1[e].equals("H")){e=e-1;};
                        for (int i=s ;i< e;i++){
                            if(!typeList1.get(i+1).equals("I") && !typeList1.get(i+1).equals("P")){
                                sum+=Integer.parseInt(value1[i]);
                            }
                        }
                        if (Integer.valueOf(start1) > Integer.valueOf(start)){
                            sb.append(chr).append("\t").append(start).append("\t").append(start+sum).append("\t");
                            sb.append(chr1).append("\t").append(start1).append("\t").append(start1+sum);
                            if (output.add(sb.toString())){
                                bwS.write(sb.toString());bwS.newLine();
                            }
                        }else{
                            sb.append(chr1).append("\t").append(start1).append("\t").append(start1+sum).append("\t");
                            sb.append(chr).append("\t").append(start).append("\t").append(start+sum);
                            if (output.add(sb.toString())){
                                bwS.write(sb.toString());bwS.newLine();
                            }
                        }
                    }
                    continue;
                }
                if (!blPE){
                    chr = tem.getContig(); start = tem.getStart(); end = tem.getEnd();
                    sb.setLength(0);sb.append(chr).append("\t").append(start).append("\t").append(end+1).append("\t");
                    if (!readName.contains(name)){
                        readName.add(name); namePosMap.put(name,sb);
                    }else{
                        readName.remove(name);
                        bwPE.write(namePosMap.get(name).toString()+sb);bwPE.newLine();
                    }
                }

            }
            sr.close();
            bwS.flush();bwS.close();
            bwPE.flush();bwPE.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
