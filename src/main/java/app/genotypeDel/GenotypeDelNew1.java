package app.genotypeDel;
import htsjdk.samtools.*;
import pgl.infra.range.Range;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import utils.Methods;

import java.io.*;
import java.util.*;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
/**
 * @author xujun on 2023-08-18-4:45 PM
 * @project SVCaller
 */
public class GenotypeDelNew1 {
    /**
     * Chromosome num of this bam file.
     */
    int chrNum = 0;
    /**
     * Path of delition library
     */
    String libPath = null;
    /**
     * Path of samtools
     */
    String samPath = null;
    /**
     * Depth of this bam file.
     */
    int support = 0;
    String inputFile = null;
    String outPutFile = null;
    ArrayList[] breakPoint = null;
    ArrayList[] breakPointR = null;
    ArrayList[] breakPointS = null;
    int[] flags = null;
    HashSet hs = null;
    String prefix = null;

    public GenotypeDelNew1 (String[] parameterS) {
        flags = new int[12];
        for (int i = 0 ; i < 12 ; i++){
            flags[i]= (int)Math.pow(2,i);
        }
        hs = new HashSet();
        hs.add("D");hs.add("N");hs.add("H");hs.add("P");
        inputFile = parameterS[0];
        outPutFile = parameterS[1];
        libPath = parameterS[2];
        samPath = parameterS[3];
        chrNum = Integer.parseInt(parameterS[4]);
        support = Integer.valueOf(parameterS[5]);
        breakPoint = this.getBreakPoint();
        this.breakRange(breakPoint);
    }
    private ArrayList<String>[] getBreakPoint () {
        int frq = 0;
        breakPoint = new ArrayList [chrNum];breakPointR = new ArrayList [chrNum];
        for(int i = 0 ; i < chrNum ; i++){
            breakPoint[i] = new ArrayList();
            breakPointR[i] = new ArrayList();
        }
        String cigar = null; String cigarSA = null; String SA = null;
        String [] value = null; String [] type = null;
        String [] value1 = null; String [] type1 = null;
        int quality = 0; int quality1 = 0;
        String chr = null; String chrSA = null;
        int start = 0; int end = 0;int start1 = 0; int end1 = 0;
        int pos = 0; int pos1 = 0; int num = 0; int num1 = 0;
        String seq = ""; boolean bl = true;
        String reads = null; int line = 0;
        try{
            SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(new File(inputFile).getAbsolutePath()));
            SAMRecordIterator r = sr.iterator();
            while(r.hasNext()){
                SAMRecord tem=r.next();
                line++;
                reads = tem.getReadName();
//                if(!reads.equals("A00358:51:H5JNCDMXX:2:1478:12057:10347")) continue;
                if (tem.getFlags() > 2047)continue;
                if (!tem.getAttributes().get(0).tag.equals("SA"))continue;
                seq = tem.getReadString();
                chr = tem.getContig();
//                bw.write(chr+"\t"+tem.getStart()+"\t"+seq);bw.newLine();
                bl = Pattern.matches("[a-zA-Z]",chr);
                if (bl){
                    num = Integer.valueOf(chr.replaceAll("[a-zA-Z]",""))-1;
                    prefix = chr.replaceAll("[0-9]","");
                }else{
                    num = Integer.parseInt(chr)-1;
                }
                if (num<0){continue;}
                cigar = tem.getCigarString();
                value = cigar.replaceAll("[a-zA-Z=]"," ").split(" ");
                type = cigar.replaceAll("[0-9]+"," ").split(" ");
                int length = 0;
                for (int i = 0; i < value.length; i++){
                    if(!hs.contains(type[i+1])){
                        length+=Integer.parseInt(value[i]);
                    }
                }
                if (length!=seq.length())continue;
                quality = tem.getMappingQuality();
                if (quality >= 20){
                    start = tem.getStart(); end = tem.getEnd(); pos = value.length;
                    if (type[1].equals("H") || type[1].equals("S")){
                        breakPoint[num].add(start);
                    }
                    if (type[pos].equals("H") || type[pos].equals("S")){
                        breakPoint[num].add(end);
                    }
                }
                SA = tem.getAttributes().get(0).value.toString();
                String[] SAArray = SA.split(";");
                for (int i = 0; i < SAArray.length; i++){
                    String temp = SAArray[i].replaceAll("SA:Z:","");
                    quality1 = Integer.parseInt(temp.split(",")[4]);
                    if (quality1 >= 20){
                        chrSA = temp.split(",")[0];
                        if ( Integer.parseInt(chrSA) == 0) continue;
                        num1 = Integer.valueOf(chrSA.replaceAll("[a-zA-Z]",""))-1;
                        cigarSA = temp.split(",")[3];
                        start1 = Integer.valueOf(temp.split(",")[1]);
                        value1 = cigarSA.replaceAll("[a-zA-Z=]"," ").split(" ");
                        type1 = cigarSA.replaceAll("[0-9]+"," ").split(" ");
                        pos1 = value1.length; int s = 0; int sum =0;
                        if (type1[1].equals("S") || type1[1].equals("H")){
                            s = 1;
                            breakPoint[num1].add(start1);
                        }
                        if (type1[pos1].equals("S") || type1[pos1].equals("H")){
                            for (int j = s; j < value1.length-1; j++){
                                if(!hs.contains(type1[j+1])){
                                    sum+=Integer.parseInt(value1[j]);
                                }
                            }
                            end1 = start1 + sum - 1;
                            breakPoint[num1].add(end1);
                        }
                    }
                }
            }
            sr.close();
            for (int i = 0 ;i < breakPoint.length ; i++){
                ArrayList tem = Methods.removeDuplicates(breakPoint[i]);
                for (int j = 0 ; j < tem.size() ; j++){
                    int temp = Integer.parseInt(tem.get(j).toString());
                    frq = Collections.frequency(breakPoint[i],temp);
                    if (frq < support)continue;
                    breakPointR[i].add(temp);
                }

            }

//            bw.flush();bw.close();
        }
        catch (Exception e) {
            System.out.println(reads+"\t"+chr+"\t"+num+"\t"+start+"\t"+end);
            System.out.println(line);
            e.printStackTrace();
        }
        return breakPointR;
    }
    public void breakRange (ArrayList[] infor){
        HashMap[] posInfor = new HashMap [chrNum];
        HashMap<Integer,Integer>[] posInforS = new HashMap [chrNum];
        List[] inforS = new ArrayList [chrNum];
        HashMap<Integer,String>[] hs = new HashMap[chrNum];
        int posLibI [][] = new int[chrNum][]; List [] posLibL = new ArrayList[chrNum];
        breakPointS = new ArrayList[chrNum];
        for (int i = 0 ; i < chrNum; i++ ){
            breakPointS[i] = new ArrayList();
            posInfor[i] = new HashMap<String, Integer>();
            hs[i] = new HashMap<Integer,String>();
//            if (infor[i].size() ==0 ) continue;
            posLibI[i] = new int[infor[i].size()];
            posLibL[i] = new ArrayList();
            for(int j = 0 ; j < infor[i].size(); j++ ){
                int pos  = Integer.valueOf(infor[i].get(j).toString());
                posInfor[i].put( Integer.valueOf(j), pos);
            }
        }

        for (int i = 0 ; i < chrNum; i++ ){
            posInforS[i] = new HashMap(); inforS[i] = new ArrayList();
//            if (posInfor[i].isEmpty()) continue;
            posInforS[i] = Methods.sortByValue(posInfor[i]);
            for (Integer url : posInforS[i].keySet()){
                inforS[i].add(infor[i].get(url));
            }
            int rowNum = 0;
            for (Integer url : posInforS[i].values()){
                posLibI[i][rowNum] = url;
                rowNum++;
            }
            posLibL[i] = Arrays.stream(posLibI[i]).boxed().collect(Collectors.toList());
        }
        BufferedReader brL = IOUtils.getTextGzipReader((new File(libPath).getAbsolutePath()));
        String temp = null;  List<String> tem = new ArrayList<>();
        StringBuilder sb = new StringBuilder();
        Range rg = null; int chr = 0; int start = 0; int end =0;
        try{
            while ((temp = brL.readLine()) != null){
                tem = PStringUtils.fastSplit(temp);
                boolean bl = Pattern.matches("[a-zA-Z]",tem.get(0));
                if (bl){
                    chr = Integer.valueOf(tem.get(0).replaceAll("[a-zA-Z]",""))-1;
                }else{
                    chr = Integer.parseInt(tem.get(0))-1;
                }
                start = Integer.parseInt(tem.get(1)); end = Integer.parseInt(tem.get(2))-1;
                rg = new Range(chr,start,end);
                if ( posLibL[chr].isEmpty()){
                    sb.setLength(0); sb.append(start).append("\t").append(end+1).append("\t0");
                    breakPointS[chr].add(sb.toString());
                    continue;
                }
                int nearst = Methods.search(rg.getRangeStart(),posLibI[chr]);
                int index = posLibL[chr].indexOf(nearst);
                for (int i = index ; i < posLibI[chr].length; i++){
                    if (posLibI[chr][i] > rg.getRangeEnd()){
                        sb.setLength(0); sb.append(rg.getRangeStart()).append("\t").append(end+1).append("\t1");
                        breakPointS[chr].add(sb.toString());
                        break;
                    }
                    int point = posLibI[chr][i];
                    if (rg.isContain(chr,point)){
                        sb.setLength(0); sb.append(rg.getRangeStart()).append("\t").append(point).append("\t1");
                        breakPointS[chr].add(sb.toString());
                        if ( (point+1) == end )break;
                        rg.setRangeStart(point+1);
                    }
                }
            }
            brL.close();
            BufferedWriter bw = IOUtils.getTextGzipWriter(outPutFile);
            for (int i = 0 ;i < breakPointS.length ; i++){
                breakPointS[i] = Methods.removeDuplicates(breakPointS[i]);
                for (int j = 0; j < breakPointS[i].size(); j++){
                    sb.setLength(0);
                    sb.append(i+1).append("\t").append(breakPointS[i].get(j));
                    bw.write(sb.toString()); bw.newLine();
                }
            }
            bw.flush();bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
