package app.Insertion;

import htsjdk.samtools.*;
import org.checkerframework.checker.units.qual.A;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import utils.Methods;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

/**
 * @author xujun on 2022-08-30-5:45 PM
 * @project SVCaller
 */

public class MakeLib {

    HashSet refHS = null;
    HashSet queHS = null;
    String inputFile = null;
    String outputFile = null;
    String libPath = null;
    int chrNum = 0;
    HashSet[] nameHSC = null; HashMap<String, ArrayList>[] nameInforC = null;
    HashSet[] nameHSS = null; HashMap<String, ArrayList>[] nameInforS = null;
    ArrayList[] insInfor = null; ArrayList[] insInforS = null;

    public MakeLib (String[] parameters) {
        long startTime = System.currentTimeMillis();
        refHS = new HashSet();
        refHS.add("M"); refHS.add("D"); refHS.add("N"); refHS.add("="); refHS.add("X");
        queHS = new HashSet();
        queHS.add("M"); queHS.add("I"); queHS.add("S"); queHS.add("="); queHS.add("X");
        inputFile = parameters[0];
        outputFile = parameters[1];
        libPath = parameters[2];
        chrNum = Integer.parseInt(parameters[3]);
        this.running();
        long endTime   = System.currentTimeMillis();
        NumberFormat formatter = new DecimalFormat("#0.00000");
        System.out.println(formatter.format((endTime - startTime) / 1000d));
    }
    public void running(){
        HashSet[] nameHSC = new HashSet[chrNum]; HashSet[] nameHSS = new HashSet[chrNum];
        HashMap<String, ArrayList>[] nameInforC = new HashMap[chrNum];
        HashMap<String, ArrayList>[] nameInforS = new HashMap[chrNum];
        ArrayList[] insInfor = new ArrayList[chrNum]; ArrayList[] insInforS = new ArrayList[chrNum];
        for (int i = 0 ; i < chrNum ; i++){
            nameHSC[i] = new HashSet(); nameInforC[i] = new HashMap();
            nameHSS[i] = new HashSet(); nameInforS[i] = new HashMap();
            insInfor[i] = new ArrayList(); insInforS[i] = new ArrayList();
        }
        try{
            BufferedReader brL = IOUtils.getTextGzipReader(new File(libPath).getAbsolutePath());
            String temp = null; String infor = null;
            String chr = null; int num = 0;
            List<String> tempT = new ArrayList<>(); List<String> tempt = new ArrayList<>();
            StringBuilder sb = new StringBuilder();
            while ((temp = brL.readLine()) != null){
                tempT = PStringUtils.fastSplit(temp);
                infor = tempT.get(5);
                chr = tempT.get(0);
                num = Integer.valueOf(chr.replaceAll("[a-zA-z]",""))-1;
                tempt = PStringUtils.fastSplit(infor.replaceAll("\\]\\[","\t").replaceAll("\\[","").replaceAll("\\]",""));
                for (int i = 0 ; i < tempt.size() ; i++){
                    tempT = PStringUtils.fastSplit(tempt.get(i).replaceAll("\\|","\t").replaceAll("\\;","\t"));
                    sb.setLength(0);
                    int length = Integer.parseInt(tempT.get(2))-Integer.parseInt(tempT.get(1));
                    sb.append(tempT.get(1)).append("\t").append(length);
                    if (tempT.get(4).equals("cigar")){
                        if (nameHSC[num].contains(tempT.get(5))){
                            ArrayList al = nameInforC[num].get(tempT.get(5));
                            al.add(sb); nameInforC[num].put(tempT.get(5),al);
                        }else{
                            nameHSC[num].add(tempT.get(5));
                            nameInforC[num].put(tempT.get(5), new ArrayList(Arrays.asList(sb.toString())));
                        }
                    }else{
                        if (nameHSS[num].contains(tempT.get(5))){
                            ArrayList al = nameInforS[num].get(tempT.get(5));
                            al.add(sb); nameInforS[num].put(tempT.get(5),al);
                        }else{
                            nameHSS[num].add(tempT.get(5));
                            nameInforS[num].put(tempT.get(5), new ArrayList(Arrays.asList(sb.toString())));
                        }
                    }
                }
            }
            brL.close();

            SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(new File(inputFile).getAbsolutePath()));
            SAMRecordIterator r = sr.iterator();
            String name = null; String seq = null; String cigar = null;
            int pos1 = 0; int pos2 = 0; int flag = 0; int length = 0; int point = 0; int sum = 0 ;
            String [] value = null; String [] type = null;
            ArrayList al = new ArrayList();
            while(r.hasNext()){
                SAMRecord tem=r.next();
                name = tem.getReadName();
                chr = tem.getContig(); num = Integer.valueOf(chr.replaceAll("[a-zA-z]",""))-1;
                seq = tem.getReadString();
                flag = tem.getFlags();
                if (flag > 2047) continue;
                if (!nameHSC[num].contains(name))continue;
                al = nameInforC[num].get(name);
                for (int i = 0 ; i < al.size() ; i++){
                    sb.setLength(0);
                    infor = al.get(i).toString();
                    pos1 = Integer.parseInt(infor.split("\t")[0]); pos2 = tem.getStart();
                    length = Integer.parseInt(infor.split("\t")[1]);
                    cigar = tem.getCigarString();
                    value = cigar.replaceAll("[a-zA-Z=]"," ").split(" ");
                    type = cigar.replaceAll("[0-9]+"," ").split(" ");
                    for (int j = 0; j < value.length; j++){
                        if(refHS.contains(type[j+1])){
                            pos2+=Integer.parseInt(value[j]);
                            if (pos1==pos2-1) point = j;
                        }
                    }
                    sum = 0;
                    for (int j = 0 ; j < point+1 ; j++){
                        if (queHS.contains(type[j+1])){
                            sum += Integer.parseInt(value[j]);
                        }
                    }
                    sb.append(pos1).append("\t").append(seq.substring(sum,sum+length));
                    insInfor[num].add(sb.toString());
                }
                if(!nameHSS[num].contains(name))continue;
                al = nameInforS[num].get(name);
                for (int i = 0 ; i < al.size() ; i++){
                    sb.setLength(0);
                    infor = al.get(i).toString();
                    pos1 = Integer.parseInt(infor.split("\t")[0]); pos2 = tem.getStart();
                    if (pos2-pos1!=1)continue;
                    length = Integer.parseInt(infor.split("\t")[1]);
                    cigar = tem.getCigarString();
                    value = cigar.replaceAll("[a-zA-Z=]"," ").split(" ");
                    type = cigar.replaceAll("[0-9]+"," ").split(" ");
                    sb.append(pos1).append("\t").append(seq.substring(Integer.parseInt(value[1])-length,length));
                    insInfor[num].add(sb.toString());
                }
            }
            sr.close();
            HashMap<Integer,Integer>[] linePos = new HashMap[chrNum];
            int row = 0;
            BufferedWriter bw = IOUtils.getTextGzipWriter(new File(this.outputFile).getAbsolutePath());
            for (int i = 0 ; i < chrNum ; i++){
                linePos[i] = new HashMap();
                row = 0;
                for (int j = 0 ; j < linePos[i].size() ; j++){
                    linePos[i].put(row,Integer.parseInt(insInfor[i].get(j).toString().split("\t")[0]));
                    row++;
                }
                HashMap<Integer,Integer> hs = Methods.sortByValue(linePos[i]);
                for (Integer url : hs.keySet()){
                    insInforS[i].add(insInfor[i].get(url));
                }
                for (int j = 0; j < insInforS[i].size(); j++){
                    bw.write(num+1); bw.write("\t"); bw.write(insInforS[i].get(j).toString());
                    bw.newLine();
                }
            }

        }
        catch (Exception ex){
            ex.getStackTrace();
        }
    }
}
