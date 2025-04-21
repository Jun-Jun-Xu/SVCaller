package app.Insertion;

import htsjdk.samtools.*;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import utils.Methods;
import java.io.BufferedWriter;
import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author xujun on 2022-09-01-11:22 PM
 * @project SVCaller
 */
public class MakeInsLib {
    HashSet refHS = null;
    HashSet queHS = null;
    String inputFile = null;
    String outputFile = null;
    int chrNum = 0;
    int quality = 0;
    int length = 0;
    int range = 0;
    int[] flags = null;
    HashSet[] nameHS = null; HashMap<String, Integer>[] nameFreqHM = null;
    HashMap<String, ArrayList>[] nameInfor = null;
    ArrayList[] result = null; ArrayList[] resultS = null;

    private int chr;
    private int start;

    class chrStart {
        private int chr;
        private int start;
        public chrStart(int chr, int start) {
            this.chr = chr;
            this.start = start;
        }
        public int getChr() {
            return chr;
        }
        public int getStart() {
            return start;
        }
    }

    public MakeInsLib (String[] parameters) {
        long startTime = System.currentTimeMillis();
        refHS = new HashSet();
        refHS.add("M"); refHS.add("D"); refHS.add("N"); refHS.add("="); refHS.add("X");
        queHS = new HashSet();
        queHS.add("M"); queHS.add("I"); queHS.add("S"); queHS.add("="); queHS.add("X");
        flags = new int[12];
        for (int i = 0 ; i < 12 ; i++){
            flags[i]= (int)Math.pow(2,i);
        }
        inputFile = parameters[0];
        outputFile = parameters[1];
        chrNum = Integer.parseInt(parameters[2]);
        quality = Integer.parseInt(parameters[3]);
        length = Integer.parseInt(parameters[4]);
        range = Integer.valueOf(parameters[5]);
        this.running();
        long endTime   = System.currentTimeMillis();
        NumberFormat formatter = new DecimalFormat("#0.00000");
        System.out.println(formatter.format((endTime - startTime) / 1000d));
    }

    public void running (){
        result = new ArrayList[chrNum]; resultS = new ArrayList[chrNum];
        for (int i = 0 ; i< chrNum ; i++){
            result[i] = new ArrayList<String>(); resultS[i] = new ArrayList<String>();
        }
        String name = null; String chr = null; int num = 0;
        ArrayList al = new ArrayList(); ArrayList alt = new ArrayList();
        String cigar = null ; String seq = null; int q = 0; String SA = null;
        int start = 0; int flag = 0; int tempFlag = 0; int strandMain = 0 ;int strand = 0;
        StringBuilder sb = new StringBuilder();
        try{
            SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(new File(inputFile).getAbsolutePath()));
            SAMRecordIterator r = sr.iterator();
            while(r.hasNext()){
                SAMRecord tem=r.next();
                cigar = tem.getCigarString();
                if(cigar.equals("*")) continue;
                if (tem.getFlags() > 2047) continue;
                name = tem.getReadName();
                chr = tem.getContig();
                num = Integer.valueOf(chr.replaceAll("[a-zA-z]",""))-1;
                if (!tem.getAttributes().get(0).tag.equals("SA")) {
                    sb.setLength(0);
                    start = tem.getStart(); q = tem.getMappingQuality(); seq = tem.getReadString();
                    sb.append(name+"\t").append(start+"\t").append(cigar+"\t").append(seq);
                    if (q < this.quality) continue;
                    alt = this.findInCigarNew(sb.toString());
                    if (alt.isEmpty())continue;
                    for (int i =0 ; i< alt.size(); i++){
                        result[num].add(alt.get(i));
                    }
                }else{
                    flag = tem.getFlags();
                    tempFlag = Methods.binarySearch(flag,flags);
                    boolean bl = true;
                    while (flag>0){
                        if (tempFlag==16){
                            bl=false;break;
                        }
                        flag = flag - tempFlag;
                        tempFlag = Methods.binarySearch(flag,flags);
                    }
                    if (bl){
                        strandMain = 0; //0 is positive.
                    }else {
                        strandMain = 1; //1 is negative.
                    }
                    sb.setLength(0); al.clear();
                    start = tem.getStart(); q = tem.getMappingQuality(); seq = tem.getReadString();
                    sb.append(num+"\t").append(start+"\t").append(q+"\t").append(cigar+"\t").append(strandMain);
                    al.add(sb.toString());
                    SA = tem.getAttributes().get(0).value.toString();
                    SA.replaceAll("SA:Z:","");
                    String[] SAArray = SA.split(";");
                    for (int i = 0; i < SAArray.length; i++){
                        String [] SATem = SAArray[i].split(",");
                        int numT = Integer.valueOf(SATem[0].replaceAll("[a-zA-z]",""))-1;
                        if (SATem[2].equals("+")) {
                            strand = 0;
                        }else{
                            strand = 1;
                        }
                        sb.setLength(0);
                        sb.append(numT+"\t").append(SATem[1]+"\t").append(SATem[4]+"\t").append(SATem[3]+"\t").append(strand);
                        al.add(sb.toString());
                    }
                    alt = this.findInSupplNew(name, strandMain, quality, seq, al);
                    if (alt.isEmpty())continue;
                    for (int i =0 ; i< alt.size(); i++){
                        num = Integer.parseInt(alt.get(i).toString().split("\t")[0]);
                        result[num].add(alt.get(i).toString().substring(1).replaceFirst("\t","").replaceFirst("[0-9]m","m"));
                    }
                }
            }
            r.close();

            BufferedWriter bw = IOUtils.getTextGzipWriter((new File(outputFile).getAbsolutePath()));
            for (int i = 0 ; i < chrNum ; i++){
                int row = 0; HashMap<Integer, Integer> sortHM = new HashMap();
                for (int j = 0 ; j < result[i].size() ; j++){
                    sortHM.put(row, Integer.parseInt(result[i].get(j).toString().split("\t")[1].replaceAll("\\[","")));
                    row++;
                }
                sortHM = Methods.sortByValue(sortHM);
                for (Integer url : sortHM.keySet()){
                    resultS[i].add(result[i].get(url));
                }
                for (int j = 0 ; j < resultS[i].size() ; j++){
                    bw.write(i+1+"\t");bw.write(resultS[i].get(j).toString().replaceAll("\\[","").replaceAll("\\]",""));
                    bw.newLine();
                }
            }
            bw.flush(); bw.close();
        }
        catch (Exception ex){
            System.out.println(name);
            ex.getStackTrace();
        }
    }

    public ArrayList findInCigarNew (String infor){ //name start cigar seq
        String cigar = null; String [] value = null; String [] type = null;
        int rSum = 0; int qSum = 0; int start = 0; StringBuilder sb = new StringBuilder();
        ArrayList al = new ArrayList(); ArrayList result = new ArrayList();
        List<String> tem = new ArrayList<>(); tem = PStringUtils.fastSplit(infor);
        cigar = tem.get(2);
        value = cigar.replaceAll("[a-zA-Z=]"," ").split(" ");
        type = cigar.replaceAll("[0-9]+"," ").split(" ");
        ArrayList typeList = new ArrayList<String>(Arrays.asList(type));
        ArrayList<Integer> allIndexes = new ArrayList(); ArrayList indexs = new ArrayList();
        allIndexes =
                (ArrayList<Integer>) IntStream.range(1, typeList.size()).boxed()
                        .filter(i -> typeList.get(i).equals("I"))
                        .collect(Collectors.toList());
        for (int i = 0; i < allIndexes.size() ; i++){
            int length = Integer.parseInt(value[allIndexes.get(i).intValue()-1]);
            if (length >= this.length){
                indexs.add(allIndexes.get(i));
            }
        }
        if (!indexs.isEmpty()){
            start = Integer.parseInt(infor.split("\t")[1])-1;
            for (int i =0 ; i < indexs.size() ; i++){
                qSum = 0 ; rSum = 0;
                for (int j = 0; j < Integer.parseInt(indexs.get(i).toString())-1 ; j++){
                    if (refHS.contains(type[j+1])) rSum += Integer.valueOf(value[j]);
                    if (queHS.contains(type[j+1])) qSum += Integer.valueOf(value[j]);
                }
                sb.setLength(0); int pos = start+rSum;
                sb.append(tem.get(0)+"\t").append(pos+"\t").append(tem.get(3).substring(qSum,qSum+Integer.valueOf(value[Integer.parseInt(indexs.get(i).toString())-1]))+"\tIns_Cigar");
                result.add(sb.toString());
            }
        }

        allIndexes.clear(); indexs.clear();
        allIndexes =
                (ArrayList<Integer>) IntStream.range(1, typeList.size()).boxed()
                        .filter(i -> typeList.get(i).equals("D"))
                        .collect(Collectors.toList());
        for (int i = 0; i < allIndexes.size() ; i++){
            int length = Integer.parseInt(value[allIndexes.get(i).intValue()-1]);
            if (length >= this.length){
                indexs.add(allIndexes.get(i));
            }
        }
        if (!indexs.isEmpty()){
            start = Integer.parseInt(infor.split("\t")[1])-1;
            for (int i =0 ; i < indexs.size() ; i++){
                rSum = 0; int tmp =0;
                for (int j = 0; j < Integer.parseInt(indexs.get(i).toString())-1 ; j++){
                    if (refHS.contains(type[j+1])) rSum += Integer.valueOf(value[j]);
                    tmp = j;
                }
                sb.setLength(0); int pos = start+rSum; int pos1 = pos+ Integer.parseInt(value[tmp+1]);
                sb.append(tem.get(0)+"\t").append(pos+"\t").append(pos1+"\t").append("Del_Cigar");
                result.add(sb.toString());
            }
        }
        return result;
    }

    public ArrayList findInSupplNew (String name, int strand, Integer quality, String seq, ArrayList infor){
        ArrayList result = new ArrayList(); List<String> tem = new ArrayList<>();
        HashMap<Integer, Integer> hs = new HashMap();
        ArrayList subInfor = new ArrayList(); ArrayList subInforS = new ArrayList();
        StringBuilder sb = new StringBuilder(); ArrayList al = new ArrayList();
        for (int i = 0 ; i < infor.size() ; i++){
            tem = PStringUtils.fastSplit(infor.get(i).toString());
            if (Integer.parseInt(tem.get(2)) < quality) continue;
            if (Integer.parseInt(tem.get(4))==strand){
                subInfor.add(infor.get(i));
            }
        }
        if (subInfor.size()==0) return result;
        String cigar[] = new String[subInfor.size()];
        int start[] = new int[subInfor.size()];
        int chr[] = new int[subInfor.size()];
        String value[][] = new String[subInfor.size()][];
        String type[][] = new String[subInfor.size()][];
        List<chrStart> chrStarts = new ArrayList();
        for (int i = 0 ; i < subInfor.size();i++){
            tem = PStringUtils.fastSplit(subInfor.get(i).toString());
            chr[i] = Integer.parseInt(tem.get(0));
            start[i] = Integer.parseInt(tem.get(1));
            chrStarts.add(new chrStart(chr[i], start[i]));
        }
        Comparator<chrStart> compareByVariable = Comparator
                .comparing(chrStart::getChr)
                .thenComparing(chrStart::getStart);
        List<chrStart> sort = chrStarts.stream()
                .sorted(compareByVariable)
                .collect(Collectors.toList());
        for (int i = 0; i < sort.size(); i++){
            int index = chrStarts.indexOf(sort.get(i));
            subInforS.add(subInfor.get(index));
        }
        tem = PStringUtils.fastSplit(subInforS.get(0).toString());
        chr[0] = Integer.parseInt(tem.get(0));
        ArrayList group[] = new ArrayList[subInforS.size()];
        ArrayList groupFinal[] = new ArrayList[subInforS.size()];
        for (int i = 0 ; i < subInforS.size(); i++){
            group[i] = new ArrayList();
            groupFinal[i] = new ArrayList();
        }
        group[0].add(subInforS.get(0)); int num = 0;
        for (int i = 1; i < subInforS.size(); i++){
            tem = PStringUtils.fastSplit(subInforS.get(i).toString());
            chr[i] = Integer.parseInt(tem.get(0));
            if(chr[i-1]!=chr[i]){
                num++;
                group[num].add(subInforS.get(i));
            }else{
                group[num].add(subInforS.get(i));
            }
        }
        num = 0;
        for (int i = 0 ; i < group.length; i++){
            if (group[i].isEmpty()) continue;
            if (group[i].size()==1){
                groupFinal[num].add(group[i].get(0));
                num++;continue;
            }
            int FS[] = new int[group[i].size()];
            tem = PStringUtils.fastSplit(group[i].get(0).toString());
            chr[0] = Integer.parseInt(tem.get(0));
            start[0] = Integer.parseInt(tem.get(1));
            cigar[0] = tem.get(3);
            value[0] = cigar[0].replaceAll("[a-zA-Z=]"," ").split(" ");
            FS[0] = Integer.parseInt(value[0][0]);
            groupFinal[num].add(group[i].get(0));
            for (int j = 1; j < group[i].size(); j++){
                tem = PStringUtils.fastSplit(group[i].get(j).toString());
                chr[j] = Integer.parseInt(tem.get(0));
                start[j] = Integer.parseInt(tem.get(1));
                if (chr[j] == chr[j-1] && start[j] == start[j-1])break;
                cigar[j] = tem.get(3);
                value[j] = cigar[j].replaceAll("[a-zA-Z=]"," ").split(" ");
                FS[j] = Integer.parseInt(value[j][0]);
                if (FS[j] - FS[j-1] < 0) continue;
                if (start[j]-start[j-1] > 150000){
                    num++;
                    groupFinal[num].add(group[i].get(j));
                }else{
                    groupFinal[num].add(group[i].get(j));
                }
            }
        }

        int refBreakpoint = 0; int refSumT = 0; int length = 0; String subSeq = null; int startPos = 0;
        int queSumT = 0 ; int queSum = 0;
        al.clear();
        for (int i = 0; i< groupFinal.length; i++){
            for(int j = 0 ; j < groupFinal[i].size() ; j++){
                tem = PStringUtils.fastSplit(groupFinal[i].get(j).toString());
                chr[i] = Integer.parseInt(tem.get(0));
                start[i] = Integer.parseInt(tem.get(1));
                cigar[i] = tem.get(3);
                sb.setLength(0);
                sb.append(name+"\t").append(start[i]+"\t").append(cigar[i]+"\t").append(seq);
                al = this.findInCigarNew(sb.toString());
                if (al.isEmpty())continue;
                for (int k =0 ; k< al.size(); k++){
                    result.add(chr[i]+"\t"+al.get(k));
                }
            }
        }
        for (int i = 0; i < groupFinal.length; i++){
            if (groupFinal[i].size()>1){
                tem = PStringUtils.fastSplit(groupFinal[i].get(0).toString());
                chr[0] = Integer.parseInt(tem.get(0));
                start[0] = Integer.parseInt(tem.get(1));
                cigar[0] = tem.get(3);
                value[0] = cigar[0].replaceAll("[a-zA-Z=]"," ").split(" ");
                type[0] = cigar[0].replaceAll("[0-9]+"," ").split(" ");
                for(int j = 1; j < groupFinal[i].size(); j++){
                    tem = PStringUtils.fastSplit(groupFinal[i].get(j).toString());
                    chr[j] = Integer.parseInt(tem.get(0));
                    start[j] = Integer.parseInt(tem.get(1));
                    cigar[j] = tem.get(3);
                    value[j] = cigar[j].replaceAll("[a-zA-Z=]"," ").split(" ");
                    type[j] = cigar[j].replaceAll("[0-9]+"," ").split(" ");
                    if (chr[j-1] != chr[j]) continue;
                    refSumT = 0; queSumT = 0;
                    for (int k = 0 ; k < value[j-1].length-1; k++){
                        if (refHS.contains(type[j-1][k+1])){
                            refSumT += Integer.parseInt(value[j-1][k]);
                        }
                        if (queHS.contains(type[j-1][k+1])){
                            queSumT += Integer.parseInt(value[j-1][k]);
                        }
                    }
                    queSum = queSumT + Integer.parseInt(value[j-1][value[j-1].length-1]);
                    refBreakpoint = start[j-1] + refSumT -1;
                    int difR = start[j]-refBreakpoint;
                    int difQ = queSumT- Integer.parseInt(value[j][0]);
                    if ( Math.abs(difR) < range){
                        if ( difQ > this.length){
                            length = Integer.parseInt(value[j-1][value[j-1].length-1])+refBreakpoint-start[j]+1-queSum+Integer.parseInt(value[j][0]);
                            if ( length < 0) continue;
                            startPos = Integer.parseInt(value[j][0])+refBreakpoint-start[j];
                            subSeq = seq.substring(startPos-length, startPos);
                            sb.setLength(0);
                            sb.append(chr[j]+"\t"+name+"\t"+refBreakpoint).append("\t").append(subSeq+"\tIns_Suppl");
                            result.add(sb.toString());
                        }
                    }else{
                        if (difR > this.length){
                            if(difQ > this.length) continue;
                            if (Math.abs(difQ) < range){
                                sb.setLength(0);
                                sb.append(name+"\t"+refBreakpoint+"\t"+(start[j]-1)+"\tDel_Suppl");
                                result.add(chr[j]+"\t"+sb.toString());
                            }
                        }
                    }
                }
            }
        }
        return result;
    }
}
