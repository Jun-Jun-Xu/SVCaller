package app.Insertion;

import htsjdk.samtools.*;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import utils.Methods;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author xujun on 2022-06-28-8:41 PM
 * @project SVCaller
 */

public class Insertion1 {
    String inputFile = null;
    String outputFile = null;
    String libPath = null;
    int chrNum = 0;
    String quality = null ;
    int range = 0;
    double ratio = 0;
    HashMap<String, SAMRecord>[] nameInforMap = null;
    HashSet[] spliteReadsName = null;
    ArrayList[] spliteSeq = null;

    public Insertion1(String[] parameters) {
        inputFile = parameters[0];
        outputFile = parameters[1];
        quality = parameters[2];
        chrNum = Integer.valueOf(parameters[3]);
        libPath = parameters[4];
        range =  Integer.valueOf(parameters[5]);
        ratio = Double.parseDouble(parameters[6]);
        spliteSeq = this.getBreakPointAndSeq();
        this.alignment(spliteSeq);
    }
    private ArrayList<String>[] getBreakPointAndSeq () {
        HashMap<String,SAMRecord>[] nameInforMap = new HashMap[chrNum];
        HashSet[] spliteReadsName = new HashSet[chrNum];
        ArrayList[] spliteSeq = new ArrayList [chrNum];
        for(int i = 0 ; i < chrNum ; i++){
            nameInforMap[i] = new HashMap<String,SAMRecord>();
            spliteReadsName[i] = new HashSet();
            spliteSeq[i] = new ArrayList();
        }
        try{
            String cigar = null; String cigarSA = null; String SA = null;
            String [] value = null; String [] type = null;
            String [] value1 = null; String [] type1 = null;
            String name = null ; String chr = null;
            int start = 0; int end = 0;int start1 = 0; int end1 = 0;
            int pos = 0; int pos1 = 0; int num = 0;
            StringBuilder sb = new StringBuilder();
            String seq = ""; String seqR = "";
            SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(new File(inputFile).getAbsolutePath()));
            SAMRecordIterator r = sr.iterator();
            while(r.hasNext()){
                SAMRecord tem=r.next();
                if (tem.getMappingQuality() < Integer.parseInt(this.quality))continue;
                if (!tem.getReferenceName().equals(tem.getMateReferenceName()))continue;
                if (tem.getFlags() < 2047 )continue;
                chr = tem.getContig(); name = tem.getReadName(); SA = tem.getAttributes().get(0).value.toString();
                num = Integer.parseInt(chr.replaceAll("[a-zA-Z]",""))-1;
                if (SA.split(";").length==1){
                    if (Integer.valueOf(SA.split(",")[4]) < Integer.parseInt(this.quality))continue;
                    String chrSA = SA.split(",")[0].replaceAll("SA:Z:","");
                    if (!chr.equals(chrSA))continue;
                    spliteReadsName[num].add(name);nameInforMap[num].put(name,tem);
                }
            }
            SamReader sr1 = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(new File(inputFile).getAbsolutePath()));
            SAMRecordIterator r1 = sr1.iterator();
            while(r1.hasNext()) {
                SAMRecord tem=r1.next();
                if (tem.getMappingQuality() < Integer.parseInt(this.quality))continue;
                if (!tem.getReferenceName().equals(tem.getMateReferenceName()))continue;
                if (tem.getFlags() > 2047 )continue;
                chr = tem.getContig(); name = tem.getReadName();
                num = Integer.parseInt(chr.replaceAll("[a-zA-Z]",""))-1;
                if (!spliteReadsName[num].contains(name))continue;
                SAMRecord temS = nameInforMap[num].get(name);
                SA = temS.getAttributes().get(0).value.toString();
                cigarSA = SA.split(",")[3];
                int tempTem = tem.getInferredInsertSize();
                int tempTemS = temS.getInferredInsertSize();
                if (tempTem*tempTemS<0)continue;
                cigar = temS.getCigarString();
                byte[] seqB = tem.getReadBases();
                value = cigar.replaceAll("[a-zA-Z=]"," ").split(" ");
                type = cigar.replaceAll("[0-9]+"," ").split(" ");
                start = temS.getStart(); end = temS.getEnd(); pos = value.length;
                value1 = cigarSA.replaceAll("[a-zA-Z=]"," ").split(" ");
                type1 = cigarSA.replaceAll("[0-9]+"," ").split(" ");
                start1 = tem.getStart(); end1 = tem.getEnd(); pos1 = value1.length;
                String typeTem = cigar.replaceAll("[0-9]","").replaceAll("[SH]","*");
                String typeTemS = cigarSA.replaceAll("[0-9]","").replaceAll("[SH]","*");
                boolean bl = false;
                if (typeTem.startsWith("*") && typeTemS.startsWith("*") && !typeTem.endsWith("*") && !typeTemS.endsWith("*"))bl=true;
                if (!typeTem.startsWith("*") && !typeTemS.startsWith("*") && typeTem.endsWith("*") && typeTemS.endsWith("*"))bl=true;
                if (!bl){
                    sb.setLength(0);
                    if (type[1].equals("S") || type[1].equals("H")){
                        seq = "";
                        for (int i = 0 ; i < Integer.valueOf(value[0]); i++){
                            seq += (char)seqB[i];
                        }
                        sb.append(start).append("\t").append(seq).append("\t").append("1").append("\n");
                    }
                    if (type[pos].equals("S") || type[pos].equals("H")){
                        seq = "";
                        int a = 150-Integer.valueOf(value[pos-1])-2;
                        for (int i = a ; i < 150; i++){
                            seq += (char)seqB[i];
                        }
                        sb.append(end).append("\t").append(seq).append("\t").append("0").append("\n");
                    }
                    if (type1[1].equals("S") || type1[1].equals("H")){
                        seq = "";
                        for (int i = 0 ; i < Integer.valueOf(value1[0]); i++){
                            seq += (char)seqB[i];
                        }
                        sb.append(start1).append("\t").append(seq).append("\t").append("1").append("\n");
                    }
                    if (type1[pos1].equals("S") || type1[pos1].equals("H")){
                        seq = "";
                        int a = 150-Integer.valueOf(value1[pos1-1])-2;
                        for (int i = a ; i < 150; i++){
                            seq += (char)seqB[i];
                        }
                        sb.append(end1).append("\t").append(seq).append("\t").append("0").append("\n");
                    }
                    if (sb.length()==0)continue;
                    String[] spliteR = sb.toString().split("\n");
                    for (int i =0 ; i < spliteR.length; i++){
                        spliteSeq[num].add(spliteR[i]);
                    }
                    //                    bw.write(sb.toString());
                }else{
                    seq = "";
                    for (int i =0;i < seqB.length ; i++){
                        seq += (char)seqB[i];
                    }
                    seqR = Methods.reverseComplement(seq);
                    sb.setLength(0);
                    if (type[1].equals("S") || type[1].equals("H")){
                        seq = "";
                        seq = seqR.substring(0,Integer.valueOf(value[0]));
                        sb.append(start).append("\t").append(seq).append("\t").append("1").append("\n");
                    }
                    if (type[pos].equals("S") || type[pos].equals("H")){
                        seq = "";
                        int a = 150-Integer.valueOf(value[pos-1])-1;
                        seq = seqR.substring(a,150);
                        sb.append(end).append("\t").append(seq).append("\t").append("0").append("\n");
                    }
                    if (type1[1].equals("S") || type1[1].equals("H")){
                        seq = "";
                        for (int i = 0 ; i < Integer.valueOf(value1[0]); i++){
                            seq += (char)seqB[i];
                        }
                        sb.append(start1).append("\t").append(seq).append("\t").append("1").append("\n");
                    }
                    if (type1[pos1].equals("S") || type1[pos1].equals("H")){
                        seq = "";
                        int a = 150-Integer.valueOf(value1[pos1-1])-1;
                        for (int i = a ; i < 150; i++){
                            seq += (char)seqB[i];
                        }
                        sb.append(end1).append("\t").append("\t").append(seq).append("\t").append("0").append("\n");
                    }
                    if (sb.length()==0)continue;
                    //                    bw.write(sb.toString());
                    String[] spliteR = sb.toString().split("\n"); HashMap posInfor = new HashMap();
                    List al = new ArrayList(); ArrayList a = new ArrayList();
                    for (int i =0 ; i < spliteR.length; i++){
                        al = PStringUtils.fastSplit(spliteR[i].toString());
                        posInfor.put(al.get(0), spliteR[i]);
                        a.add(al.get(0));
                    }
                    Collections.sort(a);
                    for (int i = 0 ; i < a.size(); i++){
                        spliteSeq[num].add(posInfor.get(a.get(i)));
                    }
                    seqR = "";
                }
            }
            sr.close(); sr1.close();
            //            bw.flush();bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return spliteSeq;
    }
    public void alignment(ArrayList[] infor){
        HashMap[] posInfor = new HashMap [chrNum];
        HashMap<Integer,Integer>[] posInforS = new HashMap [chrNum];
        List[] inforS = new ArrayList [chrNum];
        List<String> reads = new ArrayList();
        HashMap<Integer,String>[] hs = new HashMap[chrNum];
        int posLibI [][] = new int[chrNum][]; List [] posLibL = new ArrayList[chrNum];
        for (int i = 0 ; i < chrNum; i++ ){
            posInfor[i] = new HashMap<String, Integer>();
            hs[i] = new HashMap<Integer,String>();
            posLibI[i] = new int[infor[i].size()];
            posLibL[i] = new ArrayList();
            for(int j = 0 ; j < infor[i].size(); j++ ){
                reads = PStringUtils.fastSplit(infor[i].get(j).toString());
                posInfor[i].put( Integer.valueOf(j), Integer.valueOf(reads.get(0)));
            }
        }
        for (int i = 0 ; i < chrNum; i++ ){
            posInforS[i] = new HashMap(); inforS[i] = new ArrayList();
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
        String chr = null; int breakPointH = 0; int breakPointQ = 0; int num = -1 ;
        SimpleGapPenalty gapPen = new SimpleGapPenalty((short) 5, (short) 2);
        SubstitutionMatrix<NucleotideCompound> subMatrix = SubstitutionMatrixHelper.getNuc4_4();
        SequencePair<DNASequence, NucleotideCompound> psa = null;
        BufferedWriter bw = IOUtils.getTextGzipWriter(new File(this.outputFile).getAbsolutePath());
        DNASequence query = null; DNASequence hit = null; DNASequence hitS = null;
        List al = new ArrayList();
        BufferedReader brL = IOUtils.getTextGzipReader((new File(libPath).getAbsolutePath()));
        String temp = null; List<String> tem = new ArrayList<>(); List<String> temQ = new ArrayList<>();
        int errPos = -1; int position = 0;
        int nearstPos = 0; int index = 0; List indexS = new ArrayList();
        try{
            while ((temp = brL.readLine()) != null){
                tem = PStringUtils.fastSplit(temp);
                chr = tem.get(0);
                breakPointH = Integer.valueOf(tem.get(1));
                num = Integer.valueOf(chr.replaceAll("[a-zA-z]",""))-1;
                int tempPos = breakPointH-range;
                if ( tempPos < 0){
                    tempPos = 0;
                }
                nearstPos = Methods.search(tempPos,posLibI[num]);
                if(Math.abs(nearstPos - breakPointH) > range){
                    bw.write(temp+"\t0\t0");bw.newLine();continue;
                }
                hit = new DNASequence(tem.get(2));
                index = posLibL[num].indexOf(nearstPos);
                int [] result = new int[2];
                for (int m = index ; m < posLibI[num].length ; m++){
                    if (Math.abs(posLibI[num][m]-breakPointH) > range)break;
                    temQ = PStringUtils.fastSplit(inforS[num].get(m).toString());
                    query = new DNASequence(temQ.get(1));
                    position = Integer.valueOf(temQ.get(2));
                    if ( hit.getLength() > 150){
                        if (position==0){
                            hitS = new DNASequence(hit.getSequenceAsString().substring(0,150));
                        }else{
                            int length = hit.getLength();
                            hitS = new DNASequence(hit.getSequenceAsString().substring(length-150,length));
                        }
                    }else {
                        hitS = hit ;
                    }
                    psa = Alignments.getPairwiseAlignment(query, hitS, Alignments.PairwiseSequenceAlignerType.LOCAL, gapPen, subMatrix);
                    double a = (double)psa.getNumIdenticals()/(double) query.getLength();
                    if( a < this.ratio)continue;
                    errPos = m;
                    result[position] += 1;
                }
                bw.write(temp+"\t"+result[0]+"\t"+result[1]);bw.newLine();
                //                System.out.println(temp+"\t"+result[0]+"\t"+result[1]);
            }
            brL.close();
            bw.flush();bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.out.println(hit);
            System.out.println(query);
            System.out.println(errPos);
        }
    }
}