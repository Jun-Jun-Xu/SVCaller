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
 * @author xujun on 2022-07-26-2:46 PM
 * @project SVCaller
 */
public class MateInsertion {
    int[] flags = null;
    String inputFile = null;
    String outputFile = null;
    String libPath = null;
    int chrNum = 0;
    int range = 0;
    double ratio = 0;
    ArrayList[] unmappedSeq = null;
    HashSet hs = null;

    public MateInsertion (String[] parameters) {
        flags = new int[12];
        for (int i = 0 ; i < 12 ; i++){
            flags[i]= (int)Math.pow(2,i);
        }
        hs = new HashSet();
        hs.add("D");hs.add("N");hs.add("H");hs.add("P");
        inputFile = parameters[0];
        outputFile = parameters[1];
        chrNum = Integer.valueOf(parameters[2]);
        libPath = parameters[3];
        range =  Integer.valueOf(parameters[4]);
        ratio = Double.parseDouble(parameters[5]);
        unmappedSeq = this.getBreakPointAndSeq();
        this.alignment(unmappedSeq);
    }
    private ArrayList<String>[] getBreakPointAndSeq () {
        unmappedSeq = new ArrayList [chrNum];
        for(int i = 0 ; i < chrNum ; i++){
            unmappedSeq[i] = new ArrayList();
        }
        try{
            String chr = null; String cigar = null;
            String [] value = null; String [] type = null;
            int start = 0; int end = 0;
            int num = 0;
            StringBuilder sb = new StringBuilder();
            String seq = "";
            SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(new File(inputFile).getAbsolutePath()));
            SAMRecordIterator r = sr.iterator();
//            BufferedWriter bw = IOUtils.getTextGzipWriter(new File("/Users/xujun/Desktop/Tempory/lib.txt.gz").getAbsolutePath());
            while(r.hasNext()){
                SAMRecord tem=r.next();
                if (!tem.getReadUnmappedFlag())continue;
                chr = tem.getMateReferenceName();
                num = Integer.valueOf(chr.replaceAll("[a-zA-Z]",""))-1;
                start = tem.getMateAlignmentStart();
                seq = tem.getReadString();
                cigar = tem.getAttributes().get(0).value.toString();
                value = cigar.replaceAll("[a-zA-Z=]"," ").split(" ");
                type = cigar.replaceAll("[0-9]+"," ").split(" ");
                sb.setLength(0); int sum = 0;
                int s = 0 ; int e = value.length;
                if (type[1].equals("H") || type[1].equals("S")) s = 1;
                if (type[e].equals("H") || type[e].equals("S")) e = value.length-1;
                for (int i = s ; i < e ; i++){
                    if (!hs.contains(type[i+1])) sum += Integer.valueOf(value[i]);
                }
                end = start + sum;
                sb.append(end+"\t").append(seq);
                unmappedSeq[num].add(sb.toString());
            }
            sr.close();
//            bw.flush();bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return unmappedSeq;
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
        String chr = null; int breakPointH = 0; int num = -1 ;
        SimpleGapPenalty gapPen = new SimpleGapPenalty((short) 5, (short) 2);
        SubstitutionMatrix<NucleotideCompound> subMatrix = SubstitutionMatrixHelper.getNuc4_4();
        SequencePair<DNASequence, NucleotideCompound> psa = null;
        BufferedWriter bw = IOUtils.getTextGzipWriter(new File(this.outputFile).getAbsolutePath());
        DNASequence query = null; DNASequence hit = null; DNASequence hitS = null;
        List al = new ArrayList();
        BufferedReader brL = IOUtils.getTextGzipReader((new File(libPath).getAbsolutePath()));
        String temp = null; List<String> tem = new ArrayList<>(); List<String> temQ = new ArrayList<>();
        int errPos = -1;
        int nearstPos = 0; int index = 0; double a =0;
//        try{
//            BufferedWriter bwL = IOUtils.getTextGzipWriter(new File("/Users/xujun/Desktop/Tempory/lib2.txt.gz").getAbsolutePath());
//            for (int i =0 ;i < chrNum ; i++){
//                for (int j =0 ;j < inforS[i].size();j++){
//                    bwL.write((i+1)+"\t"+inforS[i].get(j).toString());bwL.newLine();
//                }
//            }
//            bwL.flush();bwL.close();
//        }
//        catch (Exception e) {
//            e.printStackTrace();
//        }
        try{
            while ((temp = brL.readLine()) != null){
                tem = PStringUtils.fastSplit(temp);
                chr = tem.get(0);
                breakPointH = Integer.valueOf(tem.get(1));
//                if(breakPointH!=137349)continue;
                num = Integer.valueOf(chr.replaceAll("[a-zA-z]",""))-1;
                if(posLibL[num].isEmpty()){
                    bw.write(temp+"\t0\t0");bw.newLine();continue;
                }
                if(breakPointH+range < posLibI[num][0]) {
                    bw.write(temp+"\t0\t0");bw.newLine();continue;
                }
                if (breakPointH-range > posLibI[num][posLibI[num].length-1]){
                    bw.write(temp+"\t0\t0");bw.newLine();continue;
                }
                int tempPos = breakPointH-range;
                if ( tempPos < 0){
                    tempPos = 0;
                }
                nearstPos = Methods.search(tempPos,posLibI[num]);
                errPos = tempPos;
//                if(Math.abs(nearstPos - breakPointH) > range){
////                    System.out.println(temp+"\t0\t0");
//                    bw.write(temp+"\t0\t0");bw.newLine();continue;
//                }
                hit = new DNASequence(tem.get(2));
                index = posLibL[num].indexOf(nearstPos);
                int result = 0;
                for (int m = index ; m < posLibI[num].length ; m++){
                    if ((posLibI[num][m]-breakPointH) > range)break;
                    temQ = PStringUtils.fastSplit(inforS[num].get(m).toString());
                    query = new DNASequence(temQ.get(1));
                    hitS = hit;
                    psa = Alignments.getPairwiseAlignment(query, hitS, Alignments.PairwiseSequenceAlignerType.LOCAL, gapPen, subMatrix);
                    a = (double)psa.getNumIdenticals()/(double) psa.getLength();
                    if ( a >= this.ratio){
                        result += 1; continue;
                    }
                    hitS = new DNASequence(Methods.reverseComplement(hit.toString()));
                    psa = Alignments.getPairwiseAlignment(query, hitS, Alignments.PairwiseSequenceAlignerType.LOCAL, gapPen, subMatrix);
                    a = (double)psa.getNumIdenticals()/(double) psa.getLength();
                    if ( a >= this.ratio){
                        result += 1; continue;
                    }
                }
                bw.write(temp+"\t"+result);bw.newLine();
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
