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
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author xujun on 2022-08-19-3:50 PM
 * @project SVCaller
 */
public class InsertionNew {
    int[] flags = null;
    String inputFile = null;
    String outputFile = null;
    String libPath = null;
    int chrNum = 0;
    int range = 0;
    double ratio = 0;
    int minl = 0;
    int quality = 0;
    ArrayList spliteSeq = null;
    HashSet hs = null;


    public InsertionNew(String[] parameters) {
        long startTime = System.currentTimeMillis();
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
        minl = Integer.valueOf(parameters[6]);
        quality = Integer.valueOf(parameters[7]);
        this.running();
        long endTime   = System.currentTimeMillis();
        NumberFormat formatter = new DecimalFormat("#0.00000");
        System.out.println(formatter.format((endTime - startTime) / 1000d));
    }
    public void running(){
        HashMap<Integer,Integer>[] linePosHM = new HashMap [chrNum];
        ArrayList[] infor = new ArrayList[chrNum];ArrayList[] inforS = new ArrayList[chrNum];
        List [] posLibL = new ArrayList[chrNum];
        for (int i =0 ; i < chrNum ; i++){
            linePosHM[i] = new HashMap();
            posLibL[i] = new ArrayList();
            infor[i] = new ArrayList();inforS[i] = new ArrayList();
        }
        try{
            BufferedReader brL = IOUtils.getTextGzipReader(new File(libPath).getAbsolutePath());
            String temp = null; String chr = null; int num = 0;
            List<String> tempT = new ArrayList<>(); StringBuilder sb = new StringBuilder();
            int row = 0;
            while ((temp = brL.readLine()) != null){
                tempT = PStringUtils.fastSplit(temp);
                chr = tempT.get(0);
                num = Integer.valueOf(chr.replaceAll("[a-zA-z]",""))-1;
                posLibL[num].add(tempT.get(1));
                infor[num].add(temp);
                linePosHM[num].put(row,Integer.parseInt(tempT.get(1)));
                row++;
            }
            brL.close();
            int posLibI [][] = new int[chrNum][];
            for (int i = 0 ; i < chrNum; i++ ){
                posLibI[i] = new int[posLibL[i].size()];
                linePosHM[i] = Methods.sortByValue(linePosHM[i]);
                for (Integer url : linePosHM[i].keySet()){
                    inforS[i].add(infor[i].get(url));
                }
                int rowNum = 0;
                for (Object url : linePosHM[i].values()){
                    posLibI[i][rowNum] = Integer.parseInt(url.toString());
                    rowNum++;
                }
                posLibL[i] = Arrays.stream(posLibI[i]).boxed().collect(Collectors.toList());
            }

            String cigar = null; String cigarSA = null; String SA = null;
            String [] value = null; String [] type = null;
            String [] value1 = null; String [] type1 = null;
            String chrSA = null;
            String strand = null; int length = 0; int q = 0;
            int start = 0; int end = 0;int start1 = 0; int end1 = 0;
            int pos = 0; int pos1 = 0; int num1 = 0;
            int flag =0; int tempFlag = 0;
            String seq = "";  String subSeq = null; String seqSA = null;
            SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(new File(inputFile).getAbsolutePath()));
            SAMRecordIterator r = sr.iterator();

            SimpleGapPenalty gapPen = new SimpleGapPenalty((short) 5, (short) 2);
            SubstitutionMatrix<NucleotideCompound> subMatrix = SubstitutionMatrixHelper.getNuc4_4();
            SequencePair<DNASequence, NucleotideCompound> psa = null;
            DNASequence query = null; DNASequence hit = null; DNASequence hitS = null;
            int nearstPos = 0; int index = 0; double a =0; double deno = 0;
            int breakPointQ = 0; int result[][][] =  new int[chrNum][][];
            spliteSeq = new ArrayList();
            for (int i =0;i < chrNum;i ++){
                result[i] = new int[posLibL[i].size()][3];
            }
            while(r.hasNext()){
                SAMRecord tem=r.next();
                if (tem.getFlags() > 2047)continue;
                cigar = tem.getCigarString();
                if(cigar.equals("*")) continue;
                if (cigar.contains("H"))continue;
                q = tem.getMappingQuality();
                if (q<quality)continue;
                boolean blS = false; boolean blSA = false;
                if (cigar.contains("S") ) blS = true;
                if (tem.getAttributes().get(0).tag.equals("SA")) blSA = true;
                if (!blS && !blSA) continue;
                spliteSeq.clear();
                chr = tem.getContig();
                num = Integer.valueOf(chr.replaceAll("[a-zA-Z]",""))-1;
                cigar = tem.getCigarString();
                value = cigar.replaceAll("[a-zA-Z=]"," ").split(" ");
                type = cigar.replaceAll("[0-9]+"," ").split(" ");
                start = tem.getStart(); end = tem.getEnd(); pos = value.length;
                seq = tem.getReadString();
                length = seq.length();
                if (blS && !blSA){
                    if (type[1].equals("S")){
                        sb.setLength(0);
                        int tl = Integer.valueOf(value[0]);
                        if (tl >= minl){
                            subSeq = seq.substring(0,tl);
                            sb.append(num).append("\t").append(start).append("\t").append(subSeq);
                            spliteSeq.add(sb.toString());
                        }
                    }
                    if (type[pos].equals("S")){
                        sb.setLength(0);
                        int tl = Integer.valueOf(value[pos-1]);
                        if (tl >= minl){
                            subSeq = seq.substring(length-tl,length);
                            sb.append(num).append("\t").append(end).append("\t").append(subSeq);
                            spliteSeq.add(sb.toString());
                        }
                    }
                }
                if (blSA){
                    length = 0;
                    for (int i = 0; i < value.length; i++){
                        if(!hs.contains(type[i+1])){
                            length+=Integer.parseInt(value[i]);
                        }
                    }
                    if (length!=seq.length())continue;
                    length = seq.length();
                    if (type[1].equals("S")){
                        sb.setLength(0);
                        int tl = Integer.valueOf(value[0]);
                        if (tl >= minl){
                            subSeq = seq.substring(0,tl);
                            sb.append(num).append("\t").append(start).append("\t").append(subSeq);
                            spliteSeq.add(sb.toString());
                        }
                    }
                    if (type[pos].equals("S")){
                        sb.setLength(0);
                        int tl = Integer.valueOf(value[pos-1]);
                        if (tl >= minl){
                            subSeq = seq.substring(length-tl,length);
                            sb.append(num).append("\t").append(end).append("\t").append(subSeq);
                            spliteSeq.add(sb.toString());
                        }
                    }
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
                        strand = "+";
                    }else {
                        strand = "-";
                    }
                    SA = tem.getAttributes().get(0).value.toString();
                    String[] SAArray = SA.split(";");
                    for (int i = 0; i < SAArray.length; i++){
                        temp = SAArray[i].replaceAll("SA:Z:","");
                        q = Integer.valueOf(temp.split(",")[4]);
                        if(q<quality)continue;
                        chrSA = temp.split(",")[0];
                        num1 = Integer.valueOf(chrSA.replaceAll("[a-zA-Z]",""))-1;
                        if (!strand.equals(temp.split(",")[2])){
                            seqSA = Methods.reverseComplement(seq);
                        }else{
                            seqSA = seq;
                        }
                        cigarSA = temp.split(",")[3];
                        start1 = Integer.valueOf(temp.split(",")[1]);
                        value1 = cigarSA.replaceAll("[a-zA-Z=]"," ").split(" ");
                        type1 = cigarSA.replaceAll("[0-9]+"," ").split(" ");
                        pos1 = value1.length; int s = 0; int sum =0;
                        if (type1[1].equals("S")){
                            sb.setLength(0); s = 1;
                            int tl = Integer.valueOf(value1[0]);
                            if (tl >= minl) {
                                subSeq = seqSA.substring(0,tl);
                                sb.append(num1).append("\t").append(start1).append("\t").append(subSeq);
                                spliteSeq.add(sb.toString());
                            }
                        }
                        if (type1[pos1].equals("S")){
                            for (int j = s; j < value1.length-1; j++){
                                if(!hs.contains(type1[j+1])){
                                    sum+=Integer.parseInt(value1[j]);
                                }
                            }
                            end1 = start1 + sum;
                            sb.setLength(0);
                            int tl = Integer.valueOf(value1[pos1-1]);
                            if (tl >= minl){
                                subSeq = seqSA.substring(length-tl,length);
                                sb.append(num1).append("\t").append(end1).append("\t").append(subSeq);
                                spliteSeq.add(sb.toString());
                            }
                        }
                    }
                }
                for (int i = 0; i < spliteSeq.size(); i++){
                    tempT = PStringUtils.fastSplit(spliteSeq.get(i).toString());
                    num = Integer.parseInt(tempT.get(0));
                    breakPointQ = Integer.valueOf(tempT.get(1));
                    if(posLibL[num].isEmpty())continue;
                    if(breakPointQ+range < posLibI[num][0])continue;
                    if (breakPointQ-range > posLibI[num][posLibI[num].length-1])continue;
                    int tempPos = breakPointQ-range;
                    if ( tempPos < 0){
                        tempPos = 0;
                    }
                    nearstPos = Methods.search(tempPos,posLibI[num]);
                    query = new DNASequence(tempT.get(2));
                    index = posLibL[num].indexOf(nearstPos);
                    if (posLibI[num][index]-breakPointQ < -range) index += 1;
                    for (int m = index ; m < posLibI[num].length ; m++) {
                        if (posLibI[num][m] - breakPointQ > range) break;
                        hit = new DNASequence(inforS[num].get(m).toString().split("\t")[2]);
                        if( hit.getLength() > 150 ){
                            hitS = new DNASequence(hit.getSequenceAsString().substring(0,150));
                            psa = Alignments.getPairwiseAlignment(query, hitS, Alignments.PairwiseSequenceAlignerType.LOCAL, gapPen, subMatrix);
                            if (hitS.getLength() < query.getLength()) {
                                deno = (double)hitS.getLength();
                            }else{
                                deno = (double)query.getLength();
                            }
                            a = (double)psa.getNumIdenticals()/deno;
                            if ( a >= this.ratio){
                                result[num][m][0] += 1; result[num][m][2] += 1;
                                continue;
                            }
                            hitS = new DNASequence(Methods.reverseComplement(hitS.toString()));
                            psa = Alignments.getPairwiseAlignment(query, hitS, Alignments.PairwiseSequenceAlignerType.LOCAL, gapPen, subMatrix);
                            if (hitS.getLength() < query.getLength()) {
                                deno = (double)hitS.getLength();
                            }else{
                                deno = (double)query.getLength();
                            }
                            a = (double)psa.getNumIdenticals()/deno;
                            if ( a >= this.ratio){
                                result[num][m][0] += 1; result[num][m][2] += 1;
                                continue;
                            }
                            hitS = new DNASequence(hit.getSequenceAsString().substring(hit.getLength()-150,hit.getLength()));
                            psa = Alignments.getPairwiseAlignment(query, hit, Alignments.PairwiseSequenceAlignerType.LOCAL, gapPen, subMatrix);
                            if (hitS.getLength() < query.getLength()) {
                                deno = (double)hitS.getLength();
                            }else{
                                deno = (double)query.getLength();
                            }
                            a = (double)psa.getNumIdenticals()/deno;
                            if ( a >= this.ratio){
                                result[num][m][1] += 1; result[num][m][2] += 1;
                                continue;
                            }
                            hitS = new DNASequence(Methods.reverseComplement(hitS.toString()));
                            psa = Alignments.getPairwiseAlignment(query, hitS, Alignments.PairwiseSequenceAlignerType.LOCAL, gapPen, subMatrix);
                            if (hitS.getLength() < query.getLength()) {
                                deno = (double)hitS.getLength();
                            }else{
                                deno = (double)query.getLength();
                            }
                            a = (double)psa.getNumIdenticals()/deno;
                            if ( a >= this.ratio){
                                result[num][m][1] += 1; result[num][m][2] += 1;
                                continue;
                            }
                            result[num][m][2] += 1;
                        }else{
                            hitS = hit;
//                            if (!hitS.getSequenceAsString().equals("ATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATA"))continue;
                            psa = Alignments.getPairwiseAlignment(query, hitS, Alignments.PairwiseSequenceAlignerType.LOCAL, gapPen, subMatrix);
                            if (hitS.getLength() < query.getLength()) {
                                deno = (double)hitS.getLength();
                            }else{
                                deno = (double)query.getLength();
                            }
                            a = (double)psa.getNumIdenticals()/deno;
                            if ( a >= this.ratio){
                                result[num][m][0] += 1; result[num][m][2] += 1;
                                continue;
                            }
                            hitS = new DNASequence(Methods.reverseComplement(hitS.toString()));
                            psa = Alignments.getPairwiseAlignment(query, hitS, Alignments.PairwiseSequenceAlignerType.LOCAL, gapPen, subMatrix);
                            if (hitS.getLength() < query.getLength()) {
                                deno = (double)hitS.getLength();
                            }else{
                                deno = (double)query.getLength();
                            }
                            a = (double)psa.getNumIdenticals()/deno;
                            if ( a >= this.ratio){
                                result[num][m][0] += 1; result[num][m][2] += 1;
                                continue;
                            }
                            result[num][m][2] += 1;
                        }

                    }
                }
            }
            sr.close();
            BufferedWriter bw = IOUtils.getTextGzipWriter(new File(this.outputFile).getAbsolutePath());
            for (int i = 0 ; i < chrNum ; i++) {
                for (int j = 0; j < result[i].length; j++) {
                    sb.setLength(0);
                    sb.append(inforS[i].get(j));
                    sb.append("\t").append(result[i][j][0]).append("\t").append(result[i][j][1]).append("\t").append(result[i][j][2]);
                    bw.write(sb.toString());bw.newLine();
                }
            }
            bw.flush();bw.close();
        }
        catch (Exception ex){
            ex.getStackTrace();
        }
    }
}
