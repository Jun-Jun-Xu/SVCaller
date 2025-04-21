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
 * @author xujun on 2022-07-27-8:10 AM
 * @project SVCaller
 */
public class UnmapInsertion {
    String inputFile = null;
    String outputFile = null;
    String libPath = null;
    double ratio = 0;
    ArrayList unmappedSeq = null;
    HashSet readName = null; HashMap nameSeqMap = null;

    public UnmapInsertion (String[] parameters) {
        inputFile = parameters[0];
        outputFile = parameters[1];
        libPath = parameters[2];
        ratio = Double.parseDouble(parameters[3]);
        long startTime = System.currentTimeMillis();
        unmappedSeq = this.getBreakPointAndSeq();
        this.alignment(unmappedSeq);
        long endTime   = System.currentTimeMillis();
        NumberFormat formatter = new DecimalFormat("#0.00000");
        System.out.println(formatter.format((endTime - startTime) / 1000d));
    }
    private ArrayList<String> getBreakPointAndSeq () {
        unmappedSeq = new ArrayList();
        readName = new HashSet(); nameSeqMap = new HashMap();
        try{
            String chr = null; String cigar = null; String name = null;
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
                if (tem.getReadUnmappedFlag() && tem.getMateUnmappedFlag()){
                    name = tem.getReadName();
                    seq = tem.getReadString();
                    if (!readName.contains(name)){
                        readName.add(name); nameSeqMap.put(name,seq);
                    }else{
                        sb.setLength(0); sb.append(nameSeqMap.get(name).toString()).append("\t").append(seq);
                        nameSeqMap.put(name,sb.toString());
                        readName.remove(name);
                    }
                }
            }
            sr.close();
            for (Object url : nameSeqMap.keySet()){
                unmappedSeq.add(nameSeqMap.get(url).toString());
            }
//            bw.flush();bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return unmappedSeq;
    }
    public void alignment(ArrayList infor){
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
        int nearstPos = 0; int index = 0;
        double [] a = new double[2]; boolean ab = true;
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
            while ((temp = brL.readLine()) != null) {
                tem = PStringUtils.fastSplit(temp);
                hit = new DNASequence(tem.get(2));
                int result = 0;
                for (int i = 0; i < infor.size(); i++) {
                    temQ = PStringUtils.fastSplit(unmappedSeq.get(i).toString());
                    for (int j = 0; j < 2; j++) {
                        query = new DNASequence(temQ.get(j));
                        psa = Alignments.getPairwiseAlignment(hit, query, Alignments.PairwiseSequenceAlignerType.LOCAL, gapPen, subMatrix);
                        a[j] = (double) psa.getNumIdenticals() / (double) psa.getLength();
                        if (a[j] < this.ratio) {
                            ab = false; break;
                        }
                    }
                    if (ab) {
                        result += 1;
                    }
                }
                bw.write(temp + "\t" + result);
                bw.newLine();
            }
            brL.close();
            bw.flush();bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
