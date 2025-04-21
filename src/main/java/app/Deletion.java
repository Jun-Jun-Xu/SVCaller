package app;

/**
 * @author xujun on 2022-06-06-10:23 PM
 * @project SVCaller
 */

import htsjdk.samtools.*;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class Deletion {
    String inputFile = null;
    String outputFile = null;
    /**
     * The minimum length of the SV
     */
    String length = null ;
    /**
     * The minimum quality of alignment
     */
    String quality = null ;
    /**
     * Number of threads
     */
    int threadNum = 16;
    String inputFileSDir = null;
    String outputFileSDir = null;
    SamReader sr = null;
    SamReader srS [] = null;
    String fileS[] = null;
    BufferedWriter bw =  null;
    BufferedWriter bwS [] = null;

    public Deletion(String[] parameters) {
        if (parameters[0].contains("bam")){
            inputFile = parameters[0];
            outputFile = parameters[1];
        }else{
            inputFileSDir = parameters[0];
            outputFileSDir = parameters[1];
        }
        length = parameters[2];
        quality = parameters[3];
        threadNum = Integer.parseInt(parameters[4]);
        this.callDeletion();
    }
    private void callDeletion () {
        try{
            if (!inputFile.equals(null)){
                sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(inputFile));
                bw = IOUtils.getTextWriter((new File(outputFile).getAbsolutePath()));
                bw.write("chr\tstart\tend\tname");bw.newLine();
                SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(inputFile));
                SAMRecordIterator r = sr.iterator();
                String cigar = null;String SA = null;
                String [] value = null; String [] type = null;
                String [] value1 = null; String [] type1 = null;
                String name = null ; String chr = null; int start = 0;
                while(r.hasNext()) {
                    SAMRecord tem=r.next();
                    cigar=tem.getCigarString();
                    if (tem.getMappingQuality() < Integer.parseInt(this.quality))continue;
                    if (!cigar.contains("D")) continue;
                    value = cigar.replaceAll("[a-zA-Z=]"," ").split(" ");
                    type = cigar.replaceAll("[0-9]+"," ").split(" ");
                    ArrayList typeList = new ArrayList<String>(Arrays.asList(type));
                    ArrayList<Integer> allIndexes =
                            (ArrayList<Integer>) IntStream.range(1, typeList.size()).boxed()
                                    .filter(i -> typeList.get(i).equals("D"))
                                    .collect(Collectors.toList());
                    boolean bl = false;
                    for (int i = 0; i < allIndexes.size() ; i++){
                        if (Integer.parseInt(value[allIndexes.get(i).intValue()-1])>=Integer.valueOf(this.length)){
                            bl=true;break;
                        }
                    }
                    if (bl==false)continue;
                    name = tem.getReadName();
                    chr = tem.getContig();
                    start=tem.getStart()-1;
                    for (int i =0 ; i < allIndexes.size() ; i++){
                        int num = Integer.parseInt(value[allIndexes.get(i)-1]);
                        if (num<Integer.parseInt(this.length))continue;
                        int sum=0;StringBuilder sb = new StringBuilder();
                        int s=0;if (type[1].equals("S")){s=1;};
                        for (int j = s; j < allIndexes.get(i)-1; j++){
                            if(!typeList.get(j+1).equals("I") && !typeList.get(j+1).equals("P")){
                                sum+=Integer.parseInt(value[j]);
                            }
                        }
                        int startPos=start+sum; int endPos=startPos+num;
                        sb.append(chr).append("\t").append(startPos).append("\t").append(endPos).append("\t").append(name).append("\tCIGAR");
                        bw.write(sb.toString());bw.newLine();
                    }
                    if (tem.getFlags()>2047){
                        SA = tem.getAttributes().get(0).value.toString();
                        String chrSA = SA.split(",")[0].replaceAll("SA:Z:","");
                        if (!chr.equals(chrSA))continue;
                        String cigarSA = SA.split(",")[3];
                        value1 = cigarSA.replaceAll("[a-zA-Z=]"," ").split(" ");
                        type1 = cigarSA.replaceAll("[0-9]+"," ").split(" ");
                        ArrayList typeList1 = new ArrayList<String>(Arrays.asList(type1));
                        int sum = 0;int s=0;if (type1[1].equals("S")){s=1;};
                        for (int i=s ;i< value1.length-1;i++){
                            if(!typeList1.get(i+1).equals("I") && !typeList1.get(i+1).equals("P")){
                                sum+=Integer.parseInt(value1[i]);
                            }
                        }
                        int pos1 = Integer.valueOf(SA.split(",")[1])+sum-1;
                        int pos2 = start;
                        if ((pos2-pos1)<Integer.valueOf(this.length))continue;
                        bw.write(chr+"\t"+pos1+"\t"+pos2+"\t"+name+"\tSUPP");bw.newLine();

                    }
                }
                sr.close();
                bw.flush();bw.close();
            }else{
                File[] fs = new File(this.inputFileSDir).listFiles();
                fs = IOUtils.listFilesEndsWith(fs, ".bam");
                List<File> fList = new ArrayList(Arrays.asList());
                for (int i = 0; i < fs.length; i++) {
                    srS[i] = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(fs[i].toString()));
                    bwS[i] = IOUtils.getTextWriter(new File(fs[i].toString().replace(".bam",".txt.gz")).getAbsolutePath());
                    fList.add(fs[i]);
                }
                ExecutorService pool = Executors.newFixedThreadPool(this.threadNum);
                fList.parallelStream().forEach(p -> {
                    try{
                        int index = fList.indexOf(p);
                        SAMRecordIterator r = srS[index].iterator();
                        String cigar = null;String SA = null;
                        String [] value = null; String [] type = null;
                        String [] value1 = null; String [] type1 = null;
                        String name = null ; String chr = null; int start = 0;
                        while(r.hasNext()) {
                            SAMRecord tem=r.next();
                            cigar=tem.getCigarString();
                            if (tem.getMappingQuality() < Integer.parseInt(this.quality))continue;
                            if (!cigar.contains("D")) continue;
                            value = cigar.replaceAll("[a-zA-Z=]"," ").split(" ");
                            type = cigar.replaceAll("[0-9]+"," ").split(" ");
                            ArrayList typeList = new ArrayList<String>(Arrays.asList(type));
                            ArrayList<Integer> allIndexes =
                                    (ArrayList<Integer>) IntStream.range(1, typeList.size()).boxed()
                                            .filter(i -> typeList.get(i).equals("D"))
                                            .collect(Collectors.toList());
                            boolean bl = false;
                            for (int i = 0; i < allIndexes.size() ; i++){
                                if (Integer.parseInt(value[allIndexes.get(i).intValue()-1])>=Integer.valueOf(this.length)){
                                    bl=true;break;
                                }
                            }
                            if (bl==false)continue;
                            name = tem.getReadName();
                            chr = tem.getContig();
                            start=tem.getStart()-1;
                            for (int i =0 ; i < allIndexes.size() ; i++){
                                int num = Integer.parseInt(value[allIndexes.get(i)-1]);
                                if (num<Integer.parseInt(this.length))continue;
                                int sum=0;StringBuilder sb = new StringBuilder();
                                int s=0;if (type[1].equals("S")){s=1;};
                                for (int j = s; j < allIndexes.get(i)-1; j++){
                                    if(!typeList.get(j+1).equals("I") && !typeList.get(j+1).equals("P")){
                                        sum+=Integer.parseInt(value[j]);
                                    }
                                }
                                int startPos=start+sum; int endPos=startPos+num;
                                sb.append(chr).append("\t").append(startPos).append("\t").append(endPos).append("\t").append(name).append("\tCIGAR");
                                bwS[index].write(sb.toString());bwS[index].newLine();
                            }
                            if (tem.getFlags()>2047){
                                SA = tem.getAttributes().get(0).value.toString();
                                String chrSA = SA.split(",")[0].replaceAll("SA:Z:","");
                                if (!chr.equals(chrSA))continue;
                                String cigarSA = SA.split(",")[3];
                                value1 = cigarSA.replaceAll("[a-zA-Z=]"," ").split(" ");
                                type1 = cigarSA.replaceAll("[0-9]+"," ").split(" ");
                                ArrayList typeList1 = new ArrayList<String>(Arrays.asList(type1));
                                int sum = 0;int s=0;if (type1[1].equals("S")){s=1;};
                                for (int i=s ;i< value1.length-1;i++){
                                    if(!typeList1.get(i+1).equals("I") && !typeList1.get(i+1).equals("P")){
                                        sum+=Integer.parseInt(value1[i]);
                                    }
                                }
                                int pos1 = Integer.valueOf(SA.split(",")[1])+sum-1;
                                int pos2 = start;
                                if ((pos2-pos1)<Integer.valueOf(this.length))continue;
                                bwS[index].write(chr+"\t"+pos1+"\t"+pos2+"\t"+name+"\tSUPP");bwS[index].newLine();

                            }
                        }
                        srS[index].close();
                        bwS[index].flush();bwS[index].close();
                    }
                    catch (Exception e) {
                        e.printStackTrace();
                    }
                });
                pool.shutdown();
                pool.awaitTermination(Long.MAX_VALUE, TimeUnit.MICROSECONDS);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
