package app.genotypeDel;

/**
 * @author xujun on 2022-06-08-3:30 PM
 * @project SVCaller
 */

import pgl.infra.utils.Dyad;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;
import utils.AppUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

public class  GenotypeDel {
    /**
     * File path of taxa and there corresponding bams and depth
     */
    String taxaBamFileS = null;
    /**
     * Current chromosome for depth profiling
     */
    String chromosome = null;
    /**
     * Path of delition library
     */
    String libPath = null;
    /**
     * Path of samtools
     */
    String samPath = null;
    /**
     * Number of threads
     */
    int threadNum = 16;
    String outfileS = null;
    String[] taxa = null;
    String[] references = null;
    HashMap<String, String> taxaBamPathsMap = new HashMap<>();
    HashMap<String, Double> taxaDepthMap = new HashMap<>();
    /**
     * Estimated range of mode
     */
    int maxDepthRange = 200;
    /**
     * Current pipeline step
     */
    int step = 0;
    /**
     * Sampling size for estimating mode
     */
    int samplingSize = 10000;
    /**
     * Window size to profile depth
     */
    int windowSize = 500_000;

    double[][] depth = null;
    ArrayList<Double>[] depthR = null;
    ArrayList<Integer>[] depthIndex = null;

    public GenotypeDel(String parameterFileS) {
        this.parseParameters(parameterFileS);
        this.profileDeletionRegion();
    }

    public void profileDeletionRegion() {
        long startTime = System.currentTimeMillis();
        String mark =null; int length =0;List<String> m = new ArrayList<>();
        int startPos = 0; int endPos = 0;
        StringBuilder sb = new StringBuilder();
        StringBuilder sbIndex = new StringBuilder();
        try{
            BufferedReader brLib = null;
            if (this.libPath.endsWith(".gz")) {
                brLib = IOUtils.getTextGzipReader(this.libPath);
            }
            else {
                brLib = IOUtils.getTextReader(this.libPath);
            }
            BufferedWriter bw = IOUtils.getTextGzipWriter(this.outfileS);
            bw.write("chr\tstart\tend\ttaxa");bw.newLine();
            while((mark = brLib.readLine()) != null){
                m = PStringUtils.fastSplit(mark);
                if (!m.get(0).equals(chromosome))continue;
                startPos = Integer.parseInt(m.get(1)); endPos = Integer.parseInt(m.get(2));
                length=endPos-startPos-1;
                int[][] windows = PArrayUtils.getSubsetsIndicesBySubsetSize(length, windowSize);
                int[][] subIndices = PArrayUtils.getSubsetsIndicesBySubsetSize(taxa.length, this.threadNum);
                try {
                    for (int i = 0; i < windows.length; i++) {
                        String[] commands = this.getSamCommands(windows[i][0]+startPos-1, windows[i][1]+startPos);
                        depth = new double[windows[i][1]-windows[i][0]+1][taxa.length];
                        depthR = new ArrayList[taxa.length];depthIndex = new ArrayList[taxa.length];
                        for (int k = 0 ;k < taxa.length;k++){
                            depthR[k] = new ArrayList<>();
                            depthIndex[k] = new ArrayList<>();
                        }
                        int startIndex = windows[i][0]+startPos;
                        int border = windows[i][1]+startPos;
                        for (int u = 0; u < subIndices.length; u++) {
                            List<Integer> indices = PArrayUtils.getIndexList(subIndices[u][0], subIndices[u][1]);
                            indices.parallelStream().forEach(j -> {
                                int v = 0;int pos = -1;
                                try {
                                    Runtime rt = Runtime.getRuntime();
                                    Process p = rt.exec(commands[j]);
                                    BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
                                    String temp = null;List<String> l = new ArrayList<>();
//                                    int v = 0;int pos = -1;
                                    while ((temp = br.readLine()) != null) {
                                        v = 0;
                                        l = PStringUtils.fastSplit(temp);
                                        for (int k = 0; k < l.size()-2; k+=2) {
                                            v+=Integer.parseInt(l.get(k+2));
                                        }
                                        pos = Integer.parseInt(l.get(1));
                                        if (pos > border) break;
                                        depth[pos-startIndex][j] = v;
                                    }
                                    br.close();
                                    p.waitFor();
                                }
                                catch (Exception e) {
                                    System.out.println(pos+"\t"+startIndex);
                                    e.printStackTrace();
                                }
                            });
                        }
                        double ratio = 0;
                        for (int j = 0; j < length; j++) {
                            for (int k = 0 ; k < taxa.length ; k++){
                                if (depth[j][k]==0){
                                    depthIndex[k].add(j);
                                    continue;
                                }
                                ratio = depth[j][k]/this.taxaDepthMap.get(taxa[k]);
                                if (ratio<0.1){
                                    depthIndex[k].add(j);
                                }
                            }
                        }
                        for (int j = 0 ; j < taxa.length; j++){
                            if (depthIndex[j].isEmpty())continue;
                            sbIndex.setLength(0);
                            sbIndex.append(" "+depthIndex[j].get(0));
                            for(int k = 1 ; k < depthIndex[j].size(); k++){
                                if (depthIndex[j].get(k) == (depthIndex[j].get(k-1) + 1)){
                                    sbIndex.append(" "+depthIndex[j].get(k));
                                }else{
                                    sbIndex.append(","+depthIndex[j].get(k));
                                }
                            }
                            String sub1[] = new String[sbIndex.toString().split(",").length];
                            for (int k =0 ; k<sub1.length ;k++){
                                sub1[k]=sbIndex.toString().split(",")[k];
                            }
                            for (int k = 0 ; k < sub1.length; k++){
                                int num = sub1[k].split(" ").length-1;
                                if(num > 40){
                                    sb.setLength(0);
                                    int start = startIndex+Integer.parseInt(sub1[k].split(" ")[1]);
                                    int end = startIndex+Integer.parseInt(sub1[k].split(" ")[num])+2;
                                    sb.append(this.chromosome+ "\t").append(start + "\t").append(end + "\t"+taxa[j]);
                                    bw.write(sb.toString());bw.newLine();
                                }
                            }

                        }
                    }
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
            brLib.close();
            bw.flush();bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        long endTime   = System.currentTimeMillis();
        NumberFormat formatter = new DecimalFormat("#0.00000");
        System.out.println(formatter.format((endTime - startTime) / 1000d));
    }

    private String[] getSamCommands(int startIndex, int endIndex) {
        String[] commands = new String[taxa.length];
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < commands.length; i++) {
            sb.setLength(0);
            sb.append(this.samPath).append(" depth -Q 20 -r ").append(this.chromosome).append(":").append(startIndex+1).append("-").append(endIndex);
            sb.append(" ").append(this.taxaBamPathsMap.get(taxa[i]));
            commands[i] = sb.toString();
        }
        return commands;
    }

    public void parseParameters(String parameterFileS) {
        Dyad<List<String>, List<String>> d = AppUtils.getParameterList(parameterFileS);
        List<String> pLineList = d.getFirstElement();
        this.taxaBamFileS = pLineList.get(0);
        this.chromosome = pLineList.get(1);
        this.libPath = pLineList.get(2);
        this.samPath = pLineList.get(3);
        this.threadNum = Integer.parseInt(pLineList.get(4));
        this.outfileS = pLineList.get(5);
        try {
            BufferedReader br = IOUtils.getTextReader(this.taxaBamFileS);
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                this.taxaBamPathsMap.put(l.get(0), l.get(1));
                this.taxaDepthMap.put(l.get(0),Double.parseDouble(l.get(2)));
            }
            Set<String> tSet= this.taxaBamPathsMap.keySet();
            this.taxa = tSet.toArray(new String[tSet.size()]);
            Arrays.sort(taxa);
            this.references = new String[taxa.length];
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}

